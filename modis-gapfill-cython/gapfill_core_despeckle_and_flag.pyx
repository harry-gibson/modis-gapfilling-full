import numpy as np
from libc.math cimport fabs, sqrt
cimport cython
cimport openmp
from cython.parallel import prange, parallel


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef setSpeckleFlags (
                   dict DataStacks,
                   dict FlagValues,
                   dict SpiralSearchConfig,
                   dict Limits,
                   dict Margins,
                   float _NDV,
                   float _CORRECTION_OFFSET = 0,
                   ):

    '''
    Cython implementation of the "despeckle" algorithm; also applies land-sea mask and any systematic correction.

    Runs the "despeckle" algorithm to identify outlier cells in the data stack,
    based on the number of standard deviations they fall from the mean, and on the similarity with their
    spatial neighbours.

    Two thresholds are provided specifying a number of standard deviations.
    Cells that differ by more than N1 standard deviations from the local mean are determined to be
    definite outliers and are replaced with nodata. Furthermore, cells that differ by more than N2 standard
    deviations (where N1 > N2) from the local mean are also replaced with nodata iif
    their z nearest neighbours don't on average also do so.

    Locations where these replacements have occurred are recorded with bit flags in the flags object that is
    returned. Also recorded here are the absence of any data, and the absence of land.
    '''

    # data, means, stds should have margins allowing the neighbour search to run
    # flags and landMask should be of the shape that is data shape - margins
    # data and flags will be modified in place ready for passing to A1 routine
    cdef:
        # 2d arrays for the upper and lower validity / possible speckle limits
        float [:,::1] validLower, validUpper, dodgyLower, dodgyUpper
        # 2d array of the zscores for the valid cells in the current day
        float[:,::1] zScoresDay
        int[:,::1] nbrIntCoords
        # loop vars etc
        Py_ssize_t zShape, xShapeTotal, yShapeTotal
        Py_ssize_t xShapeToDespeckle, yShapeToDespeckle
        # mnemonic, yT = total (outer coords), yD is shifted into the output depeckled coords
        Py_ssize_t yT, xT_prv, xD_prv, yD_prv, xDNbr_prv, yDNbr_prv
        float value_prv
        long long speckleCount_Glob, extremeCount_Glob, goodCount_Glob, oceanCount_Glob, clearedSpeckleCount_Glob
        Py_ssize_t nbrIndex_prv
        float nbrZscore_prv, zscoreDiff_prv
        int nbrZscoreCount_prv
        float nbrZscoreTotal_prv

        float _BLANK_ZSCORE

    # unpack input dictionaries to typed variables
    cdef:
        # thresholds / consts
        float hardLowerLimit = Limits["HARD_LOWER_LIMIT"]
        float hardUpperLimit = Limits["HARD_UPPER_LIMIT"]
        float _SPECKLE_ZSCORE_THRESHOLD = Limits["ZSCORE_ACCEPTANCE_STDS"]
        float stDevValidityThreshold = Limits["INVALID_BEYOND_STDS"]
        float speckleDevThreshold = Limits["SPECKLE_BEYOND_STDS"]

        int _SPECKLE_NBR_MIN_THRESHOLD = SpiralSearchConfig["MIN_NBRS_REQUIRED"]
        int _SPECKLE_NBR_MAX_THRESHOLD = SpiralSearchConfig["MAX_NBRS_REQUIRED"]
        int _MAX_NEIGHBOURS_TO_CHECK = SpiralSearchConfig["MAX_NBRS_TO_SEARCH"]

        short _OCEAN_FLAG = FlagValues["OCEAN"]
        short _NEVERDATA_FLAG = FlagValues["FAILURE"]
        short _EXTREME_FLAG = FlagValues["EXTREME"]
        short _SPECKLE_FLAG = FlagValues["SPECKLE"]

        int marginT = Margins["TOP"]
        int marginB = Margins["BOTTOM"]
        int marginL = Margins["LEFT"]
        int marginR = Margins["RIGHT"]

        float[:,:,::1] data = DataStacks["Data"]
        unsigned char[:,:,::1] flags = DataStacks["Flags"]
        float[:,::1] means = DataStacks["Means"]
        float[:,::1] stds = DataStacks["Stds"]
        char[:,::1] landMask = DataStacks["LandMask"]

        float[:,:,::1] outputData

    zShape = data.shape[0]
    yShapeTotal = data.shape[1]
    xShapeTotal = data.shape[2]
    # we will only iterate thru cells that are not in the margins
    yShapeToDespeckle = data.shape[1] - (marginT + marginB)
    xShapeToDespeckle = data.shape[2] - (marginL + marginR)
    # and we will only set flags within the cells that A1 cares about
    #yShapeActual = data.shape[1] - (marginTTot + marginBTot)
    #xShapeActual = data.shape[2] - (marginLTot + marginRTot)

    assert means.shape[0] == stds.shape[0] == yShapeTotal
    assert means.shape[1] == stds.shape[1] == xShapeTotal
    assert landMask.shape[0] == yShapeTotal
    assert landMask.shape[1] == xShapeTotal

    assert flags.shape[1] == yShapeTotal
    assert flags.shape[2] == xShapeTotal

    validLower = np.empty(shape=(yShapeTotal, xShapeTotal), dtype='float32', order='c')
    validUpper = np.empty(shape=(yShapeTotal, xShapeTotal), dtype='float32', order='c')
    dodgyLower = np.empty(shape=(yShapeTotal, xShapeTotal), dtype='float32', order='c')
    dodgyUpper = np.empty(shape=(yShapeTotal, xShapeTotal), dtype='float32', order='c')
    zScoresDay = np.empty(shape=(yShapeTotal, xShapeTotal), dtype='float32', order='c')

    _BLANK_ZSCORE = -999.0
    speckleCount_Glob = 0
    extremeCount_Glob = 0
    goodCount_Glob = 0
    oceanCount_Glob = 0
    clearedSpeckleCount_Glob = 0

    print ("Despeckle: Rejecting data beyond {0!s}s.d. of mean. Nbr search on data beyond {1!s} s.d. of mean.".
           format(stDevValidityThreshold, speckleDevThreshold))
    print ("Nbr searching for {0!s} - {1!s} nbrs within {2!s} spiral steps for z-score tolerance of {3!s}".
           format( _SPECKLE_NBR_MIN_THRESHOLD, _SPECKLE_NBR_MAX_THRESHOLD, _MAX_NEIGHBOURS_TO_CHECK,
                  _SPECKLE_ZSCORE_THRESHOLD))

    # Generate the neighbour spiral search table out to "a bit" further than needed
    _SEARCH_RADIUS =  <int> ((sqrt(_MAX_NEIGHBOURS_TO_CHECK / 3.14)) + 5)
    diam = _SEARCH_RADIUS * 2 + 1
    print "Despeckle diam = " + str(diam)
    inds = np.indices([diam, diam]) - _SEARCH_RADIUS
    distTmp = np.sqrt((inds ** 2).sum(0))
    npTmpTable = ((inds.T).reshape(diam ** 2, 2))
    npTmpTable = np.append(npTmpTable, distTmp.ravel()[:, None], axis=1)

    # sort the table by distance then x then y (the arguments are last-sort-first)
    order = np.lexsort((npTmpTable[:, 1],npTmpTable[:, 0],npTmpTable[:, 2]))
    npTmpTable = np.take(npTmpTable, order, axis=0)

    # transfer to a C-side object transposed to have three rows and many columns and in
    # C-contiguous layout, so that cython can access individual nbr coord sets more quickly
    nbrTable = np.copy((npTmpTable[npTmpTable[:,2] <= _SEARCH_RADIUS]).T, order='c')
    # cast the columns that will be used as array indices to int type once here, rather
    # than casting repeatedly inside the inner loop
    nbrIntCoords = np.asarray(nbrTable[0:2,:]).astype(np.int32)

    # We can't modify the input in this algorithm as we have parallelised it
    # and we look at adjacent results in setting the cell's value - we only want to use
    # original data for this check, so we have to keep a pristine copy
    outputData = np.empty_like(data, order='c')
    outputData[:] = _NDV

    # We need to run speckle check IN the margins so the data in those margins can be eliminated if necessary
    # and not used in checking other values not in the margins
    # However the speckle search radius is smaller than the A1 radius
    # So we give speckle check a smaller margin so it can be happy but everything that A1 needs
    # has been despeckled

    # Apply any correction offset (to correct messed-up celsius-kelvin conversion!!)
    if _CORRECTION_OFFSET != 0:
        for z in range (zShape):
            with nogil, parallel(num_threads=20):
                for yT in prange (yShapeTotal, schedule='static', chunksize=200):
                    xT_prv = -1
                    for xT_prv in range (xShapeTotal):
                        if data[z, yT, xT_prv] != _NDV:
                            data[z, yT, xT_prv] = data[z, yT, xT_prv] + _CORRECTION_OFFSET

    # Precalculate the thresholds. These are the same for every day in the stack so spend
    # some memory to avoid doing it 15 times (this is the only economy in calling this routine
    # with a 3d array, so if memory is an issue just modify to run on a 2d array one day at a time)
    # run for total shape i.e. including speckle margins
    with nogil, parallel(num_threads=20):
        for yT in prange (yShapeTotal, schedule='static', chunksize=200):
            xT_prv = -1
            for xT_prv in range (xShapeTotal):
                if landMask[yT,xT_prv] == 0:
                    continue
                if means[yT, xT_prv] == _NDV:
                    continue
                if _CORRECTION_OFFSET != 0:
                    means[yT, xT_prv] = means[yT, xT_prv] + _CORRECTION_OFFSET
                validLower[yT,xT_prv] = means[yT,xT_prv] - stds[yT,xT_prv] * stDevValidityThreshold
                validUpper[yT,xT_prv] = means[yT,xT_prv] + stds[yT,xT_prv] * stDevValidityThreshold
                dodgyLower[yT,xT_prv] = means[yT,xT_prv] - stds[yT,xT_prv] * speckleDevThreshold
                dodgyUpper[yT,xT_prv] = means[yT,xT_prv] + stds[yT,xT_prv] * speckleDevThreshold

    # now identify cells which _may_ be speckles by comparison to the limits just calculated
    # also apply land-sea mask and calculate z scores in this pass
    # Run on full dataset incl margins
    for z in range (zShape):
        #re initialise zscoresday
        zScoresDay[:] = _BLANK_ZSCORE
        with nogil,parallel(num_threads=20):
            # do for everything, including speckle margins
            for yT in prange (yShapeTotal, schedule='static', chunksize=200):
                # assign to the inner loop's variables so that Cython will make them thread-private
                xT_prv = -1
                value_prv = _NDV
                nbrIndex_prv = -1
                nbrZscoreCount_prv = 0
                nbrZscoreTotal_prv = 0.0
                nbrZscore_prv = 0.0
                zscoreDiff_prv = 0.0

                for xT_prv in range(xShapeTotal):
                    if landMask[yT, xT_prv] == 0:
                        # do not do anything where it is not in the land
                        # just set the flags image to record this
                        flags[z, yT, xT_prv] = flags[z, yT, xT_prv] | _OCEAN_FLAG

                        # Do not copy across any data existing in the sea. It's likely to be squiffy
                        # (EVI in particular) and if we leave it, the fill process inland could use these
                        # dodgy values
                        # (Output data is already NDV from initialisation)

                        # Tracking
                        oceanCount_Glob += 1
                        continue

                    # if no mean then it means we never have data in the stack and will not be able to
                    # produce a fill
                    if means[yT, xT_prv] == _NDV:
                        flags[z, yT, xT_prv] = flags[z, yT, xT_prv] | _NEVERDATA_FLAG
                        continue

                    value_prv = data[z, yT, xT_prv]
                    # is the value outside extreme thresholds? delete if so
                    if (value_prv < hardLowerLimit or value_prv > hardUpperLimit or
                        value_prv < validLower[yT, xT_prv] or value_prv > validUpper[yT, xT_prv]
                        or value_prv == _NDV):
                        if value_prv != _NDV:
                            flags[z, yT, xT_prv] = flags[z, yT, xT_prv] | _EXTREME_FLAG
                            extremeCount_Glob += 1
                        # no flag to set if it IS nodata because the fact it is nodata is all the flag
                        # we need...
                        continue

                    # implicit else
                    # pre-calculate the z-scores for all cells that are not definitely invalid
                    zScoresDay[yT, xT_prv] = (value_prv - means[yT, xT_prv]) / stds[yT, xT_prv]
                    outputData [z, yT, xT_prv] = value_prv
                    # see if the location is a possible speckle
                    if value_prv >= dodgyLower[yT, xT_prv] and value_prv <= dodgyUpper[yT, xT_prv]:
                        # it's within the local limits, so it counts as good data, so pass through
                        goodCount_Glob += 1
                    else:
                        # it's outside the speckle thresholds, but within the validity thresholds
                        # so set the speckle flag (we may remove it again later after the nbr check)
                        flags[z, yT, xT_prv] = flags[z, yT, xT_prv] | _SPECKLE_FLAG

        # now go through again, the z-scores are all populated wherever possible so no need
        # to re-calculate every time a nbr is checked as in orginal implementation
        # do not run for speckle margins - use inner shape - but this will still include the A1 margins
        # if everything's been set up right!
        with nogil,parallel(num_threads=20):
            for yT in prange (yShapeToDespeckle, schedule='dynamic', chunksize=200):

                xT_prv = -1
                xD_prv = -1
                yD_prv = yT + marginT
                xDNbr_prv = -1
                yDNbr_prv = -1

                nbrIndex_prv = -1
                nbrZscoreCount_prv = 0
                nbrZscoreTotal_prv = 0.0
                nbrZscore_prv = 0.0
                zscoreDiff_prv = 0.0

                for xT_prv in range(xShapeToDespeckle):
                    xD_prv = xT_prv + marginL
                    nbrIndex_prv = -1
                    nbrZscoreCount_prv = 0
                    nbrZscoreTotal_prv = 0.0
                    nbrZscore_prv = 0.0
                    zscoreDiff_prv = 0.0
                    if (flags[z, yD_prv, xD_prv] & _SPECKLE_FLAG) != _SPECKLE_FLAG:
                        continue

                    # else we are on a possible speckle
                    nbrZscoreCount_prv = 0
                    nbrZscoreTotal_prv = 0.0
                    # spiral search of the nbrs to see if they too are far from
                    # the mean in which case this isn't a speckle after all
                    for nbrIndex_prv in range(1, _MAX_NEIGHBOURS_TO_CHECK + 1):
                        # This loop may run 100s or 1000s of times for every potential
                        # speckle cell, i.e. potentially trillions of times.
                        # So it's really important to get this loop fast.
                        if nbrZscoreCount_prv == _SPECKLE_NBR_MAX_THRESHOLD:
                            # we've found sufficient neighbours so stop looking
                            break
                        # use int-type coords array to avoid cast op in tight loop
                        xDNbr_prv = xD_prv + nbrIntCoords[0, nbrIndex_prv]
                        yDNbr_prv = yD_prv + nbrIntCoords[1, nbrIndex_prv]

                        # is the requested neighbour cell within data bounds and with data?
                        if (xDNbr_prv >= 0 and xDNbr_prv < xShapeTotal and
                                yDNbr_prv >= 0 and yDNbr_prv < yShapeTotal and
                                zScoresDay[yDNbr_prv, xDNbr_prv] != _BLANK_ZSCORE):
                            nbrZscoreCount_prv = nbrZscoreCount_prv + 1 # tracking
                            nbrZscoreTotal_prv = nbrZscoreTotal_prv + zScoresDay[yDNbr_prv, xDNbr_prv]

                    # did we find enough neighbours?
                    if nbrZscoreCount_prv >= _SPECKLE_NBR_MIN_THRESHOLD:
                        # were the neighbours similar to this cell in terms of distance from mean?
                        nbrZscore_prv = nbrZscoreTotal_prv / nbrZscoreCount_prv
                        zscoreDiff_prv = fabs(nbrZscore_prv - zScoresDay[yD_prv, xD_prv])
                        if zscoreDiff_prv < _SPECKLE_ZSCORE_THRESHOLD:
                            # this pixels neighbours are also far from the mean, it's probably
                            # a volcano or a country that caught fire or something, so we'll believe it,
                            # clear the flag, and carry through the original value
                            flags[z, yD_prv, xD_prv] = flags[z, yD_prv, xD_prv] ^ _SPECKLE_FLAG
                            clearedSpeckleCount_Glob += 1
                            continue
                    # insufficient neighbours or the neighbours are less extreme.
                    # so this is a speckle; delete the data
                    speckleCount_Glob += 1
                    outputData[z, yD_prv, xD_prv] = _NDV

    print "Speckle count:    " + str(speckleCount_Glob)
    print "Extreme count:    " + str(extremeCount_Glob)
    print "Good count:       " + str(goodCount_Glob)
    print "Cleared Speckle:  " + str(clearedSpeckleCount_Glob)
    print "Ocean count:      " + str(oceanCount_Glob)
    return (outputData, flags)