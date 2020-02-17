cimport cython
import numpy as np
from libc.stdlib cimport abs
from cython.parallel import prange, parallel
from libc.math cimport sqrt
from gapfill_config import FlagItems, DataCharacteristicsConfig, A1SearchConfig, A1Diagnostics
from gapfill_utils import A1DataStack, PixelMargins


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
# pass in the whole data that should be checked, and optionally the margins that are provided
# - we will only fill gaps within the portion of the data that isn't excluded in the margin,
# but we will use the margin data to provide neighbour data for those pixels
# So for a tiled run data can (should) be provided with a padding of _SEARCH_RADIUS and the
# margin set to this value also. That way the algorithm can get values for "edge pixels" of its
# input data
cpdef a1_core(A1DataStack dataStacks,  # has items Data, Flags, DistTemplate (optional), KnownUnfillable (optional)
              PixelMargins margins,
              FlagItems flagValues,
              DataCharacteristicsConfig dataConfig,
              A1SearchConfig a1Config,
              Py_ssize_t nCores
              ):

    """
    Optimised, multithreaded Cython (C) implementation of the A1 gapfilling algorithm.

    a1_core(A1DataStack DataStacks, PixelMargins Margins, FlagItems flagValues,
            DataCharacteristicsConfig dataConfig, A1SearchConfig a1Config, 
            char RunFillFromPos = 0
            Py_ssize_t nCores)
                    -> A1Diagnostics
                    (members of dataStacks are modified in-place)

    Runs A1 on a stack of data of a given calendar day; the stack should have one layer (0th dimension)
    for each year. The algorithm can be run in a tiled or sliced fashion by passing in a stack
    that isn't the full extent. In this case the stack passed should have "surplus" data at the edges
    so that the neighbourhood search can proceed. The size of these edges should be specified with the
    margin parameters. These should be at least as big as the neighbourhood search radius which is set
    to 3142 cells area = ~31km radius.
    The flags should be the same shape as the data EXCLUDING the margins.
    Returns an object containing diagnostics from the fill process - the data, distance, and flags images within 
    the A1DataStack are modified in-place.
    """
    cdef:
        # define all arrays (cython memoryview objects) as being c-contiguous
        # so cython doesn't treat as strided (slower access)
        # this gives major speedup in access in tight loops as the cpu caching can
        # work more effectively

        # intermediate arrays
        double [:,::1] nbrTable
        int[:,::1] nbrIntCoords
        char [:,::1] neverDataLocs
        char  noTemplate
        # OUTPUTS
        float [:,:,::1] outputData
        float [:,:,::1] outputDists
        unsigned char [:,:,::1] outputFlags
        # in-loop working vars
        Py_ssize_t z, y
        Py_ssize_t zShape, yShape, xShape
        Py_ssize_t yShapeTotal, xShapeTotal
        int [::1] deltas
        double posInf,negInf


        #thread-private ones:
        Py_ssize_t x_prv
        Py_ssize_t newZ_prv, xi_prv, yi_prv, xNbr_prv, yNbr_prv
        int spiralStart_prv
        Py_ssize_t max_idx_prv, min_idx_prv
        double ws_prv, sw_prv,  pfv_prv, weight_prv, wfv_prv
        double altValue_prv, currentValue_prv
        int nfound_prv
        int delta_prv, deltaidx_prv, nbrIndex_prv

        # intermediate values for ratio-based calculations
        # (now just using the same vars whichever approach is used)
        #double max_Ratio_prv, maxR_Dist_prv, maxR_wpfv_prv, maxR_weight_prv
        #double min_Ratio_prv, minR_Dist_prv, minR_wpfv_prv,  minR_weight_prv
        #double ratio_prv
        # intermediate values for offset (addition) based calculations
        double max_Diff_prv, maxD_Dist_prv, maxD_wpfv_prv, maxD_weight_prv
        double min_Diff_prv, minD_Dist_prv, minD_wpfv_prv, minD_weight_prv
        double valueDiff_prv

        double sumDist_prv
        unsigned char flag_prv

        # First loop stuff
        Py_ssize_t xi, yi, x
        int nfound

        # Unpack the parameter objects to typed variables
        # (the dictionaries were just to clean up the signature a bit in the absence of motivation
        # to make some specific type to pass the data more cleanly)
        char _OCEAN_FLAG = flagValues.OCEAN
        char _FAILURE_FLAG = flagValues.FAILURE
        char _SUCCESS_FLAG = flagValues.A1_FILLED
        char _SUCCESS_WAS_FULL_FLAG = flagValues.A1_FULL
        int marginT = margins.Top
        int marginB = margins.Bottom
        int marginL = margins.Left
        int marginR = margins.Right
        float _NDV = dataConfig.NODATA_VALUE
        unsigned char RunFillFromPos = dataStacks.FillFromZPosition
        unsigned char FillByRatios = A1SearchConfig.FILL_GENERATION_METHOD ==" RATIO"
        unsigned char _TRIM_MIN_MAX = A1SearchConfig.TRIM_FULL_OUTLIERS

        float [:,:,::1] dayDataStack = dataStacks.DataArray3D
        unsigned char[:,:,::1] inputFlags = dataStacks.FlagsArray3D # embeds the land-sea mask (sea=1 land =0)
        unsigned char[:,:,::1] dataDistTemplate = None
        unsigned char[:,::1] knownUnfillableLocs = None
        float _AbsZeroPoint = dataConfig.ABSOLUTE_ZERO_OFFSET
        float _MaxAllowableRatio = A1SearchConfig.MAX_ALLOWABLE_RATIO
        float _MinAllowableRatio = 1.0 / _MaxAllowableRatio

        #  how many locations should be checked in spiral search (gives radius), historic defailt is 3142
        int _MAX_NEIGHBOURS_TO_CHECK = A1SearchConfig.MAX_NBRS_TO_SEARCH

        # Only use the values gleaned from up to this number of cells (Even if more are avail within radius) 640
        int _FILL_THRESHOLD = A1SearchConfig.MAX_USED_NBRS

        # min number of values that must be found to have a valid fill, historic default is 320
        int _FILL_MIN = A1SearchConfig.MIN_REQUIRED_NBRS

        # calc the distance that is implied by the max spiral search length
        int _SEARCH_RADIUS = <int> (sqrt((_MAX_NEIGHBOURS_TO_CHECK*2.0) / 3.14))  + 1

        # calculation / tracking vars: will be reduction variables - incremented but not read by threads
        # for the whole globe there are 933,120,000 pixels per image so this exceeds
        # max val of uint after 4 images!
        long long totalProcessedGapCells,  totalCells, oceanCells, neverDataCells
        long long scannedLevels, scannedNeighbours, usedNeighbours
        long long filledBelowThreshold, noPairsFound, insufficientPairsFound, filledToThreshold
        long long gapsAtUnfillableLocs, gapsTooBig, dataGood



    # can precalc dist from gap to nearest data for more  efficient spiral search
    # (no need to start spiral search closer than the known nearest data pixel)
    if dataStacks.DistanceTemplate3D is not None:
        dataDistTemplate = dataStacks.DistanceTemplate3D
    # can precalc locations where no alternate years exist (so no fill will be possible)
    if dataStacks.KnownUnfillable2D is not None:
        knownUnfillableLocs = dataStacks.KnownUnfillable2D

    zShape = dayDataStack.shape[0]
    # size of the input data including search margins
    yShapeTotal = dayDataStack.shape[1]
    xShapeTotal = dayDataStack.shape[2]
    # size of the data that needs to be filled.
    # we will only iterate thru cells that are not in the margins
    # (they are equal (zero margin) if filling a global image in one go, and the margin is also zero
    # for the slice edges at the edge of the global images)
    yShape = dayDataStack.shape[1] - (marginT + marginB)
    xShape = dayDataStack.shape[2] - (marginL + marginR)

    # cython doesn't have inf defined
    posInf = np.inf
    negInf = -np.inf

    # generate nbr distance table (offset coordinates for the spiral search steps) using numpy
    diam = _SEARCH_RADIUS * 2 + 1
    inds = np.indices([diam,diam]) - _SEARCH_RADIUS
    distTmp = np.sqrt((inds ** 2).sum(0))
    npTmpTable = ((inds.T).reshape(diam**2, 2))
    npTmpTable = np.append(npTmpTable, distTmp.ravel()[:,None],axis=1)
    # sort the table by distance then x then y (the arguments are last-sort-first)
    order = np.lexsort((npTmpTable[:,1], npTmpTable[:,0], npTmpTable[:,2]))
    npTmpTable = np.take(npTmpTable, order, axis=0)

    # the C-side result of the distance table calculations
    # transpose it to have three rows and many columns and take a C contiguous copy of this
    # so that cython can access individual nbr coord sets more quickly
    nbrTable = np.copy((npTmpTable[npTmpTable[:,2] <= _SEARCH_RADIUS]).T,order='c')

    # the distance table is stored with a float type (for the distances) but we need ints
    # for indexing based on coords. we can cast at the time we get them out, but as this happens
    # in the innermost loop it is done trillions of times and so the time penalty of that
    # becomes signficant. So, store an int version of the coords array
    nbrIntCoords = np.asarray(nbrTable[0:2,:]).astype(np.int32)

    # initialise C arrays each from an appropriate empty numpy object
    # these will have the same size as the area we're actually filling i.e. not the (buffered)
    # input section
    outputData = np.empty(shape=(zShape, yShape, xShape), dtype='float32', order='c')
    outputDists = np.empty(shape=(zShape, yShape, xShape), dtype='float32', order='c')
    outputFlags = np.zeros(shape=(zShape, yShape, xShape), dtype='uint8', order='c')
    outputData[:] = _NDV
    outputDists[:] = _NDV

    # handle optional parameters
    if dataDistTemplate is None:
        noTemplate = 1
    else:
        noTemplate = 0
        assert (dataDistTemplate.shape[0] == zShape and
                dataDistTemplate.shape[1] == yShape and
                dataDistTemplate.shape[2] == xShape
                )

    # diagnostics; TODO use logging
    print ("Running A1 (Full Spiral Search).")
    print ("No data template: {0!s}. Using fill generation method: {1!s}. Searching for {2!s} - {3!s} nbrs within {4!s} spiral steps".
           format(noTemplate, A1SearchConfig.FILL_GENERATION_METHOD, _FILL_MIN, _FILL_THRESHOLD, _MAX_NEIGHBOURS_TO_CHECK))
    print ("Calculating nbr table out to radius of {0!s}.".format(_SEARCH_RADIUS))
    print ("Filling from stack position {0!s}.".format(RunFillFromPos))

    assert inputFlags.shape[1] == yShapeTotal
    assert inputFlags.shape[2] == xShapeTotal

    if knownUnfillableLocs is None: # i.e. has not been precalced with np.all
        # knownUnfillableLocs will a 2D map of cells where there isk no data in any year (z axis)
        knownUnfillableLocs = np.ones(shape=(yShape, xShape), dtype=np.uint8, order='c')
        # locate cells that CAN be filled (negating it like this means we can keep
        # the loop order cache-friendly - as opposed to iterating through Z at each
        # y,x location until/unless we find a value > 0.
        # This process is still slower on 1 thread than numpy.all, but not significant overall,
        # and multithreading makes it way faster to do here
        with nogil:
            for z in range(zShape):
                # work is trivial and fairly well balanced so use static schedule
                for y in prange(yShape, schedule='static', num_threads=nCores):
                    x_prv = -1
                    yi_prv = y + marginT
                    for x_prv in range(xShape):
                        xi_prv = x_prv + marginL
                        if dayDataStack[z, yi_prv, xi_prv] != _NDV:
                            knownUnfillableLocs[y, x_prv] = 0
    else:
        assert (knownUnfillableLocs.shape[0] == yShape and
                knownUnfillableLocs.shape[1] == xShape
                )

    if FillByRatios:
        # we need to work with ratios between cells but will be running this on non-ratio scale variables
        # (i.e. where 0 is an arbitrary point) such as temp in celsius. This could lead to div/0 giving infinite
        # ratios, and it's not valid anyway (e.g. 20/10 != 10/-10 even though the interval is the same).
        # The latter point doesn't matter _much_ because we only use the ratio to multiply back against an alternate
        # value which will be similar. Use a large "absolute zero" relative to the values and it matters even less
        # (but not too large, to avoid FP errors)
        with nogil:
            for z in range(zShape):
                for y in prange(yShapeTotal, schedule='static', num_threads=nCores):
                    x_prv = -1
                    for x_prv in range(xShapeTotal):
                        if dayDataStack[z, y, x_prv] != _NDV:
                            dayDataStack[z, y, x_prv] = dayDataStack[z, y, x_prv] -  _AbsZeroPoint
    else:
        _AbsZeroPoint = 0 # to avoid need for further check in loop below

    # initialise metrics
    totalProcessedGapCells = 0 # for testing
    totalCells = 0 # for testing
    scannedLevels = 0 # for testing
    scannedNeighbours = 0 # for testing
    filledToThreshold = 0
    filledBelowThreshold=0
    noPairsFound = 0
    insufficientPairsFound = 0

    gapsAtUnfillableLocs = 0
    gapsTooBig=0
    usedNeighbours = 0
    gapsInKnownUnfillable = 0
    dataGood = 0
    oceanCells = 0
    neverDataCells = 0

    # a single call to produce the alternating forward/back year offsets
    deltas = alternates_cy(zShape)

    # now do A1 gap fill! Iterate through the entire data stack
    # The arrays are C-ordered so it's fastest in terms of CPU cache to have
    # X on the innermost loop. Parallelise over the y axis.

    # Telling cython how to parallelise the variables is done implicitly by how you access them.
    # - assign to a variable if you want it to be thread-private i.e. x = x + 1
    # - increase it in-place if you want it to be shared i.e. x += 1
    # The latter case turns it into a "reduction" variable which cannot be read in the parallel
    # loop.
    # Therefore all variables that need to be thread private must be made so by artifially assigning
    # (something) to them within the parallel block.
    # Check the generated C code to be sure it's worked as intended!
    with nogil, parallel(num_threads=nCores): # change num cores here
        for z in range(zShape):
            if z >= RunFillFromPos:
                # experiments with chunksizes 500, 50, 20, 2, and parallelising on z axis instead,
                # showed 20 to be the quickest (tradeoff between allocation overhead vs one thread doing
                # a "gappy" area all by itself). It's likely to vary by machine and num threads due to
                # CPU cache differences i.e. the spiral search needs a lot of cache to be efficient
                # as it crosses rows so the access cannot be completely cache-economical.
                # However 6 threads was still slower than 12 on a 6-physical 12-virtual machine
                for y in prange(yShape, schedule='dynamic', chunksize=20):
                    # assign something to all variables that are used in filling one cell and need
                    # to be private (i.e. they get modified but we don't want the iteration that fills
                    # another cell to know about that).
                    # It's verbose but this tricks cython into making them thread-private
                    nfound_prv = 0
                    ws_prv = 0
                    sw_prv = 0
                    weight_prv=0
                    pfv_prv=0
                    #ratio_prv=0
                    valueDiff_prv = 0
                    delta_prv=0
                    deltaidx_prv = -1
                    currentValue_prv=0
                    altValue_prv=0
                    xNbr_prv=-1
                    yNbr_prv=-1
                    xi_prv=-1
                    yi_prv=-1
                    x_prv=-1
                    newZ_prv=-1
                    nbrIndex_prv=-1
                    sumDist_prv=0

                    max_Diff_prv = negInf
                    maxD_Dist_prv = 0
                    maxD_weight_prv=0
                    maxD_wpfv_prv=0

                    min_Diff_prv = posInf
                    minD_Dist_prv = 0
                    minD_wpfv_prv = 0
                    minD_weight_prv = 0

                    max_idx_prv = 0
                    min_idx_prv = 0
                    flag_prv = -1

                    yi_prv = y + marginT

                    for x_prv in range(xShape):
                        xi_prv = x_prv + marginL
                        #re initialise variables for this cell
                        nfound_prv= 0
                        ws_prv = 0.0
                        sw_prv = 0
                        weight_prv = 0
                        pfv_prv = 0
                        #ratio_prv = 0
                        valueDiff_prv = 0
                        delta_prv = 0
                        deltaidx_prv = -1
                        currentValue_prv = 0
                        altValue_prv = 0
                        xNbr_prv = -1
                        yNbr_prv = -1
                        newZ_prv = -1
                        nbrIndex_prv = -1
                        sumDist_prv = 0

                        max_Diff_prv = negInf
                        maxD_Dist_prv = 0
                        maxD_weight_prv = 0
                        maxD_wpfv_prv = 0

                        min_Diff_prv = posInf
                        minD_Dist_prv = 0
                        minD_wpfv_prv = 0
                        minD_weight_prv = 0

                        max_idx_prv = 0
                        min_idx_prv = 0
                        flag_prv = -1

                        totalCells += 1 # testing log

                        if ((inputFlags[z, yi_prv, xi_prv] & _OCEAN_FLAG) == _OCEAN_FLAG):
                            oceanCells += 1
                            outputFlags[z, y, x_prv] = inputFlags[z, yi_prv, xi_prv]
                            #in the ocean do not copy across (even if there is data; MODIS may not match shore cleanly)
                            #output is already nodata and flag is alraady set so just
                            continue

                        if ((inputFlags[z, yi_prv, xi_prv] & _FAILURE_FLAG) == _FAILURE_FLAG):
                            # also do this if failure flag (2) is set (by despeckle algorithm indicating that mean
                            # was ND thus never any data on any calendar day)
                            outputFlags[z, y, x_prv] = inputFlags[z, yi_prv, xi_prv]
                            neverDataCells +=1
                            continue

                        currentValue_prv = dayDataStack[z, yi_prv, xi_prv]

                        #if value is valid just copy it across
                        if not currentValue_prv == _NDV: # 0.0 is valid!!
                            # we have data, copy it across - output does not have margins
                            outputData[z, y, x_prv] = currentValue_prv + _AbsZeroPoint
                            outputFlags[z, y, x_prv] = inputFlags[z, yi_prv, xi_prv]
                            dataGood += 1
                            # flag is already non-ocean so just
                            continue

                        # if it's a location that has no data in any alternate year of this calendar day
                        # then we also cannot fill
                        if knownUnfillableLocs[y, x_prv] == 1:
                            totalProcessedGapCells += 1
                            gapsAtUnfillableLocs += 1
                            # for testing to check failure reason set this to 3?
                            outputFlags[z, y, x_prv] = outputFlags[z, y, x_prv] | _FAILURE_FLAG
                            # leave output at _NDV and just
                            continue

                        # if we have a precalculated "large gaps" template does this show that
                        # the cell will be unfillable due to lack of neighbours?
                        # The template gives locations that are more than _SEARCH_RADIUS
                        # away from any data point in the same year, i.e. where the spiral search will not
                        # be able to find anything. Precalculating this in scipy is generally somewhat faster
                        # than leaving the code to iterate through everything here IF we only use 1 core
                        # NB sc evaluation makes this safe whether or not dataDistTemplate exists
                        if (noTemplate == 0 and dataDistTemplate[z, yi_prv, xi_prv] > _SEARCH_RADIUS + 1):
                            # TODO set a flag for this
                            outputFlags[z, y, x_prv] = outputFlags[z, y, x_prv] | _FAILURE_FLAG
                            totalProcessedGapCells += 1
                            gapsTooBig += 1
                            continue

                        # else attempt to fill the gap, using A1
                        if 1: # can't be bothered to rejig indents
                            # METRICS
                            totalProcessedGapCells += 1
                            flag_prv = inputFlags[z, yi_prv, xi_prv]
                            # deltas will be (1, -1, 2, -2, ...)
                            for deltaidx_prv in range(0, deltas.shape[0]):
                                delta_prv = deltas[deltaidx_prv]
                                newZ_prv = z + delta_prv
                                # does this zDelta fall off the start or end of the stack?
                                # or is the alternate year also blank here?
                                # (or, obviously, are we checking the current year!)
                                if (newZ_prv >= zShape or newZ_prv < 0
                                    or newZ_prv == z
                                    or dayDataStack[newZ_prv, yi_prv, xi_prv] == _NDV):
                                    continue
                                # otherwise...
                                scannedLevels += 1 # tracking year-switches

                                # note that the alternate value can be zero. on a ratio based fill
                                # this would result in a fill value of zero, too, no matter what the nbr values are
                                altValue_prv = dayDataStack[newZ_prv, yi_prv, xi_prv]

                                #step thru spiral-outward coords table
                                spiralStart_prv = 1 # zero is this cell
                                # if possible, check where to start the search based on known closest data
                                if noTemplate == 0:
                                    spiralStart_prv = <int>(((dataDistTemplate[z, yi_prv, xi_prv]-1)**2) * 3.141)
                                    if spiralStart_prv < 1:
                                        spiralStart_prv = 1
                                    elif spiralStart_prv > _MAX_NEIGHBOURS_TO_CHECK:
                                        spiralStart_prv = _MAX_NEIGHBOURS_TO_CHECK

                                # spiral search for this alternate year
                                for nbrIndex_prv in range(spiralStart_prv, _MAX_NEIGHBOURS_TO_CHECK + 1):
                                    #scannedNeighbours += 1 # don't track as it's unnecessary op in inner loop

                                    # This is the inner loop that runs once for every neighbour of
                                    # every gap cell, i.e. potentially many trillions of times.
                                    # So it's really important to get this loop fast.

                                    # Total time of this loop where the if below fails is now ~4.7nS
                                    # or around 12 CPU cycles. Considering there are potentially 4 array
                                    # accesses and several comparisons this is probably as good as we can
                                    # get and it's clear that the L1D CPU cache is mostly being hit (as a
                                    # single access to L2 cache is ~14 cycles)

                                    # moving int cast outside the loop saves
                                    # ~ 2.7nS per loop i.e. >50% on an invalid iteration
                                    xNbr_prv = xi_prv + nbrIntCoords[0, nbrIndex_prv]
                                    yNbr_prv = yi_prv + nbrIntCoords[1, nbrIndex_prv]

                                    # is the requested neighbour cell within data bounds?
                                    # (allowing for the fact that the input data can have greater
                                    # extent than the cells we are searching to allow tiled processing)
                                    # and does it have data in both relevant years?
                                    # Put current year check before alternate year check as current year is more
                                    # likely to not have data (given we are already at a gap) so it's more
                                    # efficient to bail out there first
                                    if (xNbr_prv >= 0 and xNbr_prv < xShapeTotal and
                                            yNbr_prv >= 0 and yNbr_prv < yShapeTotal and
                                            dayDataStack[z, yNbr_prv, xNbr_prv] != _NDV and
                                            dayDataStack[newZ_prv, yNbr_prv, xNbr_prv] != _NDV):
                                        # we're good to go! do the maths for this contributing cell
                                        usedNeighbours +=  1 # tracking

                                        if FillByRatios == 0:
                                        #calculate difference not ratio
                                            valueDiff_prv = (dayDataStack[z, yNbr_prv, xNbr_prv]
                                                             - dayDataStack[newZ_prv, yNbr_prv, xNbr_prv])
                                            pfv_prv = altValue_prv + valueDiff_prv
                                        else:
                                            # calculate ratio, taking some precautions to avoid stupid values
                                            # when one is close to zero
                                            valueDiff_prv = (dayDataStack[z, yNbr_prv, xNbr_prv]
                                                             / dayDataStack[newZ_prv, yNbr_prv, xNbr_prv])
                                            if (dayDataStack[newZ_prv, yNbr_prv, xNbr_prv] == 0
                                                or valueDiff_prv > _MaxAllowableRatio):
                                                valueDiff_prv = _MaxAllowableRatio
                                            elif valueDiff_prv < _MinAllowableRatio:
                                                valueDiff_prv = _MinAllowableRatio
                                            pfv_prv = altValue_prv * valueDiff_prv + _AbsZeroPoint
                                        # implicit assumption that one year and one cell distance have the same weighting
                                        weight_prv = (1.0 / abs(delta_prv)) * (1.0 / nbrTable[2,nbrIndex_prv])
                                        ws_prv += pfv_prv * weight_prv
                                        sw_prv += weight_prv
                                        # the only thing the replacement values array got used for in original IDL was
                                        # calculating the mean contibuting distance. we don't need to
                                        # do all those assignments for that
                                        sumDist_prv = sumDist_prv + nbrTable[2,nbrIndex_prv]
                                        # track the things we need to trim min/max
                                        if valueDiff_prv < min_Diff_prv:
                                            min_Diff_prv = valueDiff_prv
                                            minD_Dist_prv = nbrTable[2,nbrIndex_prv]
                                            minD_wpfv_prv = pfv_prv*weight_prv
                                            minD_weight_prv = weight_prv
                                        if valueDiff_prv > max_Diff_prv:
                                            max_Diff_prv = valueDiff_prv
                                            maxD_Dist_prv = nbrTable[2,nbrIndex_prv]
                                            maxD_wpfv_prv = pfv_prv*weight_prv
                                            maxD_weight_prv = weight_prv
                                        nfound_prv = nfound_prv + 1
                                        #nfound = nfound + 1
                                        if nfound_prv == _FILL_THRESHOLD:
                                            break
                                if nfound_prv == _FILL_THRESHOLD:
                                    break

                            if nfound_prv >= _FILL_MIN:
                                # we have found between min and threshold pairs
                                flag_prv = flag_prv | _SUCCESS_FLAG
                                if nfound_prv == _FILL_THRESHOLD:
                                    # we have found the upper limit of pairs (full fill)
                                    # (set both success flags)
                                    filledToThreshold += 1
                                    flag_prv = flag_prv | _SUCCESS_WAS_FULL_FLAG
                                else:
                                    filledBelowThreshold += 1
                                if _TRIM_MIN_MAX:
                                    ws_prv = ws_prv - (minD_wpfv_prv + maxD_wpfv_prv)
                                    sw_prv = sw_prv - (minD_weight_prv + maxD_weight_prv)
                                    sumDist_prv = sumDist_prv - (minD_Dist_prv + maxD_Dist_prv)
                                    nfound_prv = nfound_prv - 2
                                    #flag = 12 (old code recorded this per-pixel, don't see why)
                                wfv_prv = ws_prv / sw_prv
                                # record the result!!
                                outputData[z, y, x_prv] = wfv_prv
                                outputDists[z, y, x_prv] = sumDist_prv / nfound_prv
                            else:
                                flag_prv = flag_prv | _FAILURE_FLAG
                                # do not do anything to output: leave it at no-data
                                if nfound_prv > 0:
                                    insufficientPairsFound += 1
                                else:
                                    noPairsFound += 1
                            outputFlags [z, y, x_prv] = flag_prv
    resultDiagnostics = A1Diagnostics(
        TotalCellsSearched=totalCells,
        CellsWithGoodData=dataGood,
        OceanCells=oceanCells,
        NeverDataLocations=neverDataCells,
        GapCellsTotal=totalProcessedGapCells,
        GapCellsTooBig=gapsTooBig,
        PermanentGapCells=gapsAtUnfillableLocs,
        GapCellsFullyFilled=filledToThreshold,
        GapCellsPartFilled=filledBelowThreshold,
        FailedInsufficientPairs=insufficientPairsFound,
        FailedNoPairs=noPairsFound,
        TotalAlternateYearsScanned=scannedLevels,
        TotalNbrsChecked=scannedNeighbours,
        TotalNbrsUsed=usedNeighbours
    )
    dataStacks.DataArray3D = outputData
    # todo check this is valid
    dataStacks.DistanceTemplate3D = outputDists
    dataStacks.FlagsArray3D = outputFlags
    #return (outputData,outputDists,outputFlags,resultDiagnostics)
    return resultDiagnostics

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[::1] alternates_cy(int lim):
    ''' yields 1, -1, 2, -2, etc up to lim '''
    cdef int x
    x = 1
    cdef int[::1] arr = np.empty(shape=(lim*2),dtype=np.int32,order='c')
    for x in range(0,lim):
        arr[2*x] = x+1
        arr[2*x+1] = -(x+1)
        x += 1
    return arr


