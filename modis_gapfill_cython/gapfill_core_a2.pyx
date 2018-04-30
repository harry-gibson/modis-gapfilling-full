import numpy as np
cimport cython
from libc.math cimport sqrt
#cdef int _SEARCH_RADIUS = 10
from cython.parallel import prange, parallel

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef a2_core (
            dict DataImages,
            dict FlagValues,
            float _NDV,
            Py_ssize_t _A2_MAX_NBRS,
            char FillByRatios = 0,
            float RatioAbsZeroPoint = 0,
            float RatioLimit = 1
            ):
    '''
    Cython (C) implementation of the A2 gapfilling algorithm main "pass" code. This function should be called 8
    times with the data structured in such a way that the 8 different directional pass fills are generated, in
    order to generate overall A2 fill results. A separate (python) function is provided for this.

    The nature of the A2 algorithm is such that it cannot be easily parallelised - it modifies the input based
    on results from neighbouring cells, so the cells must be run in a deterministic order.

    Likewise A2 "drags" fill value ratios / calculations for an unlimited distance from nearest data pixels. This
    means that it cannot be run on separate tiles like A1 and must be run on global images.

    These two factors mean that unlike in the published paper, A2 is both slower and more demanding of memory than A1,
    whilst producing less good results.

    However it is still required if we need to be assured that _all_ gaps will be filled because A1 only works out to
    a specified radius (gap size). (In future we may try iterative / repeated running of A1 until all gaps are filled)

    This implementation is optimised as far as it can be in terms of Cython optimisations, except for the approach taken
    to the 8 different directional passes. To reduce memory use, the caller function (elsewhere) does not not make a
    C-ordered copy of the data but re-strides it. This greatly slows the A2 function (by a factor of around 6) and so
    on a machine with sufficent memory the data should be copied into the right order before passing to this function.
    '''
    # We assume that the arrays passed in are strided / transposed
    # such that iterating through in the standard c-normal order
    # will go through the source data in the correct order for this
    # directional pass. This means that some passes (the ones where the data
    # are in the native order) are several times quicker than the others
    cdef:
        # intermediate arrays
        double [:,::1] nbrTable
        int[:,::1] nbrIntCoords

        float[:,::1] diffImage_Local
        float [:, :] origDistImage_LocalCopy
        # in loop working vars
        Py_ssize_t y, x, yShape, xShape, nbrIndex, xNbr, yNbr
        double nbrDiffSum, nbrDiffCount, nbrDistSum, diffValThisPass
        double local_distance, nbrDist

        # metrics
        long long gotPixelVals = 0

        # unpack inputs to typed variables
        float[:,:] dataImage_Global_R = DataImages["Data"] # global in both senses
        unsigned char[:,:] flagsImage_Global_R = DataImages["Flags"]
        float[:,:] origDistImage_Global_R = DataImages["Distances"]
        float[:,:] meanImage_Global_R = DataImages["Means"]
        # these inputs get modified, i.e. they are "out" parameters in a proper language
        # It is done like this rather than having return values because they are actually going to be
        # strided views on arrays (to change iteration order)
        float[:,:] sumDistImage_Global_W = DataImages["SumDist"]
        float[:,:] outputImage_ThisPass_W = DataImages["Output"]
        char _FILL_FAILED_FLAG = FlagValues["FAILURE"]
        char _OCEAN_FLAG = FlagValues["OCEAN"]

        float _AbsZeroPoint = RatioAbsZeroPoint
        float _MaxAllowableRatio = RatioLimit
        float _MinAllowableRatio = 1.0 / _MaxAllowableRatio

    yShape = dataImage_Global_R.shape[0]
    xShape = dataImage_Global_R.shape[1]

    # it's actually the max neigbours value that defines how far out the search runs.
    # this just makes sure that the nbr table is generated far enough out.
    _SEARCH_RADIUS = <int>sqrt(_A2_MAX_NBRS / 3.14) + 10
    diam = _SEARCH_RADIUS * 2 + 1
    inds = np.indices([diam,diam]) - _SEARCH_RADIUS
    distTmp = np.sqrt((inds ** 2).sum(0))
    npTmpTable = ((inds.T).reshape(diam**2, 2))
    npTmpTable = np.append(npTmpTable, distTmp.ravel()[:,None],axis=1)

    # sort the table by distance then x then y (the arguments are last-sort-first)
    order = np.lexsort((npTmpTable[:,1],npTmpTable[:,0],npTmpTable[:,2]))
    npTmpTable = np.take(npTmpTable,order,axis=0)

    # the C-side result of the distance table calculations
    # transpose it to have three rows and many columns and take a C contiguous copy of this
    # so that access to individual nbr coord sets is optimised
    nbrTable = np.copy((npTmpTable[npTmpTable[:,2] <= _SEARCH_RADIUS]).T,order='c')
    # the distance table is stored with a float type but we need ints for indexing
    # based on its coords. we can cast at the time we get them out, but as this happens
    # in the innermost loop it is done trillions of times and so the time penalty of that
    # becomes signficant. So, store an int version of the coords array
    nbrIntCoords = np.asarray(nbrTable[0:2,:]).astype(np.int32)

    print ("Beginning pass of A2...")
    # populate the ratio or difference image,
    # this is just (data / mean) or (data - mean) (whether the data are original or from A1)
    diffImage_Local = np.empty_like(dataImage_Global_R)
    diffImage_Local[:] = _NDV
    with nogil, parallel(num_threads=20):
        #outerIdx = -1
        for y in prange(0, yShape):
            x = -1
            for x in range(0, xShape):
                if ((flagsImage_Global_R[y, x] & _OCEAN_FLAG) == _OCEAN_FLAG
                    # or meanImage_Global_R[y, x] == 0 # we need to be able to cope with mean = 0
                    or meanImage_Global_R[y, x] == _NDV
                    or dataImage_Global_R[y, x] == _NDV
                    ):
                    continue
                if FillByRatios == 0:
                    diffImage_Local[y, x] = (dataImage_Global_R[y, x] - meanImage_Global_R[y, x])
                else:
                    diffImage_Local[y, x] =((dataImage_Global_R[y, x] - _AbsZeroPoint)
                        / (meanImage_Global_R[y, x] - _AbsZeroPoint))
                    if diffImage_Local[y, x] > _MaxAllowableRatio:
                        diffImage_Local[y, x] = _MaxAllowableRatio
                    elif diffImage_Local[y, x] < _MinAllowableRatio:
                        diffImage_Local[y, x] = _MinAllowableRatio


    # the distances image gets modified within a pass but we don't want / need
    # to store it outside of the pass
    # However we are now creating the fresh copy in the caller function for clarity
    #origDistImage_LocalCopy = np.copy(origDistImage_Global_R)
    origDistImage_LocalCopy = origDistImage_Global_R

    with nogil:
        for y in range (0, yShape): # can't do parallel here,  boooo
            for x in range (0, xShape):
                #flag = flagsImage[y,x]
                if (flagsImage_Global_R[y, x] & _FILL_FAILED_FLAG) != _FILL_FAILED_FLAG:
                    # it's already good data, or filled with A1
                    continue
                if meanImage_Global_R[y,x] == _NDV:
                    #we can't fill, but, if we are here then the flag is already set
                    #to failure (from A1)... could optionally set a separate A2 failure flag (128)
                    continue
                # else...
                nbrDiffSum = 0
                nbrDiffCount = 0
                nbrDistSum = 0
                # summarise the values in the surrounding 8 pixels
                # we cannot precalculate this elementwise because the
                # offsets / ratios in diffImageThisPass are updated in the loop,
                # affecting later iterations, hence why we can't parallelise simply
                for nbrIndex in xrange(1, _A2_MAX_NBRS+1):
                    # +1 because the first row of nbr table is the current cell
                    # this was an attempt at loop reversal so we could always keep it c-contiguous -
                    # specify the order of the inner / outer loop as a parameter. Never got it working yet...
                    #if outerLoopShouldBe == 0:
                    xNbr = x + nbrIntCoords[0, nbrIndex]
                    yNbr = y + nbrIntCoords[1, nbrIndex]
                    if (xNbr >= 0 and xNbr < xShape and
                            yNbr >= 0 and yNbr < yShape and
                            diffImage_Local[yNbr, xNbr] != _NDV):
                        nbrDiffSum += diffImage_Local[yNbr, xNbr]
                        nbrDiffCount += 1
                        if origDistImage_LocalCopy[yNbr, xNbr] != _NDV:
                            nbrDist = origDistImage_LocalCopy[yNbr, xNbr]
                            # the distance of the filled cell will be the average of the distances
                            # of the neighbour cells used. Where "distance" of a neighbour cell
                            # will be the physical distance from the working cell, plus the
                            # distance already associated with the filled value in the nbr cell
                            # from the A1 algorithm, if applicable.
                            nbrDistSum += nbrDist
                        nbrDistSum += nbrTable[2, nbrIndex]

                # if any of the surrounding 8 pixels had a value then derive the cell
                # value from it
                if nbrDiffCount > 0 and (FillByRatios == 0 or nbrDiffSum > 0):
                    gotPixelVals += 1
                    # ratio / diff to use is the mean of the (up to) 8 values surrounding the cell
                    diffValThisPass = nbrDiffSum / nbrDiffCount
                    # fill in the same ratio image that we are checking in the neigbour search
                    # step such that the next cell along, in the current pass direction and if
                    # it's a gap, will have this calculated value as an available input.
                    # Hence values are "smeared" in direction of the scan.
                    # This means loop iterations are not independent (not embarrasingly
                    # parallel!) and thus this algorithm needs to be run on an entire
                    # global image.
                    diffImage_Local[y, x] = diffValThisPass
                    if FillByRatios == 0:
                        outputImage_ThisPass_W[y, x] = (diffValThisPass + meanImage_Global_R[y, x])
                    else:
                        outputImage_ThisPass_W[y, x] = (diffValThisPass *
                                                        (meanImage_Global_R[y, x] - _AbsZeroPoint)
                                                    + _AbsZeroPoint)
                else:
                    outputImage_ThisPass_W[y, x] = _NDV
                    # outside the function we can fill in flags image with _UNFILLABLE_AT_A2
                    # if all passes are nodata in finished image pass stack
    # nothing is returned, as the ratio and dist per-pass images are modified in-place
    #print "A2 filled {0!s} locations on this pass".format(gotPixelVals)
    return gotPixelVals
