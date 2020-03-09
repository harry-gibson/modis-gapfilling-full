import numpy as np
cimport cython
from libc.math cimport sqrt
#cdef int _SEARCH_RADIUS = 10
from cython.parallel import prange, parallel
from .gapfill_config_types import A2SearchConfig, DataLimitsConfig, FlagItems
from .gapfill_utils import  A2PassData

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef a2_core (
            dataStack: A2PassData,
            dataConfig: DataLimitsConfig,
            a2Config: A2SearchConfig,
            flagValues: FlagItems
            ):
    '''
    Cython (C) implementation of the A2 gapfilling algorithm main "pass" code. This function always iterates through 
    the data in C-normal order i.e. by column then row starting from top left. The function should be called 8
    times with the data transposed differently each time in such a way that the 8 different directional pass fills 
    are generated, in order to generate overall A2 fill results. A separate (python) function is provided for this.

    The nature of the A2 algorithm is such that it cannot be easily parallelised - it modifies the input based
    on results from neighbouring cells, so the cells must be run in a deterministic order.

    Likewise A2 "drags" fill value ratios / calculations for an unlimited distance from nearest data pixels, to ensure 
    that a fill value can always be generated (except on islands which have no data at all). 
    This means that it cannot be run on separate tiles like A1 and must be run on global images, hence the memory 
    requirements are very high. 

    These two factors mean that unlike in the published paper, A2 is both slower and more demanding of memory than A1,
    whilst producing less good results.

    However it is still required if we need to be assured that _all_ gaps will be filled because A1 only works out to
    a specified radius (gap size). (In future we may try iterative / repeated running of A1 until all gaps are filled)

    This implementation is optimised as far as it can be in terms of Cython optimisations, except for the approach taken
    to the 8 different directional passes. To reduce memory use, the caller function (elsewhere) does not not make a
    C-ordered copy of the data but re-strides it. This greatly slows the A2 function (by a factor of around 6) and so
    on a machine with sufficent memory the data should be copied into the right order before passing to this function.
    '''
    cdef:
        # intermediate arrays
        double [:,::1] nbrTable
        int[:,::1] nbrIntCoords

        float[:,::1] diffImage_Local
        float [:, :] origDistImage_LocalCopy
        # in loop working vars
        Py_ssize_t y, x, yShape, xShape, nbrIndex, xNbr, yNbr
        double nbrDiffSum, nbrDiffCount, nbrDistSum, diffValThisPass, distValThisPass
        double local_distance, nbrDist

        # metrics
        long long gotPixelVals = 0

        # unpack inputs to typed variables
        # nb we can use [:,::1] here (rather than [:,:] *iif* we have transformed the data with a np.copy rather than
        # just re-striding it in the caller function
        float[:,::1] dataImage = dataStack.DataArray2D
        unsigned char[:,::1] flagsImage = dataStack.FlagsArray2D
        float[:,::1] distImageIn = dataStack.DataArray2D
        float[:,::1] meanImage = dataStack.MeanArray2D
        # these inputs get modified, i.e. they are "out" parameters
        float[:,::1] outputImageDist = dataStack.SumOfPassDistancesArray2D
        float[:,::1] outputImageData = dataStack.DataArrayOutput2D

        char _FILL_FAILED_FLAG = flagValues.FAILURE
        char _OCEAN_FLAG = flagValues.OCEAN

        float _AbsZeroPoint = dataConfig.ABSOLUTE_ZERO_OFFSET
        float _MaxAllowableRatio = a2Config.MAX_ALLOWABLE_RATIO
        float _MinAllowableRatio = 1.0 / _MaxAllowableRatio
        unsigned char FillByRatios = a2Config.FILL_GENERATION_METHOD == "RATIO"
        float _NDV = dataConfig.NODATA_VALUE
        Py_ssize_t A2_MAX_NBRS = a2Config.SPIRAL_CONFIG.MAX_NBRS_TO_SEARCH

    yShape = dataImage.shape[0]
    xShape = dataImage.shape[1]

    # it's actually the max neigbours value that defines how far out the search runs.
    # this just makes sure that the nbr table is generated far enough out.
    _SEARCH_RADIUS = <int>sqrt(A2_MAX_NBRS / 3.14) + 10
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
    diffImage_Local = np.empty_like(dataImage)
    diffImage_Local[:] = _NDV
    with nogil, parallel(num_threads=20):
        #outerIdx = -1
        for y in prange(0, yShape):
            x = -1
            for x in range(0, xShape):
                if ((flagsImage[y, x] & _OCEAN_FLAG) == _OCEAN_FLAG
                    # or meanImage_Global_R[y, x] == 0 # we need to be able to cope with mean = 0
                    or meanImage[y, x] == _NDV
                    or dataImage[y, x] == _NDV
                    ):
                    continue
                if FillByRatios == 0:
                    diffImage_Local[y, x] = (dataImage[y, x] - meanImage[y, x])
                else:
                    diffImage_Local[y, x] =((dataImage[y, x] - _AbsZeroPoint)
                        / (meanImage[y, x] - _AbsZeroPoint))
                    if diffImage_Local[y, x] > _MaxAllowableRatio:
                        diffImage_Local[y, x] = _MaxAllowableRatio
                    elif diffImage_Local[y, x] < _MinAllowableRatio:
                        diffImage_Local[y, x] = _MinAllowableRatio


    # the distances image gets modified within a pass but we don't want / need
    # to store it outside of the pass
    # However we are now creating the fresh copy in the caller function for clarity
    #origDistImage_LocalCopy = np.copy(origDistImage_Global_R)
    with nogil:
        for y in range (0, yShape): # can't do parallel here,  boooo
            for x in range (0, xShape):
                #flag = flagsImage[y,x]
                if (flagsImage[y, x] & _FILL_FAILED_FLAG) != _FILL_FAILED_FLAG:
                    # it's already good data, or filled with A1
                    continue
                if meanImage[y, x] == _NDV:
                    #we can't fill, but, if we are here then the flag is already set
                    #to failure (from A1)... could optionally set a separate A2 failure flag (128)
                    continue
                # else...
                nbrDiffSum = 0
                nbrDiffCount = 0
                nbrDistSum = 0
                # summarise the values in the surrounding 8 pixels
                # we cannot precalculate this elementwise because the
                # offsets / ratios in diffImageLocal are updated in the loop,
                # affecting later iterations, hence why we can't parallelise simply
                for nbrIndex in range(1, A2_MAX_NBRS+1):
                    # +1 because the first row of nbr table is the current cell
                    xNbr = x + nbrIntCoords[0, nbrIndex]
                    yNbr = y + nbrIntCoords[1, nbrIndex]
                    if (0 <= xNbr < xShape and
                            0 <= yNbr < yShape and
                            diffImage_Local[yNbr, xNbr] != _NDV):
                        nbrDiffSum += diffImage_Local[yNbr, xNbr]
                        nbrDiffCount += 1
                        if distImageIn[yNbr, xNbr] != _NDV:
                            nbrDist = distImageIn[yNbr, xNbr]
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
                    distValThisPass = nbrDistSum / nbrDiffCount
                    # fill in the same ratio image that we are checking in the neighbour search
                    # step such that the next cell along, in the current pass direction and if
                    # it's a gap, will have this calculated value as an available input.
                    # Hence values are "smeared" in direction of the scan.
                    # This means loop iterations are not independent (not embarrasingly
                    # parallel!) and thus this algorithm needs to be run on an entire
                    # global image.
                    diffImage_Local[y, x] = diffValThisPass
                    outputImageDist[y, x] += distValThisPass
                    if FillByRatios == 0:
                        outputImageData[y, x] = (diffValThisPass + meanImage[y, x])
                    else:
                        outputImageData[y, x] = (diffValThisPass *
                                                        (meanImage[y, x] - _AbsZeroPoint)
                                                    + _AbsZeroPoint)
                else:
                    outputImageData[y, x] = _NDV
                    # outside the function we can fill in flags image with _UNFILLABLE_AT_A2
                    # if all passes are nodata in finished image pass stack
    # nothing is returned, as the ratio and dist per-pass images are modified in-place
    #print "A2 filled {0!s} locations on this pass".format(gotPixelVals)
    return gotPixelVals
