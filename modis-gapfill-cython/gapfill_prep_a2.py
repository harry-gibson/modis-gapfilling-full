import time
import numpy as np
import gc
from gapfill_core_clamp import MinMaxClip
from gapfill_core_a2 import a2_core

from MAP_Raster_Utilities.General_Raster_Funcs.TiffManagement import ReadAOI_PixelLims, SaveLZWTiff


def A2Caller(dataImage, flagsImage, distImage, meanImage):
    '''Runs the 8 directional passes of A2 for a single image and generates the output from the mean or median.

    The input parameters are modified in place and will end up as the filled arrays ready for writing'''

    if isinstance(meanImage, str):
        meanImageIn = ReadAOI_PixelLims(meanImage, None, None)
    else:
        assert isinstance(meanImage, np.ndarray)
        meanImageIn = meanImage

    if isinstance(dataImage, str):
        dataImageIn = ReadAOI_PixelLims(dataImage, None, None)
    else:
        assert isinstance(dataImage, np.ndarray)
        dataImageIn = dataImage

    if isinstance(distImage, str):
        distImageIn = ReadAOI_PixelLims(distImage, None, None)
    else:
        assert isinstance(distImage, np.ndarray)
        distImageIn = distImage

    if isinstance(flagsImage, str):
        flagsImageIn = ReadAOI_PixelLims(flagsImage, None, None)
    else:
        assert isinstance(flagsImage, np.ndarray)
        flagsImageIn= flagsImage

    global stdsPreloaded
    assert dataImageIn.shape == flagsImageIn.shape
    assert flagsImageIn.shape == distImageIn.shape

    assert distImageIn.shape == meanImageIn.shape

    start_time = time.time()

    #(dataImageIn + distImageIn + meanImageIn + dataPassStack*8 + sumDistImage + distImageLocal + ratioImageLocal)*4 + flagsImageIn
    # = 57 bytes per pixel required, 933.1M pixels, ~50Gb RAM required.
    # So A2 for 1km global Should _just_ be ok on a 64Gb machine! However in practice it gets memory errors because
    # by the time it runs the intermediate data from the last slice of A1 has not been fully flushed to disk by the OS.
    # If we have to, then we could reread data or dist images between passes, or discard the pass results and keep a single
    # running mean from the passes, if we are using a mean rather than median selector.

    dataPassStack = np.repeat(dataImageIn[np.newaxis, :, :], 8, axis=0) # gets modified in place

    sumDistImage = np.zeros_like(dataImageIn) # gets incremented in place
    global _NDV

    # we're not tinkering with this one, A2 will always run with 8 neighbours
    _A2_MAX_NBRS = 8

    # The A2 cython code just runs in the same order each time on what it gets (y outer, x inner, 0-n).
    # To get the 8 different directional passes we re-stride the inputs here, so that this causes them
    # to be passed over in each of the directions. This means that some passes are _much_ slower than
    # others as they aren't accessing the memory in the native contiguous order, but it saves us
    # making a physical copy each time
    for passNumber in xrange(0, 8):
        print "Running A2 pass "+str(passNumber)
        tDelta = time.time()
        # these will not be modified by A2, so they just get re-strided and passed in without copying
        dataGlobalImage = None
        flagsGlobalImage = None
        distGlobalImage = None
        meanGlobalImage = None
        distImageCopy = np.copy(distImageIn)
        # these ones get modified by A2, so they get restrided each time
        dataPassImage = None
        # this one gets modified by A2 and is a running total of all passes so does not get reset
        sumDistPassImage = sumDistImage
        outerLoop=0
        if passNumber == 0: # pass "a" takes ~50s for a global image
            dataGlobalImage = dataImageIn.T
            flagsGlobalImage = flagsImageIn.T
            distGlobalImage = distImageCopy.T
            meanGlobalImage = meanImageIn.T
            sumDistPassImage = sumDistImage.T
            dataPassImage = dataPassStack[passNumber].T
            outerLoop=1
        elif passNumber == 1: # pass "b" ~50s
            dataGlobalImage = dataImageIn.T[:,::-1]
            flagsGlobalImage = flagsImageIn.T[:,::-1]
            distGlobalImage = distImageCopy.T[:,::-1]
            meanGlobalImage = meanImageIn.T[:,::-1]
            sumDistPassImage = sumDistImage.T[:,::-1]
            dataPassImage = dataPassStack[passNumber].T[:,::-1]
            outerLoop=1
        elif passNumber == 2: # pass "c" ~50s
            dataGlobalImage = dataImageIn.T[::-1,:]
            flagsGlobalImage = flagsImageIn.T[::-1,:]
            distGlobalImage = distImageCopy.T[::-1,:]
            meanGlobalImage = meanImageIn.T[::-1,:]
            sumDistPassImage = sumDistImage.T[::-1,:]
            dataPassImage = dataPassStack[passNumber].T[::-1,:]
            outerLoop=1
        elif passNumber == 3: # pass "d" ~50s
            dataGlobalImage = dataImageIn.T[::-1,::-1]
            flagsGlobalImage = flagsImageIn.T[::-1,::-1]
            distGlobalImage = distImageCopy.T[::-1,::-1]
            meanGlobalImage = meanImageIn.T[::-1,::-1]
            sumDistPassImage = sumDistImage.T[::-1,::-1]
            dataPassImage = dataPassStack[passNumber].T[::-1,::-1]
            outerLoop=1
        elif passNumber == 4: # pass "e" - this one is the native c-ordering, ~7s!
            dataGlobalImage = dataImageIn
            flagsGlobalImage = flagsImageIn
            distGlobalImage = distImageCopy
            meanGlobalImage = meanImageIn
            sumDistPassImage = sumDistImage
            dataPassImage = dataPassStack[passNumber]

        elif passNumber == 5: # pass "f" ~9s
            dataGlobalImage = dataImageIn[:,::-1]
            flagsGlobalImage = flagsImageIn[:,::-1]
            distGlobalImage = distImageCopy[:,::-1]
            meanGlobalImage = meanImageIn[:,::-1]
            sumDistPassImage = sumDistImage[:,::-1]
            dataPassImage = dataPassStack[passNumber][:,::-1]

        elif passNumber == 6: # pass "g" ~8s
            dataGlobalImage = dataImageIn[::-1,:]
            flagsGlobalImage = flagsImageIn[::-1,:]
            distGlobalImage = distImageCopy[::-1,:]
            meanGlobalImage = meanImageIn[::-1,:]
            sumDistPassImage = sumDistImage[::-1,:]
            dataPassImage = dataPassStack[passNumber][::-1,:]

        else: #elif passNumber == 7: # pass "h" ~8s
            dataGlobalImage = dataImageIn[::-1,::-1]
            flagsGlobalImage = flagsImageIn[::-1,::-1]
            distGlobalImage = distImageCopy[::-1,::-1]
            meanGlobalImage = meanImageIn[::-1,::-1]
            sumDistPassImage = sumDistImage[::-1,::-1]
            dataPassImage = dataPassStack[passNumber][::-1,::-1]

        #print "a2 pass "+str(passNumber)+"...",

        a2_core (
             {
                "Data":      dataGlobalImage,
                "Flags":     flagsGlobalImage,
                "Distances": distGlobalImage,
                "Means":     meanGlobalImage,
                "SumDist":   sumDistPassImage,
                "Output":    dataPassImage
             },
             FlagValues,
             _NDV,
             _A2_MAX_NBRS,
             _FILL_BY_RATIO,
             _DATA_ZERO_OFFSET,
             _MAX_ALLOWABLE_RATIO
         )

        #print "done in "+str(time.time()-tDelta)+" seconds"

    tDelta = time.time()
    # process numpyisms on the A2 data stack tile-by-tile because it makes temp arrays,
    # and we probably can't afford that memory
    xCorners = np.linspace(0,dataImageIn.shape[1],6).astype(np.int32)
    yCorners = np.linspace(0,dataImageIn.shape[0],6).astype(np.int32)
    medianVals = np.empty_like(dataImageIn) # or mean, depending on choice
    countVals = np.empty(shape=dataImageIn.shape,dtype='byte')
    for x in xrange(len(xCorners)-1):
        for y in xrange(len(yCorners)-1):
            x0 = xCorners[x]
            x1 = xCorners[x+1]
            y0 = yCorners[y]
            y1 = yCorners[y+1]
            dataPassSlice = dataPassStack[:,y0:y1,x0:x1]
            countSlice = countVals[y0:y1,x0:x1]
            np.sum(dataPassSlice != _NDV, axis=0, out=countSlice)
            # exclude no-data from mean / median calculation by using nan-aware functions
            dataPassSlice[dataPassSlice == _NDV] = np.nan
            medianValsSlice = medianVals[y0:y1,x0:x1]
            if _A2_PASS_SELECTOR == "MEAN":
                medianValsSlice[:] = bn.nanmean(dataPassSlice, axis=0)
            else:
                medianValsSlice[:] = bn.nanmedian(dataPassSlice, axis=0)


    dataPassStack = None
    del dataPassStack # and breathe
    gc.collect()

    # A2 doesn't modify the flags image, so anywhere that the fill had failed is where it was
    # needed
    a2FillWasNeeded = np.bitwise_and(flagsImageIn,_FILL_FAILED_FLAG)==_FILL_FAILED_FLAG
    a2FillWasNeededCount = np.count_nonzero(a2FillWasNeeded)
    # we have a value for anywhere that the output isn't _NDV (now np.nan)
    a2FilledLocs = np.logical_and(a2FillWasNeeded,~np.isnan(medianVals))
    a2FilledCount = np.count_nonzero(a2FilledLocs)
    # so copy those a2 results into the "output" data image
    dataImageIn[a2FilledLocs] = medianVals[a2FilledLocs]
    flagsImageIn[a2FilledLocs] -= _FILL_FAILED_FLAG # should use bitwise ^ really
    flagsImageIn[a2FilledLocs] += _A2_SUCCESS_FLAG  # should use bitwise or
    # set distance metric to be the total A2 distance divided by the number of A2 passes it was from
    distImageIn[a2FilledLocs] = sumDistImage[a2FilledLocs] / countVals[a2FilledLocs]

    #print "Flattening done in "+str(time.time()-tDelta)+" seconds"
    print "A2 done, filled {0} locations of {1} required in {2} seconds using {3} of 8 passes".format(
                                                    a2FilledCount,a2FillWasNeededCount,
                                                    str(time.time()-start_time), _A2_PASS_SELECTOR)

