import time
import numpy as np
import gc
import bottleneck as bn
cimport cython

from gapfill_core_a2 import a2_core
from gapfill_utils import A2PassData
from gapfill_config import A2SearchConfig, FlagItems, DataCharacteristicsConfig, A2Diagnostics


def A2ImageCaller(dataImageIn, flagsImageIn, distImageIn, meanImageIn,
             A2SearchConfig a2Config, FlagItems flagValues, DataCharacteristicsConfig dataConfig):
    '''Runs the 8 directional passes of A2 for a single image and generates the output from the mean or median.

    The input parameters are modified in place and will end up as the filled arrays ready for writing'''

    assert isinstance(meanImageIn, np.ndarray)
    assert isinstance(dataImageIn, np.ndarray)
    assert isinstance(distImageIn, np.ndarray)
    assert isinstance(flagsImageIn, np.ndarray)

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

    # we need to duplicate the data image 8 times, once for each pass, because if we want to track the median of
    # the 8 passes then we can't do that with a running calculation
    dataPassStack = np.repeat(dataImageIn[np.newaxis, :, :], 8, axis=0) # gets modified in place

    sumDistImage = np.zeros_like(dataImageIn) # gets incremented in place
    global _NDV

    # we're not tinkering with this one, A2 will always run with 8 neighbours
    _A2_MAX_NBRS = a2Config.SPIRAL.MAX_NBRS_TO_SEARCH
    _FILL_FAILED_FLAG = flagValues.FAILURE
    _A2_SUCCESS_FLAG = flagValues.A2_FILLED

    # The A2 cython code just runs in the same order each time on what it gets (y outer, x inner, 0-n).
    # To get the 8 different directional passes we flip / mirror / transpose the inputs such that iterating in
    # this standard order is equivalent to going in a different direction over the original values.
    # In the notebook version of this code, we did this by simply re-striding the inputs here (basically creating
    # a "view" on the arrays that presents the data from the same memory location in a different order).
    # This means that some passes are _much_ slower than others as they aren't accessing the memory in the native
    # contiguous order, but it saves us making a physical copy each time
    # without copying, passes on a typical global image take approx:
    # 0: 50s, 1: 50s, 2: 50s, 3: 50s, 4: 7s, 5: 9s,  6: 8s, 7: 8s
    # In this version, the A2PassData object will actually create concrete copies of the data for iteration which
    # means that all passes should occur at the same speed, but more memory is needed. There are now enough years
    # in the MODIS record that the memory use of A1 has increased to the point that A2 is not really the dominant
    # memory hog any more.
    for passNumber in range(0, 8):
        print ("Running A2 pass "+str(passNumber))
        tDelta = time.time()
        dataPassImage = None
        a2Data = A2PassData(passNumber, dataImageIn, flagsImageIn, distImageIn, meanImageIn, sumDistImage)

        a2diag = a2_core(a2Data, dataConfig=dataConfig, a2Config=a2Config)
        dataPassStack[passNumber] = a2Data.getOutputData()
        sumDistImage += a2Data.getOutputDists()

        #print "done in "+str(time.time()-tDelta)+" seconds"

    tDelta = time.time()

    # process numpyisms on the A2 data stack tile-by-tile because it makes temp arrays,
    # and we probably can't afford that memory
    # create 6 * 6 tiled slices
    xCorners = np.linspace(0,dataImageIn.shape[1],6).astype(np.int32)
    yCorners = np.linspace(0,dataImageIn.shape[0],6).astype(np.int32)
    # create an output image for the average (mean or median) result overall
    passAverageImage = np.empty_like(dataImageIn) # or mean, depending on choice
    countVals = np.empty(shape=dataImageIn.shape,dtype='byte')
    for x in range(len(xCorners)-1):
        for y in range(len(yCorners)-1):
            x0 = xCorners[x]
            x1 = xCorners[x+1]
            y0 = yCorners[y]
            y1 = yCorners[y+1]
            # a view of the results stack
            dataPassSlice = dataPassStack[:,y0:y1,x0:x1]
            countSlice = countVals[y0:y1,x0:x1]
            np.sum(dataPassSlice != _NDV, axis=0, out=countSlice)
            # exclude no-data from mean / median calculation by using nan-aware functions
            dataPassSlice[dataPassSlice == _NDV] = np.nan
            passAverageSlice = passAverageImage[y0:y1,x0:x1]
            if a2Config.PASS_AVERAGE_TYPE == "MEAN":
                passAverageSlice[:] = bn.nanmean(dataPassSlice, axis=0)
            else:
                passAverageSlice[:] = bn.nanmedian(dataPassSlice, axis=0)

    dataPassStack = None
    del dataPassStack # and breathe
    gc.collect()

    # A2 doesn't modify the flags image, so anywhere that the fill had failed is where it was
    # needed

    a2FillWasNeeded = np.bitwise_and(flagsImageIn,_FILL_FAILED_FLAG)==_FILL_FAILED_FLAG
    a2FillWasNeededCount = np.count_nonzero(a2FillWasNeeded)
    # we have a value for anywhere that the output isn't _NDV (now np.nan)
    a2FilledLocs = np.logical_and(a2FillWasNeeded,~np.isnan(passAverageImage))
    a2FilledCount = np.count_nonzero(a2FilledLocs)

    # so copy those a2 results into the original data image (overwriting it)
    dataImageIn[a2FilledLocs] = passAverageImage[a2FilledLocs]
    flagsImageIn[a2FilledLocs] -= _FILL_FAILED_FLAG # should use bitwise ^ really
    flagsImageIn[a2FilledLocs] += _A2_SUCCESS_FLAG  # should use bitwise or
    # set distance metric to be the total A2 distance divided by the number of A2 passes it was from
    distImageIn[a2FilledLocs] = sumDistImage[a2FilledLocs] / countVals[a2FilledLocs]

    timeSeconds = time.time()-start_time
    # we just return diagnostics, the arrays have been modified in-place
    return A2Diagnostics(GapCellsTotal=a2FillWasNeededCount, GapCellsFilled=a2FilledCount, TimeSeconds=timeSeconds)
