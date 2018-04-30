
def class gapfiller:
    def __init__(self):


    def RunFill():


    if _CLIP_TO_LIMITS:
        # if we don't have memory to store global stds, we could load them again here
        # as they are only needed at this point so we don't need them in mem at same
        # time as the A2 stack
        # with rasterio.open(std_Synoptic_Path) as src:
        if stdsPreloaded is None:
            src = gdal.Open(std_Synoptic_Path)
            bnd = src.GetRasterBand(1)
            stds = bnd.ReadAsArray(xLims[0], yLims[0], totalFillWidth, totalFillHeight)
            bnd = None
            src = None
            # stds = src.read_band(1, window=((totalMarginT,totalMarginT+totalFillHeight),
            #                                        (totalMarginL,totalMarginL+totalFillWidth)),
            #                     masked=False)
        else:
            stds = stdsPreloaded
        MinMaxClip(dataImageIn, flagsImageIn, meanImageIn, stds,
                   _A2_SUCCESS_FLAG, _MINMAX_CLIP_FLAG, _FLOOR_CEILING_VALUE, _NDV,
                   _DATA_UPPER_LIMIT, _DATA_LOWER_LIMIT)
    # todo - clip to min/max?
