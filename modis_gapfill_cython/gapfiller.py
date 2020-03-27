# hack in the path to map raster utilities for I/O
import sys

# standard python libraries
import os
import math
import numpy as np
import itertools
import glob
from collections import defaultdict
# io management functions
import rasterio as rio
from rasterio.windows import Window
import rasterio.dtypes as rio_dt
import rasterio.transform
from rasterio.profiles import Profile as rio_prof
from osgeo import gdal
# from raster_utilities.tileProcessor import tileProcessor
from collections import namedtuple
# cython functions
from .gapfill_core_despeckle_and_flag import setSpeckleFlags
from .gapfill_core_a1 import a1_core
from .gapfill_core_clamp import MinMaxClip3D, MinMaxClip
from .gapfill_prep_a2 import A2ImageCaller

# configuration types (NamedTuples with yaml-parsing factory methods)
from .gapfill_config_types import GapfillFilePaths, GapfillJobConfig, DataLimitsConfig, \
    DespeckleConfig, A1SearchConfig, A2SearchConfig, FlagItems, RasterProps

from .gapfill_utils import PixelMargins, A1DataStack, A2PassData

class GapFiller:
    def __init__(self, gapfillFilePaths: GapfillFilePaths,
                 despeckleConfig: DespeckleConfig,
                 a1Config: A1SearchConfig,
                 a2Config: A2SearchConfig,
                 dataSpecificConfig: DataLimitsConfig,
                 flagValues: FlagItems,
                 jobDetails: GapfillJobConfig
                 ):

        assert isinstance(gapfillFilePaths, GapfillFilePaths)
        assert isinstance(despeckleConfig, DespeckleConfig)
        assert isinstance(a1Config, A1SearchConfig)
        assert isinstance(a2Config, A2SearchConfig)
        assert isinstance(dataSpecificConfig, DataLimitsConfig)
        assert isinstance(flagValues, FlagItems)
        assert isinstance(jobDetails, GapfillJobConfig)
        self._despeckleConfig = despeckleConfig
        self._a1Config = a1Config
        self._a2Config = a2Config
        self._dataSpecificConfig = dataSpecificConfig
        self._flagValues = flagValues
        self._jobDetails = jobDetails

        # initialise input files and index by calendar day and year
        self._inputFileDict = defaultdict(lambda: defaultdict(str))
        self.__InitialiseFiles(gapfillFilePaths.DATA_FILES_GLOB_PATTERN)
        self._filePaths = gapfillFilePaths
        self._outTemplate = self._jobDetails.OutputFilenameTag+"{}.{}.{}.TemporalSummary.Res.SpatialSummary.tif"
        self._intermediateFiles = defaultdict(dict)

        # initialise limits of fill area in pixel coords of the input files
        _lonLims = (jobDetails.XMin_Deg, jobDetails.XMax_Deg)
        _latLims = (jobDetails.YMax_Deg, jobDetails.YMin_Deg)

        self._setInputProperties()
        if None in _latLims and None in _lonLims:
            self.OutputProps = self.InputRasterProps
            self.xLims = (0, self.InputRasterProps.width)
            self.yLims = (0, self.InputRasterProps.height)

        else:
            self.xLims, self.yLims = self._getPixelLims(lonLims=_lonLims, latLims=_latLims)
            # TODO barring snapping just create this manually from the origin and existing resolution
            outGT = self._getClippedTransform(self.InputRasterProps.gt, xLims=self.xLims,
                                                 yLims=self.yLims)
            outW = self.xLims[1] - self.xLims[0]
            outH = self.yLims[1] - self.yLims[0]
            outProj = self.InputRasterProps.proj
            outNdv = self.InputRasterProps.ndv
            outDT = self.InputRasterProps.datatype
            self.OutputProps = RasterProps(gt=outGT, proj=outProj, ndv=outNdv, width=outW, height=outH,
                                           datatype=outDT)
        maxZSize = max([len(v) for k,v in self._inputFileDict.items()])
        # calculate slices to run a1
        self._slices = self.CalculateSliceEdges(maxZSize)

    def _getPixelLims(self, lonLims, latLims):
        # TODO consider replacing with rasterio.DatasetReader.index.
        # However this contains logic to deal with common issues in MAP around rounding of coordinates to
        # ensure properly aligned output, so may not do
        assert isinstance(lonLims, tuple) and len(lonLims) == 2
        assert isinstance(latLims, tuple) and len(latLims) == 2

        EastLimitOut = lonLims[1]
        WestLimitOut = lonLims[0]
        NorthLimitOut = latLims[0]
        SouthLimitOut = latLims[1]

        assert EastLimitOut > WestLimitOut
        assert NorthLimitOut > SouthLimitOut

        if None in latLims and None in lonLims:
            self.OutputProps = self.InputRasterProps
            self.xLims = (0, self.InputRasterProps.width)
            self.yLims = (0, self.InputRasterProps.height)
            return

        oneday = next(iter(self._inputFileDict.values()))
        onefile = next(iter(oneday.values()))
        # TODO fix incompatibility with the item positions of the Affine type
        with rio.open(onefile) as templateFile:
            inGT = templateFile.transform
            xRes = inGT[0]
            yRes = -(inGT[4])
            # check for the commonly-used imprecise resolutions and recalculate to ensure we don't end up missing a pixel
            inResXRnd = round(xRes, 8)
            inResYRnd = round(yRes, 8)
            if inResXRnd == 0.00833333:
                xRes = 1.0 / 120.0
            elif inResXRnd == 0.04166667:
                xRes = 1.0 / 24.0
            elif inResXRnd == 0.08333333:
                xRes = 1.0 / 12.0
            if inResYRnd == 0.00833333:
                yRes = 1.0 / 120.0
            elif inResYRnd == 0.04166667:
                yRes = 1.0 / 24.0
            elif inResYRnd == 0.08333333:
                yRes = 1.0 / 12.0

            OverallNorthLimit = inGT[5]
            OverallWestLimit = inGT[2]

            x0 = int((WestLimitOut - OverallWestLimit) / xRes)
            x1 = int(((EastLimitOut - OverallWestLimit) / xRes) + 0.5)
            y0 = int((OverallNorthLimit - NorthLimitOut) / yRes)
            y1 = int(((OverallNorthLimit - SouthLimitOut) / yRes) + 0.5)
            return((x0, x1),(y0, y1))

    def _setInputProperties(self):
        oneday = next(iter(self._inputFileDict.values()))
        onefile = next(iter(oneday.values()))
        with rio.open(onefile) as firstFile:
            initialProps = RasterProps(gt=firstFile.transform,  # an Affine type, not just a tuple
                    proj=firstFile.crs,  # a CRS type, not just a string
                    ndv=firstFile.nodata,  # not sure what this does if file is multiband,
                    # maybe should do get_nodata()[0]
                    width=firstFile.width, height=firstFile.height,
                    datatype=firstFile.dtypes[0]
                    )
        self.InputRasterProps = initialProps

    def _getClippedTransform(self, inGT, xLims, yLims):
        '''Returns the GDAL geotransform of a clipped subset of an existing geotransform

        Where clipping coordinates are specified as pixel limits'''
        topLeftLongIn = inGT[2]
        topLeftLatIn = inGT[5]
        resX = inGT[0]
        resY = inGT[4]

        # the origin coord of the output is simply a whole number of pixels from the input origin
        topLeftLongOut = topLeftLongIn + xLims[0] * resX
        topLeftLatOut = topLeftLatIn + yLims[0] * resY

        clippedGT = rasterio.transform.from_origin(topLeftLongOut, topLeftLatOut, resX, -resY)
        return clippedGT
    def RunFill(self):
        allDays = self._inputFileDict.keys()
        startYear = self._jobDetails.StartYear
        if startYear is None:
            startYear = 0  # just some value smaller than the first year
        onlyForDays = self._jobDetails.CalendarDaysToFill
        for calendarDay in allDays:
            if onlyForDays is not None and int(calendarDay) not in onlyForDays:
                # if we're not filling for this calendar day, we don't need to do anything
                continue
            for slice in self._slices:
                # years are slightly more complex as we always need to read the day's data for all years
                # this populates the self._intermediateFiles dictionary
                self.DespeckleAndA1SliceRunner(slice, calendarDay, startYear)
            self.A2BatchRunner()

    def __InitialiseFiles(self, globPattern):
        dayMonths = {1: 1, 9: 1, 17: 1, 25: 1, 33: 2, 41: 2, 49: 2, 57: 2, 65: 3, 73: 3, 81: 3, 89: 3, 97: 4, 105: 4,
                     113: 4, 121: 4,
                     129: 5, 137: 5, 145: 5, 153: 6, 161: 6, 169: 6, 177: 6, 185: 7, 193: 7, 201: 7, 209: 7, 217: 8,
                     225: 8, 233: 8,
                     241: 8, 249: 9, 257: 9, 265: 9, 273: 9, 281: 10, 289: 10, 297: 10, 305: 10, 313: 11, 321: 11,
                     329: 11, 337: 12,
                     345: 12, 353: 12, 361: 12}
        calendarDays = [str(i).zfill(3) for i in np.arange(1, 365, 8)]

        allFiles = glob.glob(globPattern)

        for f in allFiles:
            calendarday, year = self.__parseDailyFilename(f)
            self._inputFileDict[calendarday][year] = f

    def __parseDailyFilename(self, f):
        """Get julian day and year from a filename in old or new MAP MODIS filename formats,
        and set instance variable _outTemplate to match the non-temporally-varying parts of this
        input file (or, on subsequent calls, check these are the same as last time).
        _outTemplate will be set to e.g.
            LST_Day{}.{}.{}.Data.1km.Data"""
        base = os.path.basename(f)

        tokens = base.split('.')
        if len(tokens) < 6:
            # assume it's an old file in the format A2000089etcetc.tif i.e. ?YYYYDDD*
            yr = base[1:5]
            day = base[5:8]
        else:
            # assume it's a file in the newer format ?*.YYYY.DDD.etc format
            varname, yr, day, temporalSummary, res, spatialSummary = tokens[0:6]
            outTemplate = varname + "{}.{}.{}." + "{}.{}.{}.tif".format(temporalSummary, res, spatialSummary)
            if self._outTemplate == "FILLED-OUTPUT{}.{}.{}.TemporalSummary.Res.SpatialSummary.tif":
                self._outTemplate = outTemplate
            else:
                assert self._outTemplate == outTemplate
        return day, yr

    def __CalculateSliceSize(self, readShapeZYX):
        """Calculate the X-size of the slices we can run for the given Y, Z dimensions and mem limit"""
        # readShapeZYX is the dimension of the data we must READ to fill the required output area;
        #  i.e .the fill area plus margins. If we're filling globally it's the same thing.
        dataBPP = 4
        memLimit = self._jobDetails.MemTargetBytes
        outputsBPP = dataBPP * 2 + 1  # the output data, distances, and flags
        # approximate total number of pixels we can read for each file
        sliceSqrd = memLimit / (readShapeZYX[0] * (dataBPP + outputsBPP))
        # not implementing slicing in y dimension so xSize is total pixels / total height
        sliceXSize = sliceSqrd / readShapeZYX[1]
        return sliceXSize

    def CalculateSliceEdges(self, maxZSize):
        """
        Generate the slice boundaries, slicing only in the x dimension
        (we are using full data in y dimension for now, though this doesn't have to be so)
        The slice boundaries overlap by 2* _SEARCH_RADIUS pixels (_SEARCH_RADIUS on each slice)
        This allows the gap-filling code to run on the non-overlapping section of the slice
        while having all the data necessary to fill a gap up to the edge of that non-overlapping

        :return a list of 2-tuples where each part is a 6-tuple, the first representing the x-edges of this job's
        searches and the second the y-edges. Within each 6-tuple the items are the pixel coordinates for this dimension,
        relative to the input data files, of:
        (W/N edge of despeckle data, W/N edge of A1 data, W/N edge of A1 Output,
        E/S edge of A1 Output, E/S edge of A1 data, E/S edge of Despeckle data)
        """
        # always round up as there's no harm done in passing slightly too much data
        _A1_SEARCH_RADIUS = int(math.sqrt(self._a1Config.SPIRAL_CONFIG.MAX_NBRS_TO_SEARCH / 3.14) + 1)
        _DESPECKLE_SEARCH_RADIUS = int(math.sqrt(self._despeckleConfig.SPIRAL_CONFIG.MAX_NBRS_TO_SEARCH / 3.14) + 1)
        _TOTAL_DESIRED_MARGIN = _A1_SEARCH_RADIUS + _DESPECKLE_SEARCH_RADIUS

        reqDataL = max(0, self.xLims[0] - _TOTAL_DESIRED_MARGIN)
        reqDataR = min(self.xLims[1] + _TOTAL_DESIRED_MARGIN, self.InputRasterProps.width)
        reqDataT = max(0, self.yLims[0] - _TOTAL_DESIRED_MARGIN)
        reqDataB = min(self.yLims[1] + _TOTAL_DESIRED_MARGIN, self.InputRasterProps.height)
        reqDataHeight = reqDataB - reqDataT
        reqDataWidth = reqDataR - reqDataL

        # get an estimated maximum width for the slices based on the available memory and the number of years (vertical
        # size of the stack) and height (y dimension) of the data (not slicing in Y dimension)
        readShapeZYX = (maxZSize, reqDataHeight, reqDataWidth)
        sliceXSize = self.__CalculateSliceSize(readShapeZYX=readShapeZYX)

        # how many slices are we going to need
        totalFillWidth = self.xLims[1] - self.xLims[0]
        nChunks = int((totalFillWidth // sliceXSize) + 1)

        # generate the "chunk" boundaries that will represent the data processed in one thread's job, in terms of the
        # pixel coordinates of the source data files.
        # the actual slice boundaries for filling and writing out
        chunkEdges = np.linspace(self.xLims[0], self.xLims[1], nChunks + 1).astype(np.int32)
        leftRealEdges = chunkEdges[:-1]
        rightRealEdges = chunkEdges[1:]
        # the boundaries of these slices plus the margin of extra data for A1, but respecting the fact we can't
        # go beyond the edge of the source data
        left_A1_edges = np.clip((chunkEdges - _A1_SEARCH_RADIUS)[:-1],
                                0, np.inf).astype(np.int32)
        right_A1_edges = np.clip((chunkEdges + _A1_SEARCH_RADIUS)[1:],
                                 -np.inf, self.InputRasterProps.width).astype(np.int32)
        # the boundary of these slices plus the margin of extra data for A1 and despeckle, but respecting the fact
        # we can't go beyond the edge of the source data
        left_Despeckle_edges = np.clip((chunkEdges - _TOTAL_DESIRED_MARGIN)[:-1],
                                       0, np.inf).astype(np.int32)
        right_Despeckle_edges = np.clip((chunkEdges + _TOTAL_DESIRED_MARGIN)[1:],
                                        -np.inf, self.InputRasterProps.width).astype(np.int32)

        # the left and right task boundaries are _SEARCH_RADIUS bigger than the data that will be searched
        # within them so that all pixels can have neighbours, if possible (not at global edge)
        x_offsets_overlapping = zip(left_Despeckle_edges, left_A1_edges, leftRealEdges,
                                    rightRealEdges, right_A1_edges, right_Despeckle_edges)

        # can we pad the top and bottom? (not if we are doing a global run as y-slicing isn't implemented)
        topA1Edge = np.clip(self.yLims[0] - _A1_SEARCH_RADIUS, 0, np.inf).astype(np.int32)
        bottomA1Edge = np.clip(self.yLims[1] + _A1_SEARCH_RADIUS,
                               -np.inf, self.InputRasterProps.height).astype(np.int32)
        topDespeckleEdge = np.clip(self.yLims[0] - _TOTAL_DESIRED_MARGIN, 0, np.inf).astype(np.int32)
        bottomDespeckleEdge = np.clip(self.yLims[1] + _TOTAL_DESIRED_MARGIN,
                                      -np.inf, self.InputRasterProps.height).astype(np.int32)

        # create a list of the slices needed for each calendar day, i.e. all the x slices and for each of these,
        # all the y slices. Each dimension and slice coords are recorded as a 6-tuple which is a member of a list.
        # As there is only 1 y slice for now, these are a single-item list.
        # NB in the notebook version we also added the days themselves to the cartesian product for an overall
        # tasklist, here we are not bothering with that part
        taskListForEachDay = list(itertools.product(
            x_offsets_overlapping,
            [(topDespeckleEdge, topA1Edge, self.yLims[0], self.yLims[1], bottomA1Edge, bottomDespeckleEdge)]))
        return taskListForEachDay

    # = "despeckleAndA1Caller" - initial part
    def CalculateAvailableSliceMargins(self, sliceInfo):
        """Calculates the windows to read from the input, and to write to the output, for a given slice
        These are not the same thing, because to fill a given output area needs a double (concentric) margin around it
        from which to draw despeckle and then fill neighbour values.
        But these margins can't be achieved at the edge of the images. So we need to read what we can, then tell
        the despeckle and fill algorithms what margins to actually use on each side of the array, the final output
        after despeckle and fill will be without margins and will be the piece we actually write to the output."""

        # Warning headfuck ensues with regard to which data has which margins!...:
        #
        # A1 needs a margin outside the area it is filling, if possible
        # (i.e. if it's not a global image, which it isn't because we wouldn't have enough mem for that, and
        # where this slice isn't at the true edge i.e. 180deg E/W or 90deg N/S)
        #
        # Despeckle needs a further margin outside of that, if possible, so that A1 gets margin data that has
        # itself already been despeckled.
        #
        # So we want to read a section that covers the despeckle plus A1 margins, but then only fill the non-margin
        # parts. However at the edges of the input data we can't supply a margin.
        #
        # A sliceInfo is structured like this:
        # ((21720, 21740, 21840, 23040, 23140, 23160),
        # (9540, 9560, 9660, 10560, 10660, 10680))
        # The first 6-tuple is the x-margins of the slice in pixel coords relative to the overall image; the second
        # gives the y-margins.
        # Items 0 and 5 are the pixel coords of the despeckle margin of the slice (they define the data that should be
        # passed to despeckle). Items 1 and 4 define the data that should be passed to A1. Items 2 and 3 define the data
        # that will be filled by A1.
        #

        xInfo = sliceInfo[0]
        yInfo = sliceInfo[1]
        g_sliceDespeckleL, g_sliceA1L, g_sliceFillL, g_sliceFillR, g_sliceA1R, g_sliceDespeckleR = xInfo
        g_sliceDespeckleT, g_sliceA1T, g_sliceFillT, g_sliceFillB, g_sliceA1B, g_sliceDespeckleB = yInfo

        # the total width that will be gapfilled and written out
        sliceFillWidth = g_sliceFillR - g_sliceFillL
        # the total width of data we need to read, equal to reqDataWidth if only 1 slice
        sliceDespeckleWidth = g_sliceDespeckleR - g_sliceDespeckleL
        # the total width that will be passed to a1, not actually needed
        sliceA1Width = g_sliceA1R - g_sliceA1L

        # as above but for height, equal to reqDataHeight if only 1 slice
        sliceFillHeight = g_sliceFillB - g_sliceFillT
        sliceDespeckleHeight = g_sliceDespeckleB - g_sliceDespeckleT
        sliceA1Height = g_sliceA1B - g_sliceA1T
        assert sliceFillHeight == self.yLims[1] - self.yLims[0]  # not actually supporting y slices for now...

        # at the edges of the input, a margin is not possible... so how big a margin does despeckle
        # actually have to work with?
        sliceDespeckleMarginL = int(g_sliceA1L - g_sliceDespeckleL)
        sliceDespeckleMarginR = int(g_sliceDespeckleR - g_sliceA1R)
        sliceDespeckleMarginT = int(g_sliceA1T - g_sliceDespeckleT)
        sliceDespeckleMarginB = int(g_sliceDespeckleB - g_sliceA1B)
        assert ((sliceDespeckleMarginL >= 0) and (sliceDespeckleMarginR >= 0)
                and (sliceDespeckleMarginR >= 0) and (sliceDespeckleMarginB >= 0))

        # likewise how big a margin does A1 actually have to work with?
        sliceA1MarginL = g_sliceFillL - g_sliceA1L
        sliceA1MarginR = g_sliceA1R - g_sliceFillR
        sliceA1MarginT = g_sliceFillT - g_sliceA1T
        sliceA1MarginB = g_sliceA1B - g_sliceFillB
        assert ((sliceA1MarginL >= 0) and (sliceA1MarginR >= 0)
                and (sliceA1MarginR >= 0) and (sliceA1MarginB >= 0))

        # the coords within the input data that can be used, where the output (filled)
        # data will sit
        sliceTotalMarginT = int(sliceA1MarginT + sliceDespeckleMarginT)
        sliceTotalMarginB = int(sliceA1MarginB + sliceDespeckleMarginB)
        sliceTotalMarginL = int(sliceA1MarginL + sliceDespeckleMarginL)
        sliceTotalMarginR = int(sliceA1MarginR + sliceDespeckleMarginR)

        # if we are not running the whole globe, our output files do not have
        # global pixel coords but local ones relative to the total processing size
        # E.g. "Africa" starts at global X of 18467 but this translates to X of zero in the output image
        # Hence get the coords of this slice within the output space coordiantes, prefixed with out_
        # - i.e. where the output slice fits into the overall output image
        out_SliceFillL = g_sliceFillL - self.xLims[0]
        out_SliceFillR = out_SliceFillL + sliceFillWidth
        out_SliceFillT = g_sliceFillT - self.yLims[0]
        out_SliceFillB = out_SliceFillT + sliceFillHeight
        # and as we aren't yet supporting slices in the vertical dimension:
        assert out_SliceFillT == 0

        despeckleMargins = PixelMargins(
            top=sliceDespeckleMarginT, bottom=sliceDespeckleMarginB,
            left=sliceDespeckleMarginL, right=sliceDespeckleMarginR)

        a1Margins = PixelMargins(
            top=sliceTotalMarginT,
            bottom=sliceTotalMarginB,
            left=sliceTotalMarginL,
            right=sliceTotalMarginR
        )

        dataReadWindow = PixelMargins(
            top=g_sliceDespeckleT, bottom=g_sliceDespeckleB,
            left=g_sliceDespeckleL, right=g_sliceDespeckleR
        )

        dataWriteWindow = PixelMargins(
            top=out_SliceFillT, bottom=out_SliceFillB,
            left=out_SliceFillL, right=out_SliceFillR
        )

        return {
            "despeckleMargins": despeckleMargins,
            "a1Margins": a1Margins,
            "dataReadWindow": dataReadWindow,
            "dataWriteWindow": dataWriteWindow
        }

    # ="despeckleAndA1Caller" except for initial calcs
    def DespeckleAndA1SliceRunner(self, sliceInfo, calendarDay: int, startYear: set):
        """Runs despeckle and A1 gapfill across all required years for a given calendar day,
        for the slice specified.
        Specify a value for firstYearToFill to only fill data from that and subsequent years."""

        margins = self.CalculateAvailableSliceMargins(sliceInfo)
        despeckleMargins = margins["despeckleMargins"]
        a1Margins = margins["a1Margins"]
        sliceFillWidth = sliceInfo[0][3] - sliceInfo[0][2]
        sliceFillHeight = sliceInfo[1][3] - sliceInfo[1][2]
        dataReadWindow = margins["dataReadWindow"]
        dataWriteWindow = margins["dataWriteWindow"]
        if not calendarDay in self._inputFileDict:
            raise FileNotFoundError("No data files were identified for the requested calendar day ({})"
                                    .format(calendarDay))
        sliceGT = self._getClippedTransform(self.OutputProps.gt,
                                               xLims=(dataWriteWindow.Left, dataWriteWindow.Right),
                                               yLims=(dataWriteWindow.Top, dataWriteWindow.Bottom))
        sliceDespeckleHeight = dataReadWindow.Bottom - dataReadWindow.Top
        sliceDespeckleWidth = dataReadWindow.Right - dataReadWindow.Left
        # replaced reading with rasterio
        dataReadWindow_rio = Window.from_slices(rows=(dataReadWindow.Top, dataReadWindow.Bottom),
                                                cols=(dataReadWindow.Left, dataReadWindow.Right))
        with rio.open(self._filePaths.SYNOPTIC_MEAN_FILE) as meanReader:
            sliceMeanArr = meanReader.read(1, window=dataReadWindow_rio)
        with rio.open(self._filePaths.SYNOPTIC_SD_FILE) as sdReader:
            sliceSDArr = sdReader.read(1, window=dataReadWindow_rio)
        with rio.open(self._filePaths.COASTLINE_FILE) as coastReader:
            sliceCoastArr = coastReader.read(1, window=dataReadWindow_rio)
        yearFileDictThisDay = self._inputFileDict[calendarDay]
        sortedFiles = [f for (y, f) in sorted(yearFileDictThisDay.items())]
        if startYear != 0:
            try:
                startFillFromPos = sorted(yearFileDictThisDay.keys()).index(str(int(startYear)))
            except ValueError:
                raise FileNotFoundError("No data files were identified for the requested calendar day and year({}/{})"
                                        .format(startYear, calendarDay))
        else:
            startFillFromPos = 0

        # inputFileDict is a nested dict of dicts, maps {day : {year : file}}
        # but in this version of the code we are only called for one day at once
        # fileList = [file for yr, file in sorted(yearFileDict.items())]
        # we need to read all data regardless, as the algorithm can still draw on earlier years
        # to derive fill values
        dataStack = A1DataStack()
        dataStack.FillFromZPosition = startFillFromPos
        dataStack.DataArray3D = np.empty(shape=(len(sortedFiles), sliceDespeckleHeight, sliceDespeckleWidth),
                                         dtype='Float32')
        dataStack.FlagsArray3D = np.zeros(shape=(len(sortedFiles), sliceDespeckleHeight, sliceDespeckleWidth),
                                          dtype=np.uint8)
        dataStack.MeansArray2D = sliceMeanArr
        dataStack.SDArray2D = sliceSDArr
        dataStack.Coastline2D = sliceCoastArr
        dataStack.DistanceTemplate3D = None
        dataStack.KnownUnfillable2D = None
        for zPos in range(len(sortedFiles)):
            with rio.open(sortedFiles[zPos]) as f:
                dataStack.DataArray3D[zPos] = f.read(1, window=dataReadWindow_rio)
                assert f.nodata == self._dataSpecificConfig.NODATA_VALUE
        #assert dataStack.DataArray3D.flags.c_contiguous

        # Run the despeckle
        # dataStack.DataArray3D and dataStack.FlagsArray3D are modified in place. They will have a different
        # (smaller) shape as the speckle search margins will be removed.
        # the returned object is an instance of DespeckleDiagnostics
        # Note that despeckle runs on the whole stack, not taking account of firstYearToFill, this is
        # because we don't want extreme values to remain to be usable in deriving fill values
        print(self._despeckleConfig.GetSummaryMessage())
        despeckleResult = setSpeckleFlags(dataStacks=dataStack, margins=despeckleMargins,
                                          flagValues=self._flagValues, dataConfig=self._dataSpecificConfig,
                                          speckleConfig=self._despeckleConfig, nCores=self._jobDetails.NCores)

        print(despeckleResult.GetSummaryMessage())
        # Run A1
        # dataStack.DataArray3D, dataStack.FlagsArray3D, and dataStack.DistanceTemplate3D are modified in place.
        # They will have a different (smaller) shape as the A1 search margins will be removed.
        # The returned object is an instance of A1Diagnostics
        print(self._a1Config.GetSummaryMessage())
        a1Result = a1_core(dataStacks=dataStack, margins=a1Margins,
                           flagValues=self._flagValues, dataConfig=self._dataSpecificConfig,
                           a1Config=self._a1Config,
                           nCores=self._jobDetails.NCores)
        print(a1Result.GetSummaryMessage())
        # dataStack has been modified in place: DataArray32, FlagsArray3D and DistanceTemplate3D members
        # now contain the results

        if not ((self._dataSpecificConfig.CEILING_VALUE == self._dataSpecificConfig.NODATA_VALUE) and
                self._dataSpecificConfig.FLOOR_VALUE == self._dataSpecificConfig.NODATA_VALUE):
            # apply clip if either the floor / ceiling values are set
            # TODO add the heights
            sliceMeanArr = np.copy(sliceMeanArr[a1Margins.Top:a1Margins.Top + sliceFillHeight,
                                   a1Margins.Left:a1Margins.Left + sliceFillWidth])
            sliceSDArr = np.copy(sliceSDArr[a1Margins.Top:a1Margins.Top + sliceFillHeight,
                                 a1Margins.Left:a1Margins.Left + sliceFillWidth])
            dataStack.MeansArray2D = sliceMeanArr
            dataStack.SDArray2D = sliceSDArr
            # TODO set this to respect firstYearToFill
            MinMaxClip3D(dataStacks=dataStack,
                         flagToCheck=self._flagValues.A1_FILLED,
                         flagToSet=self._flagValues.CLIPPED,
                         dataConfig=self._dataSpecificConfig,
                         nCores=self._jobDetails.NCores)

        # TODO fix this
        if True:  # len(self._slices) > 1:
            # use intermediate files.
            # in this code we're not saving separate tile files, we're saving global tiffs and just writing
            # part of them at once
            dataWriteWindow_rio = Window.from_slices(rows=(dataWriteWindow.Top, dataWriteWindow.Bottom),
                                                     cols=(dataWriteWindow.Left,dataWriteWindow.Right))
            dataWriteProfile = self.OutputProps.GetRasterioProfile()
            for y in range(len(sortedFiles)):
                if y < startFillFromPos:
                    continue
                inputFilename = os.path.basename(sortedFiles[y])

                day, yr = self.__parseDailyFilename(inputFilename)
                outDataFile = self._outTemplate.format("-A1IntermediateData", yr, day)
                outDistFile = self._outTemplate.format("-A1IntermediateDists", yr, day)
                outFlagFile = self._outTemplate.format("-A1IntermediateFlags", yr, day)

                outNameData = os.path.join(self._filePaths.TEMP_FOLDER, outDataFile)
                outNameDist = os.path.join(self._filePaths.TEMP_FOLDER, outDistFile)
                outNameFlag = os.path.join(self._filePaths.TEMP_FOLDER, outFlagFile)
                if not os.path.exists(self._filePaths.TEMP_FOLDER):
                    os.makedirs(self._filePaths.TEMP_FOLDER)
                self._intermediateFiles[inputFilename]["Data"] = outNameData
                self._intermediateFiles[inputFilename]["Distances"] = outNameDist
                self._intermediateFiles[inputFilename]["Flags"] = outNameFlag
                with rio.open(outNameData, 'w',
                              **dataWriteProfile) as dst:
                    dst.write(np.asarray(dataStack.DataArray3D[y]), 1,
                              window=dataWriteWindow_rio)
                with rio.open(outNameDist, 'w',
                              **dataWriteProfile) as dst:
                    dst.write(np.asarray(dataStack.DistanceTemplate3D[y]), 1,
                              window=dataWriteWindow_rio)
                dataWriteProfile.update(dtype=rio_dt.uint8)
                dataWriteProfile.update(nodata=None)
                with rio.open(outNameFlag, 'w',
                              **dataWriteProfile) as dst:
                    dst.write(np.asarray(dataStack.FlagsArray3D[y]), 1,
                              window=dataWriteWindow_rio)


    def A2BatchRunner(self):
        # cf. a2Caller, but run for a provided filename
        dataReadWindow_rio = Window.from_slices(rows=self.yLims, cols=self.xLims)
        # the MAP writer code checks the array is of an equal or "higher" datatype based on their gdal numeric
        # enum values (which are in a suitable order for this to be useful). We can recreate those from rio using the
        # order of this list:
        ordered_rio_dtypes = list(rio_dt.dtype_fwd.values())
        with rio.open(self._filePaths.SYNOPTIC_MEAN_FILE) as f:
            arr_Mean = f.read(1, window=dataReadWindow_rio)
        if False: # self._CACHE_SD:
            with rio.open(self._filePaths.SYNOPTIC_SD_FILE) as f:
                arr_SD = f.read(1, window=dataReadWindow_rio)
        for inputName, intermediateDataForDate in self._intermediateFiles.items():
            dataFile = intermediateDataForDate["Data"]
            with rio.open(dataFile) as f:
                # the entirety of this intermediate file corresponds to the window read from the
                # global synoptic files above
                arr_Data = f.read(1)
                # reusing the RasterProps namedtuple type from MAP-raster-utilities just for convenience
                props_Data = RasterProps(gt=f.transform, # an Affine type, not just a tuple
                                         proj=f.crs, # a CRS type, not just a string
                                         ndv=f.nodata,  # not sure what this does if file is multiband,
                                                        # maybe should do get_nodata()[0]
                                         width=f.width, height=f.height,
                                         datatype=f.dtypes[0]
                                         )
            distsFile = intermediateDataForDate["Distances"]
            with rio.open(distsFile) as f:
                arr_Dists = f.read(1)
                # reusing the RasterProps namedtuple type from MAP-raster-utilities just for convenience
                props_Dists = RasterProps(gt=f.transform,  # an Affine type, not just a tuple
                                         proj=f.crs,  # a CRS type, not just a string
                                         ndv=f.nodata,  # not sure what this does if file is multiband,
                                         # maybe should do get_nodata()[0]
                                         width=f.width, height=f.height,
                                         datatype=f.dtypes[0]  # a string, not a numeric gdal type code
                                         )

            flagsFile = intermediateDataForDate["Flags"]
            with rio.open(flagsFile) as f:
                arr_Flags = f.read(1)
                # reusing the RasterProps namedtuple type from MAP-raster-utilities just for convenience
                props_Flags = RasterProps(gt=f.transform,  # an Affine type, not just a tuple
                                         proj=f.crs,  # a CRS type, not just a string
                                         ndv=f.nodata,  # not sure what this does if file is multiband,
                                         # maybe should do get_nodata()[0]
                                         width=f.width, height=f.height,
                                         datatype=f.dtypes[0]  # a string, not a numeric gdal type code
                                         )
            if self._jobDetails.RunA2:
                print(f"Running A2 for image {dataFile}")
                print(self._a2Config.GetSummaryMessage())
                A2Diagnostics = A2ImageCaller(dataImageIn=arr_Data, flagsImageIn=arr_Flags, distImageIn=arr_Dists,
                                              meanImageIn=arr_Mean,
                                              a2Config=self._a2Config, flagValues=self._flagValues,
                                              dataConfig=self._dataSpecificConfig)
                print(A2Diagnostics.GetSummaryMessage())
            if self._jobDetails.ClipMinMax:
                # TODO add SD caching - for now re-reading for each image to leave more mem free for A2 above
                # if not self._CACHE_SD:
                if True:
                    with rio.open(self._filePaths.SYNOPTIC_SD_FILE) as f:
                        arr_SD = f.read(1, window=dataReadWindow_rio)
                MinMaxClip(dataImage=arr_Data, flagsImage=arr_Flags, meansImage=arr_Mean, stdImage=arr_SD,
                           flagToCheck=self._flagValues.A2_FILLED, flagToSet=self._flagValues.CLIPPED,
                           floor_ceiling_value=self._dataSpecificConfig.FLOOR_CEILING_ZSCORE,
                           _NDV=self._dataSpecificConfig.NODATA_VALUE,
                           upperHardLimit=self._dataSpecificConfig.CEILING_VALUE,
                           lowerHardLimit=self._dataSpecificConfig.FLOOR_VALUE,
                           nCores=self._jobDetails.NCores)
            day, yr = self.__parseDailyFilename(inputName)
            dataWriteProfile = self.OutputProps.GetRasterioProfile()
            outDataFile = self._outTemplate.format("-FilledData", yr, day)
            outDistFile = self._outTemplate.format("-FilledDists", yr, day)
            outFlagFile = self._outTemplate.format("-FillFlags", yr, day)
            if not os.path.exists(self._filePaths.OUTPUT_FOLDER):
                os.makedirs(self._filePaths.OUTPUT_FOLDER)
            with rio.open(os.path.join(self._filePaths.OUTPUT_FOLDER, outDataFile), 'w',
                          **dataWriteProfile) as dataWriter:
                dataWriter.write(arr_Data, 1)
            with rio.open(os.path.join(self._filePaths.OUTPUT_FOLDER, outDistFile), 'w',
                          **dataWriteProfile) as distWriter:
                distWriter.write(arr_Dists, 1)

            dataWriteProfile.update(dtype=rio_dt.uint8)
            dataWriteProfile.update(nodata=None)
            with rio.open(os.path.join(self._filePaths.OUTPUT_FOLDER, outFlagFile), 'w',
                          **dataWriteProfile) as flagWriter:
                flagWriter.write(arr_Flags, 1)

            #os.remove(dataFile)
            #os.remove(distsFile)
            #os.remove(flagsFile)
        self._intermediateFiles.clear()
