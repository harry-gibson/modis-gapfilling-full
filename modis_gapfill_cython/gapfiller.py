# hack in the path to map raster utilities for I/O
import sys

sys.path.insert(0, r'E:\Data\Harry\Documents\Git\MAP-raster-utilities')
# standard python libraries
import os
import math
import numpy as np
import itertools
import glob
from collections import defaultdict
# io management functions
from raster_utilities.utils.geotransform_calcs import CalculatePixelLims, CalculateClippedGeoTransform
from raster_utilities.io.TiffFile import SingleBandTiffFile, RasterProps
# from raster_utilities.tileProcessor import tileProcessor

# cython functions
from .gapfill_core_despeckle_and_flag import setSpeckleFlags
from .gapfill_core_a1 import a1_core
from .gapfill_core_clamp import MinMaxClip3D, MinMaxClip
from .gapfill_prep_a2 import A2ImageCaller

# configuration types (NamedTuples with yaml-parsing factory methods)
from .gapfill_config_types import GapfillFilePaths, GapfillJobConfig, DataLimitsConfig, \
    DespeckleConfig, A1SearchConfig, A2SearchConfig, FlagItems

from .gapfill_utils import PixelMargins, A1DataStack, A2PassData


class GapFiller:
    def __init__(self, gapfillFilePaths: GapfillFilePaths,
                 despeckleconfig: DespeckleConfig,
                 a1Config: A1SearchConfig,
                 a2Config: A2SearchConfig,
                 dataSpecificConfig: DataLimitsConfig,
                 flagValues: FlagItems,
                 jobDetails: GapfillJobConfig
                 ):
        self.nCores = 20

        assert isinstance(gapfillFilePaths, GapfillFilePaths)
        assert isinstance(despeckleconfig, DespeckleConfig)
        assert isinstance(a1Config, A1SearchConfig)
        assert isinstance(a2Config, A2SearchConfig)
        assert isinstance(dataSpecificConfig, DataLimitsConfig)
        assert isinstance(flagValues, FlagItems)
        assert isinstance(jobDetails, GapfillJobConfig)
        self._despeckleConfig = despeckleconfig
        self._a1Config = a1Config
        self._a2Config = a2Config
        self._dataSpecificConfig = dataSpecificConfig
        self._flagValues = flagValues
        self._jobDetails = jobDetails

        # initialise input files and index by calendar day and year
        self._inputFileDict = defaultdict(defaultdict(str))
        self.__InitialiseFiles(gapfillFilePaths.DATA_FILES_GLOB_PATTERN)
        self._filePaths = gapfillFilePaths
        self._outTemplate = "FILLED-OUTPUT{}.{}.{}.TemporalSummary.Res.SpatialSummary.tif"
        self._intermediateFiles = defaultdict(dict)

        # initialise limits of fill area in pixel coords of the input files
        _latLims = (jobDetails.XMax_Deg, jobDetails.XMin_Deg)
        _lonLims = (jobDetails.YMin_Deg, jobDetails.YMax_Deg)

        if None in _latLims and None in _lonLims:
            self.OutputProps = self.InputRasterProps
            self.xLims = (0, self.InputRasterProps.width)
            self.yLims = (0, self.InputRasterProps.height)

        else:
            self.xLims, self.yLims = CalculatePixelLims(self.InputRasterProps.gt, _lonLims, _latLims)
            outGT = CalculateClippedGeoTransform(self.InputRasterProps.gt, self.xLims, self.yLims)
            outW = self.xLims[1] - self.xLims[0]
            outH = self.yLims[1] - self.yLims[0]
            outProj = self.InputRasterProps.proj
            outNdv = self.InputRasterProps.ndv
            outRes = self.InputRasterProps.res
            outDT = self.InputRasterProps.datatype
            self.OutputProps = RasterProps(gt=outGT, proj=outProj, ndv=outNdv, width=outW, height=outH,
                                           res=outRes, datatype=outDT)

        # calculate slices to run a1
        self._slices = self.CalculateSliceEdges()

    def RunFill(self):
        allDays = self._inputFileDict.keys()
        startYear = self._jobDetails.StartYear
        if startYear is None:
            startYear = 0  # just some value smaller than the first year
        onlyForDays = self._jobDetails.CalendarDaysToFill
        for calendarDay in allDays:
            if onlyForDays is not None and calendarDay not in onlyForDays:
                # if we're not filling for this calendar day, we don't need to do anything
                continue
            for slice in self._slices:
                # years are slightly more complex as we always need to read the day's data for all years
                # this populates the self._intermediateFiles dictionary
                self.DespeckleAndA1SliceRunner(slice, calendarDay, startYear)
            self.A2BatchRunner(startYear)

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
            year, calendarday = self.__parseDailyFilename(f)
            self._inputFileDict[calendarday][year] = f
        initialProps = SingleBandTiffFile(allFiles[0])
        self._InputRasterProps = initialProps.GetExistingProperties()

    def __parseDailyFilename(self, f):
        '''Get julian day and year from a filename in old or new MAP MODIS filename formats,
        and set instance variable _outTemplate to match the non-temporally-varying parts of this
        input file (or, on subsequent calls, check these are the same as last time).
        _outTemplate will be set to e.g.
            LST_Day{}.{}.{}.Data.1km.Data'''
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

    @staticmethod
    def __CalculateSliceSize(readShapeZYX, memLimit):
        '''Calculate the X-size of the slices we can run for the given Y, Z dimensions and mem limit'''
        # readShapeZYX is the dimension of the data we must READ to fill the required output area;
        #  i.e .the fill area plus margins. If we're filling globally it's the same thing.
        dataBPP = 4
        outputsBPP = dataBPP * 2 + 1  # the output data, distances, and flags
        # approximate total number of pixels we can read for each file
        sliceSqrd = memLimit / (readShapeZYX[0] * (dataBPP + outputsBPP))
        # not implementing slicing in y dimension so xSize is total pixels / total height
        sliceXSize = sliceSqrd / readShapeZYX[1]
        return sliceXSize

    def CalculateSliceEdges(self, fillShapeZYX, memLimit):
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
        reqDataR = min(self.xLims[1] + _TOTAL_DESIRED_MARGIN, self._InputRasterProps.width)
        reqDataT = max(0, self.yLims[0] - _TOTAL_DESIRED_MARGIN)
        reqDataB = min(self.yLims[1] + _TOTAL_DESIRED_MARGIN, self._InputRasterProps.height)
        reqDataHeight = reqDataB - reqDataT
        reqDataWidth = reqDataR - reqDataL

        # get an estimated maximum width for the slices based on the available memory and the number of years (vertical
        # size of the stack) and height (y dimension) of the data (not slicing in Y dimension)
        readShapeZYX = (fillShapeZYX[0], reqDataHeight, reqDataWidth)
        sliceXSize = self.__CalculateSliceSize(readShapeZYX, memLimit)

        # how many slices are we going to need
        totalFillWidth = fillShapeZYX[2]
        nChunks = (totalFillWidth / sliceXSize) + 1

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
                                 -np.inf, self._InputRasterProps.width).astype(np.int32)
        # the boundary of these slices plus the margin of extra data for A1 and despeckle, but respecting the fact
        # we can't go beyond the edge of the source data
        left_Despeckle_edges = np.clip((chunkEdges - _TOTAL_DESIRED_MARGIN)[:-1],
                                       0, np.inf).astype(np.int32)
        right_Despeckle_edges = np.clip((chunkEdges + _TOTAL_DESIRED_MARGIN)[1:],
                                        -np.inf, self._InputRasterProps.width).astype(np.int32)

        # the left and right task boundaries are _SEARCH_RADIUS bigger than the data that will be searched
        # within them so that all pixels can have neighbours, if possible (not at global edge)
        x_offsets_overlapping = zip(left_Despeckle_edges, left_A1_edges, leftRealEdges,
                                    rightRealEdges, right_A1_edges, right_Despeckle_edges)

        # can we pad the top and bottom? (not if we are doing a global run as y-slicing isn't implemented)
        topA1Edge = np.clip(self.yLims[0] - _A1_SEARCH_RADIUS, 0, np.inf).astype(np.int32)
        bottomA1Edge = np.clip(self.yLims[1] + _A1_SEARCH_RADIUS,
                               -np.inf, self._InputRasterProps.height).astype(np.int32)
        topDespeckleEdge = np.clip(self.yLims[0] - _TOTAL_DESIRED_MARGIN, 0, np.inf).astype(np.int32)
        bottomDespeckleEdge = np.clip(self.yLims[1] + _TOTAL_DESIRED_MARGIN,
                                      -np.inf, self._InputRasterProps.height).astype(np.int32)

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

        a1Margins = {
            "TOP": sliceTotalMarginT,
            "BOTTOM": sliceTotalMarginB,
            "LEFT": sliceTotalMarginL,
            "RIGHT": sliceTotalMarginR,
            "FILLHEIGHT": sliceFillHeight,
            "FILLWIDTH": sliceFillWidth
        }

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
        dataReadWindow = margins["dataReadWindow"]
        dataWriteWindow = margins["dataWriteWindow"]
        if not calendarDay in self._inputFileDict:
            raise FileNotFoundError("No data files were identified for the requested calendar day ({})"
                                    .format(calendarDay))
        sliceGT = CalculateClippedGeoTransform(self.OutputProps.gt,
                                               xPixelLims=(dataWriteWindow.Left, dataWriteWindow.Right),
                                               yPixelLims=(dataWriteWindow.Top, dataWriteWindow.Bottom))
        sliceDespeckleHeight = dataReadWindow.Bottom - dataReadWindow.Top
        sliceDespeckleWidth = dataReadWindow.Right - dataReadWindow.Left

        meanReader = SingleBandTiffFile(self._filePaths.SYNOPTIC_MEAN_FILE)
        sliceMeanArr, _, _, _ = meanReader.ReadForPixelLims(xLims=(dataReadWindow.Left, dataReadWindow.Right),
                                                            yLims=(dataReadWindow.Top, dataReadWindow.Bottom))
        sdReader = SingleBandTiffFile(self._filePaths.SYNOPTIC_SD_FILE)
        sliceSDArr, _, _, _ = sdReader.ReadForPixelLims(xLims=(dataReadWindow.Left, dataReadWindow.Right),
                                                        yLims=(dataReadWindow.Top, dataReadWindow.Bottom))

        coastReader = SingleBandTiffFile(self._filePaths.COASTLINE_FILE)
        sliceCoastArr, _, _, _ = coastReader.ReadForPixelLims(xLims=(dataReadWindow.Left, dataReadWindow.Right),
                                                              yLims=(dataReadWindow.Top, dataReadWindow.Bottom))
        yearFileDictThisDay = self._inputFileDict[calendarDay]
        sortedFiles = [f for (y, f) in sorted(yearFileDictThisDay)]
        if startYear != 0:
            try:
                startFillFromPos = sorted(yearFileDictThisDay.keys()).index(startYear)
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
        for zPos in range(len(sortedFiles)):
            f = SingleBandTiffFile(sortedFiles[zPos])
            dataStack.DataArray3D[zPos], _, _, _ = f.ReadForPixelLims(
                xLims=(dataReadWindow.Left, dataReadWindow.Right),
                yLims=(dataReadWindow.Top, dataReadWindow.Bottom))
            assert f.GetNdv() == self._dataSpecificConfig.NODATA_VALUE
        assert dataStack.DataArray3D.flags.c_contiguous

        # Run the despeckle
        # dataStack.DataArray3D and dataStack.FlagsArray3D are modified in place. They will have a different
        # (smaller) shape as the speckle search margins will be removed.
        # the returned object is an instance of DespeckleDiagnostics
        # Note that despeckle runs on the whole stack, not taking account of firstYearToFill, this is
        # because we don't want extreme values to remain to be usable in deriving fill values
        despeckleResult = setSpeckleFlags(dataStacks=dataStack, margins=despeckleMargins,
                                          flagValues=self._flagValues, dataConfig=self._dataSpecificConfig,
                                          speckleConfig=self._despeckleConfig, nCores=self.nCores)

        # Run A1
        # dataStack.DataArray3D, dataStack.FlagsArray3D, and dataStack.DistanceTemplate3D are modified in place.
        # They will have a different (smaller) shape as the A1 search margins will be removed.
        # The returned object is an instance of A1Diagnostics
        a1Result = a1_core(dataStacks=dataStack, margins=a1Margins,
                           flagValues=self._flagValues, dataConfig=self._dataSpecificConfig,
                           a1Config=self._a1Config,
                           nCores=self.nCores)
        # dataStack has been modified in place: DataArray32, FlagsArray3D and DistanceTemplate3D members
        # now contain the results

        if not ((self._dataSpecificConfig.CEILING_VALUE == self._dataSpecificConfig.NODATA_VALUE) and
                self._dataSpecificConfig.FLOOR_VALUE == self._dataSpecificConfig.NODATA_VALUE):
            # apply clip if either the floor / ceiling values are set
            # TODO add the heights
            sliceMeanArr = np.copy(sliceMeanArr[a1Margins.Top:a1Margins.Top + a1Margins["FILLHEIGHT"],
                                   a1Margins.Left:a1Margins.Left + a1Margins["FILLWIDTH"]])
            sliceSDArr = np.copy(sliceSDArr[a1Margins.Top:a1Margins.Top + a1Margins["FILLHEIGHT"],
                                 a1Margins.Left:a1Margins.Left + a1Margins["FILLWIDTH"]])
            # TODO set this to respect firstYearToFill
            MinMaxClip3D(dataStacks=dataStack,
                         flagToCheck=self._flagValues.A1_FILLED,
                         flagToSet=self._flagValues.CLIPPED,
                         dataConfig=self._dataSpecificConfig,
                         nCores=self.nCores)

        # TODO fix this
        if True:  # len(self._slices) > 1:
            # use intermediate files.
            # in this code we're not saving separate tile files, we're saving global tiffs and just writing
            # part of them at once
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
                thisOutFile = SingleBandTiffFile(outNameData)
                try:
                    thisOutFile.SetProperties(self.OutputProps)
                    self._intermediateFiles[inputFilename]["Data"] = outNameData
                except RuntimeError:
                    pass
                thisOutFile.SavePart(np.asarray(dataStack.DataArray3D[y],
                                                dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                thisOutFile = SingleBandTiffFile(outNameDist)
                try:
                    thisOutFile.SetProperties(self.OutputProps)
                    self._intermediateFiles[inputFilename]["Distances"] = outNameDist
                except RuntimeError:
                    pass
                thisOutFile.SavePart(np.asarray(dataStack.DistanceTemplate3D[y],
                                                dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                thisOutFile = SingleBandTiffFile(outNameFlag)
                try:
                    op = self.OutputProps
                    flagProps = RasterProps(gt=op.gt, proj=op.proj, ndv=None,
                                            width=op.width, height=op.height, res=op.res,
                                            datatype=gdal.GDT_BYTE)
                    thisOutFile.SetProperties(flagProps)
                    self._intermediateFiles[inputFilename]["Flags"] = outNameFlag
                except RuntimeError:
                    pass
                thisOutFile.SavePart(np.asarray(dataStack.FlagsArray3D[y],
                                                dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                thisOutFile = None

    def A2BatchRunner(self):
        # cf. a2Caller, but run for a provided filename

        arr_Mean, _, _, _ = SingleBandTiffFile(self._filePaths.SYNOPTIC_MEAN_FILE).ReadAll()
        if self._CACHE_SD:
            arr_SD, _, _, _ = SingleBandTiffFile(self._filePaths.SYNOPTIC_SD_FILE).ReadAll()
        for inputName, intermediateDataForDate in self._intermediateFiles:
            dataFile = intermediateDataForDate["Data"]
            dataReader = SingleBandTiffFile(dataFile)
            arr_Data, gt_Data, proj_Data, ndv_Data = dataReader.ReadAll()
            props_Data = dataReader.GetProperties()

            distsFile = intermediateDataForDate["Distances"]
            distsReader = SingleBandTiffFile(distsFile)
            arr_Dists, gt_Dists, proj_Dists, ndv_Dists = distsReader.ReadAll()
            props_Dists = distsReader.GetProperties()

            flagsFile = intermediateDataForDate["Flags"]
            flagsReader = SingleBandTiffFile(flagsFile)
            arr_Flags, gt_Flags, proj_Flags, ndv_Flags = flagsReader.ReadAll()
            props_Flags = flagsReader.GetProperties()

            if self._jobDetails.RunA2:
                A2ImageCaller(dataImageIn=arr_Data, flagsImageIn=arr_Flags, distImageIn=arr_Dists, meanImageIn=arr_Mean,
                              a2Config=self._a2Config, flagValues=self._flagValues,
                              dataConfig=self._dataSpecificConfig)
            if self._jobDetails.ClipMinMax:
                #if not self._CACHE_SD: # TODO add this
                if True:
                    arr_SD, _, _, _ = SingleBandTiffFile(self._filePaths.SYNOPTIC_SD_FILE).ReadAll()
                MinMaxClip(dataImage=arr_Data, flagsImage=arr_Flags, meansImage=arr_Mean, stdImage=arr_SD,
                           flagToCheck=self._flagValues.A2_FILLED, flagToSet=self._flagValues.CLIPPED,
                           floor_ceiling_value=self._dataSpecificConfig.FLOOR_CEILING_ZSCORE,
                           _NDV=self._dataSpecificConfig.NODATA_VALUE,
                           upperHardLimit=self._dataSpecificConfig.CEILING_VALUE,
                           lowerHardLimit=self._dataSpecificConfig.FLOOR_VALUE,
                           nCores=self.nCores)
            day, yr = self.__parseDailyFilename(inputName)
            outDataFile = self._outTemplate.format("-FilledData", yr, day)
            outDistFile = self._outTemplate.format("-FilledDists", yr, day)
            outFlagFile = self._outTemplate.format("-FillFlags", yr, day)
            dataWriter = SingleBandTiffFile(os.path.join(self._filePaths.OUTPUT_FOLDER, outDataFile))
            dataWriter.SetProperties(props_Data)
            dataWriter.Save(arr_Data)
            distWriter = SingleBandTiffFile(os.path.join(self._filePaths.OUTPUT_FOLDER, outDistFile))
            distWriter.SetProperties(props_Dists)
            distWriter.Save(arr_Dists)
            flagWriter = SingleBandTiffFile(os.path.join(self._filePaths.OUTPUT_FOLDER, outFlagFile))
            flagWriter.SetProperties(props_Flags)
            flagWriter.Save(arr_Flags)

            os.remove(dataFile)
            os.remove(distsFile)
            os.remove(flagsFile)
            self._intermediateFiles.clear()
