# hack in the path to map raster utilities for I/O
import sys
sys.path.insert(0, r'C:\Users\zool1301.NDPH\Documents\Code_General\MAP-raster-utilities')
# standard python libraries
import os
import math
import numpy as np
import itertools
import glob
import datetime.date
from collections import defaultdict
# io management functions
from raster_utilities.utils.geotransform_calcs import CalculatePixelLims, CalculateClippedGeoTransform
from raster_utilities.io.TiffFile import SingleBandTiffFile, RasterProps
from raster_utilities.tileProcessor import tileProcessor

# cython functions
from gapfill_core_despeckle_and_flag import setSpeckleFlags
from gapfill_core_a1 import a1_core
from gapfill_core_clamp import MinMaxClip3D
from gapfill_prep_a2 import A2ImageCaller

from gapfill_config import A1SearchConfig, A2SearchConfig, SpeckleSelectionConfig, DataSpecificConfig, FlagItems

from gapfill_utils import PixelMargins, A1DataStack

from gapfill_defaults import DefaultFlagValues, DespeckleThresholdDefaultConfig, \
    A1DefaultParams, A2DefaultParams, DataSpecificDefaultConfig

class gapfiller:
    def __init__(self, fileWildcard, meanFile, sdFile, coastFile,
                 fillForLatLims=None, fillForLonLims=None, memLimit=70e9):

        # initialise input files and index by calendar day and year


        # initialise limits of fill area in pixel coords of the input files

        # calculate slices to run a1

        # set should-use-slices property


        self.latLims = fillForLatLims
        self.lonLims = fillForLonLims
        self.xLims, self.yLims = CalculatePixelLims(fillForLonLims, fillForLatLims)

        self.inputFileDict = defaultdict(defaultdict(list))
        self.InitialiseFiles(fileWildcard)
        self.meanFile = meanFile
        self.stdFile = sdFile
        self.coastFile = coastFile

        if fillForLatLims is None and fillForLonLims is None:
            self.OutputProps = self.InputRasterProps
        else:
            outGT = CalculateClippedGeoTransform(self.InputRasterProps.gt, self.xLims, self.yLims)
            outW = self.xLims[1] - self.xLims[0]
            outH = self.yLims[1] - self.yLims[0]
            outProj = self.InputRasterProps.proj
            outNdv = self.InputRasterProps.ndv
            outRes = self.InputRasterProps.res
            outDT = self.InputRasterProps.datatype
            self.OutputProps = RasterProps(gt=outGT, proj=outProj, ndv=outNdv, width=outW, height=outH,
                                           res=outRes, datatype=outDT)
        self.flagValues = DefaultFlagValues
        self.dataSpecificConfig = DataSpecificDefaultConfig
        self.nCores = 20
        self.a1Config = A1DefaultParams
        self.a2Config = A2DefaultParams
        self.despeckleConfig = DespeckleThresholdDefaultConfig


    def RunFill(self, onlyDays=None, onlyYears=None):
        days = self.inputFileDict.keys()
        if onlyYears is None:
            onlyYears = set()
            for d in self.inputFileDict.keys():
                onlyYears.update(self.inputFileDict[d].keys())
        for calendarDay in days:
            if onlyDays is not None and calendarDay not in onlyDays:
                continue
            for slice in self._Slices:
                self.DespeckleAndA1JobRunner(slice, calendarDay, onlyYears)
            self.A2Caller(onlyYears)

    def InitialiseFiles(self, globPattern):
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
            self.inputFileDict[calendarday][year].append(f)
        initialProps = SingleBandTiffFile(allFiles[0])
        self.InputRasterProps = initialProps.GetExistingProperties()

        base = os.path.basename(f)

    def __parseDailyFilename(self, f):
        yr = base[1:5]
        day = base[5:8]
        return(day,yr)

    def __CalculateSliceSize(self, readShapeZYX, memLimit):
        # readShapeZYX is the dimension of the data we must READ to fill the required output area;
        #  i.e .the fill area plus margins. If we're filling globally it's the same thing.
        dataBPP = 4
        outputsBPP = dataBPP * 2 + 1 # the output data, distances, and flags
        # approximate total number of pixels we can read for each file
        sliceSqrd = memLimit / (readShapeZYX[0] * (dataBPP + outputsBPP))
        # not implementing slicing in y dimension so xSize is total pixels / total height
        sliceXSize = sliceSqrd / readShapeZYX[1]
        return sliceXSize

    def CalculateJobLimits(self, fillShapeZYX, memLimit):
        '''
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
        '''
        # always round up as there's no harm done in passing slightly too much data
        _A1_SEARCH_RADIUS = int(math.sqrt(self.a1Config.SPIRAL.MAX_NBRS_TO_SEARCH / 3.14) + 1)
        _DESPECKLE_SEARCH_RADIUS = int(math.sqrt(self.despeckleConfig.SPIRAL.MAX_NBRS_TO_SEARCH / 3.14) + 1)
        _TOTAL_DESIRED_MARGIN = _A1_SEARCH_RADIUS + _DESPECKLE_SEARCH_RADIUS

        reqDataL = max(0, self.xLims[0] - _TOTAL_DESIRED_MARGIN)
        reqDataR = min(self.xLims[1] + _TOTAL_DESIRED_MARGIN, self.InputRasterProps.width)
        reqDataT = max(0, self.yLims[0] - _TOTAL_DESIRED_MARGIN)
        reqDataB = min(self.yLims[1] + _TOTAL_DESIRED_MARGIN, self.InputRasterProps.height)
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
        left_A1_edges = np.clip((chunkEdges - _A1_SEARCH_RADIUS)[:-1], 0, np.inf).astype(np.int32)
        right_A1_edges = np.clip((chunkEdges + _A1_SEARCH_RADIUS)[1:], -np.inf, self.InputRasterProps.width).astype(np.int32)
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

        taskList = list(itertools.product(
            x_offsets_overlapping,
            [(topDespeckleEdge, topA1Edge, self.yLims[0], self.yLims[1], bottomA1Edge, bottomDespeckleEdge)]))
        return taskList

    def CalculateAvailableJobMargins(self, sliceInfo):
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

    def DespeckleAndA1JobRunner(self, sliceInfo, calendarDay, years):

        margins = self.CalculateAvailableJobMargins(sliceInfo)
        despeckleMargins = margins["despeckleMargins"]
        a1Margins = margins["a1Margins"]
        dataReadWindow = margins["dataReadWindow"]
        dataWriteWindow = margins["dataWriteWindow"]

        sliceGT = CalculateClippedGeoTransform(self.OutputProps.gt,
                                               xPixelLims=(dataWriteWindow.Left, dataWriteWindow.Right),
                                               yPixelLims=(dataWriteWindow.Top, dataWriteWindow.Bottom))
        sliceDespeckleHeight = dataReadWindow.Bottom - dataReadWindow.Top
        sliceDespeckleWidth = dataReadWindow.Right - dataReadWindow.Left

        meanReader = SingleBandTiffFile(self.meanFile)
        sliceMeanArr , _, _, _ = meanReader.ReadForPixelLims(xLims=(dataReadWindow.Left, dataReadWindow.Right),
                                              yLims=(dataReadWindow.Top, dataReadWindow.Bottom))
        sdReader = SingleBandTiffFile(self.stdFile)
        sliceSDArr , _, _, _ = sdReader.ReadForPixelLims(xLims=(dataReadWindow.Left, dataReadWindow.Right),
                                              yLims=(dataReadWindow.Top, dataReadWindow.Bottom))

        coastReader = SingleBandTiffFile(self.coastFile)
        sliceCoastArr , _, _, _ = coastReader.ReadForPixelLims(xLims=(dataReadWindow.Left, dataReadWindow.Right),
                                              yLims=(dataReadWindow.Top, dataReadWindow.Bottom))

        for calendarDay, fileList in self.inputFileDict.iteritems():

            dataStack = A1DataStack()
            dataStack.DataArray3D = np.empty(shape=(len(fileList), sliceDespeckleHeight, sliceDespeckleWidth),
                                 dtype='Float32')
            dataStack.FlagsArray3D = np.zeros(shape=(len(fileList), sliceDespeckleHeight, sliceDespeckleWidth),
                                    dtype=np.uint8)
            dataStack.MeansArray2D = sliceMeanArr
            dataStack.SDArray2D = sliceSDArr
            dataStack.Coastline2D = sliceCoastArr
            for y in range(len(fileList)):
                f = SingleBandTiffFile(fileList[y])
                dataStack.DataArray3D[y], _, _, _ = f.ReadForPixelLims(xLims=(dataReadWindow.Left, dataReadWindow.Right),
                                                       yLims=(dataReadWindow.Top, dataReadWindow.Bottom))
                assert f.GetNdv() == DataSpecificConfig.NODATA_VALUE
            assert dataStack.DataArray3D.flags.c_contiguous

            #data.DistanceTemplate3D = None
            #data.KnownUnfillable2D = None

            # Run the despeckle
            despeckleResult = setSpeckleFlags(dataStacks = dataStack, margins = despeckleMargins,
                                              flagValues=self.flagValues, dataConfig=self.dataSpecificConfig,
                                              speckleConfig=self.despeckleConfig, nCores=self.nCores)

            # Replace the data with the despeckled data, for A1
            dataStack.DataArray3D = despeckleResult[0]
            dataStack.FlagsArray3D = despeckleResult[1]

            # todo enable a different start position
            # todo replace data  inplace on dataStacks and just return the diagnostics
            a1Result = a1_core(dataStacks=dataStack, margins=a1Margins,
                               flagValues=self.flagValues, dataConfig=self.dataSpecificConfig,
                               a1Config = self.a1Config, nCores=self.nCores)
            # result is tuple of (output, dists, new flags, log info)
            del dataStack

            if not ((self.dataSpecificConfig.CEILING_VALUE == self.dataSpecificConfig.NODATA_VALUE) and
                self.dataSpecificConfig.DATA_LOWER_LIMIT == self.dataSpecificConfig.NODATA_VALUE):
                # apply clip if either the floor / ceiling values are set
                #TODO add the heights
                sliceMeanArr = np.copy(sliceMeanArr[a1Margins.Top:a1Margins.Top+a1Margins["FILLHEIGHT"],
                                     a1Margins.Left:a1Margins.Left+a1Margins["FILLWIDTH"]])
                sliceSDArr = np.copy(sliceSDArr[a1Margins.Top:a1Margins.Top + a1Margins["FILLHEIGHT"],
                                       a1Margins.Left:a1Margins.Left + a1Margins["FILLWIDTH"]])
                MinMaxClip3D(dataStacks=dataStack,
                             flagToCheck=self.flagValues.A1_FILLED,
                             flagToSet=self.flagValues.CLIPPED,
                             dataConfig=self.dataSpecificConfig,
                             nCores=self.nCores)


            if len(self.slices) > 1:
                # use intermediate files.
                # in this code we're not saving separate tile files, we're saving global tiffs and just writing
                # part of them at once
                for y in range(len(fileList)):
                    if y < dataStack.FillFromZPosition:
                        continue
                    inputFilename = os.path.basename(fileList[y])
                    outNameData = os.path.join(self.tempDir, "A1IntermediateData_" + inputFilename)
                    outNameDist = os.path.join(self.tempDir, "A1IntermediateDists_" + inputFilename)
                    outNameFlag = os.path.join(self.tempDir, "A1IntermediateFlags_" + inputFilename)
                    d = SingleBandTiffFile(outNameData)
                    d.SavePart(np.asarray(a1Result[0][y], dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                    d = SingleBandTiffFile(outNameDist)
                    d.SavePart(np.asarray(a1Result[1][y], dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                    d = SingleBandTiffFile(outNameFlag)
                    d.SavePart(np.asarray(a1Result[2][y], dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                    d = None

    def A2BatchRunner(self):
        for img in self.GetImages():
            A2ImageCaller(img, fla)
