import sys
sys.path.insert(0, r'C:\Users\zool1301.NDPH\Documents\Code_General\MAP-raster-utilities')
from raster_utilities.utils.geotransform_calcs import CalculatePixelLims, CalculateClippedGeoTransform
from raster_utilities.io.TiffFile import SingleBandTiffFile
from gapfill_config import A1SearchConfig, A2SearchConfig, DespeckleSearchConfig, DataSpecificConfig, FlagValues
from gapfill_core_despeckle_and_flag import setSpeckleFlags
from gapfill_core_a1 import a1_core
from gapfill_core_clamp import MinMaxClip3D
import math
import numpy as np
import itertools
import glob
import datetime.date


class gapfiller:
    def __init__(self, fileWildcard, fillForLatLims, fillForLonLims, memLimit=70e9):

        # intialise files

        # initialise limits of fill area in pixel coords of the input files

        # calculate slices to run a1

        # set should-use-slices property

        self.sourceDataHeight = 21600
        self.sourceDataWidth = 43200
        self.sliceXSize = 0
        self.latLims = fillForLatLims
        self.lonLims = fillForLonLims
        self.xLims, self.yLims = CalculatePixelLims(fillForLonLims, fillForLatLims)
        self.inputFileDict = {}
        self.meanFile = ""
        self.stdFile = ""

    def __InitialiseSlices(self, fillShapeZYX, memLimit):


    def RunFill(self):
        self.__InitialiseSlices()
        f

    def InitialiseFiles(self, globPattern, filenameDateParser):
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
            fileDate = filenameDateParser(f)

    def __CalculateSliceSize(self, readShapeZYX, memLimit):
        # readShapeZYX is the dimension of the data we must READ to fill the required output area;
        #  i.e .the fill area plus margins. If we're filling globally it's the same thing.
        dataBPP = 4
        outputsBPP = dataBPP * 2 + 1 # the output data, distances, and flags
        # approximate total number of pixels we can read for each file
        sliceSqrd = memLimit / (readShapeZYX[0] * (dataBPP + outputsBPP))
        # not implementing slicing in y dimension so xSize is total pixels / total height
        sliceXSize = sliceSqrd / readShapeZYX[1]

    def CalculateSearchSlices(self, fillShapeZYX, memLimit):
        '''
        Generate the slice boundaries, slicing only in the x dimension
        (we are using full data in y dimension for now, though this doesn't have to be so)
        The slice boundaries overlap by 2* _SEARCH_RADIUS pixels (_SEARCH_RADIUS on each slice)
        This allows the gap-filling code to run on the non-overlapping section of the slice
        while having all the data necessary to fill a gap up to the edge of that non-overlapping

        :return:
        '''
        # always round up as there's no harm done in passing slightly too much data
        _A1_SEARCH_RADIUS = int(math.sqrt(A1SearchConfig["MAX_NBRS_TO_SEARCH"] / 3.14) + 1)
        _DESPECKLE_SEARCH_RADIUS = int(math.sqrt(DespeckleSearchConfig["MAX_NBRS_TO_SEARCH"] / 3.14) + 1)
        _TOTAL_REQ_MARGIN = _A1_SEARCH_RADIUS + _DESPECKLE_SEARCH_RADIUS

        reqDataL = max(0, self.xLims[0] - _TOTAL_REQ_MARGIN)
        reqDataR = min(self.xLims[1] + _TOTAL_REQ_MARGIN, self.sourceDataWidth)
        reqDataT = max(0, self.yLims[0] - _TOTAL_REQ_MARGIN)
        reqDataB = min(self.yLims[1] + _TOTAL_REQ_MARGIN, self.sourceDataHeight)
        reqDataHeight = reqDataB - reqDataT
        reqDataWidth = reqDataR - reqDataL
        readShape

        # section

        totalFillWidth = fillShapeZYX[2]

        nchunks = (self.totalFillWidth / self.sliceXSize) + 1
        # generate the "chunk" boundaries that will represent the data processed in one thread's job.
        chunkedges = np.linspace(self.xLims[0], self.xLims[1], nchunks + 1).astype(np.int32)
        leftRealEdges = chunkedges[:-1]
        rightRealEdges = chunkedges[1:]
        left_A1_edges = np.clip((chunkedges - _A1_SEARCH_RADIUS)[:-1], 0, np.inf).astype(np.int32)
        right_A1_edges = np.clip((chunkedges + _A1_SEARCH_RADIUS)[1:], -np.inf, self.sourceDataWidth).astype(np.int32)
        left_Despeckle_edges = np.clip((chunkedges - _TOTAL_REQ_MARGIN)[:-1], 0, np.inf).astype(np.int32)
        right_Despeckle_edges = np.clip((chunkedges + _TOTAL_REQ_MARGIN)[1:], -np.inf, self.sourceDataWidth).astype(np.int32)

        # the left and right task boundaries are _SEARCH_RADIUS bigger than the data that will be searched
        # within them so that all pixels can have neighbours, if possible (not at global edge)
        x_offsets_overlapping = zip(left_Despeckle_edges, left_A1_edges, leftRealEdges,
                                    rightRealEdges, right_A1_edges, right_Despeckle_edges)

        # can we pad the top and bottom? (not if we are doing a global run as y-slicing isn't implemented)
        topA1Edge = np.clip(self.yLims[0] - _A1_SEARCH_RADIUS, 0, np.inf).astype(np.int32)
        bottomA1Edge = np.clip(self.yLims[1] + _A1_SEARCH_RADIUS, -np.inf, self.sourceDataHeight).astype(np.int32)
        topDespeckleEdge = np.clip(self.yLims[0] - _TOTAL_REQ_MARGIN, 0, np.inf).astype(np.int32)
        bottomDespeckleEdge = np.clip(self.yLims[1] + _TOTAL_REQ_MARGIN, -np.inf, self.sourceDataHeight).astype(np.int32)

        taskList = list(itertools.product(
            x_offsets_overlapping,
            [(topDespeckleEdge, topA1Edge, self.yLims[0], yLims[1], bottomA1Edge, bottomDespeckleEdge)]))
        return taskList

    def CalculateAvailableMargins(self, sliceInfo):
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
        assert sliceFillHeight == self.totalFillHeight  # not actually supporting y slices for now...

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
        out_SliceFillL = g_sliceFillL - xLims[0]
        out_SliceFillR = out_SliceFillL + sliceFillWidth
        out_SliceFillT = g_sliceFillT - yLims[0]
        out_SliceFillB = out_SliceFillT + sliceFillHeight
        # and as we aren't yet supporting slices in the vertical dimension:
        assert out_SliceFillT == 0

        # The the means / stds we read in have a different coordinate space again: it's the same as the
        # output, but with the overall (total) margins added. These give the coordinates of the output image
        # within the space of the data that is read in
        in_SliceDataL = out_SliceFillL + totalMarginL
        in_SliceDataR = in_SliceDataL + sliceFillWidth
        in_SliceDataT = out_SliceFillT + totalMarginT
        in_SliceDataB = in_SliceDataT + sliceFillHeight

        despeckleMargins = {
            "TOP": sliceDespeckleMarginT,
            "BOTTOM": sliceDespeckleMarginB,
            "LEFT": sliceDespeckleMarginL,
            "RIGHT": sliceDespeckleMarginR
        }

        a1Margins = {
            "TOP": sliceTotalMarginT,
            "BOTTOM": sliceTotalMarginB,
            "LEFT": sliceTotalMarginL,
            "RIGHT": sliceTotalMarginR,
            "FILLHEIGHT": sliceFillHeight,
            "FILLWIDTH": sliceFillWidth
        }

        dataReadWindow = {
            "TOP": g_sliceDespeckleT,
            "BOTTOM": g_sliceDespeckleB,
            "LEFT": g_sliceDespeckleL,
            "RIGHT": g_sliceDespeckleR
        }

        dataWriteWindow = {
            "TOP": out_SliceFillT,
            "BOTTOM": out_SliceFillB,
            "LEFT": out_SliceFillL,
            "RIGHT": out_SliceFillR
        }

        return {
            "despeckleMargins": despeckleMargins,
            "a1Margins": a1Margins,
            "dataReadWindow": dataReadWindow,
            "dataWriteWindow": dataWriteWindow
        }

    def DespeckleAndA1Caller(self, sliceInfo):

        res = self.CalculateAvailableMargins(sliceInfo)
        despeckleMargins = res["despeckleMargins"]
        a1Margins = res["a1Margins"]
        dataReadWindow = res["dataReadWindow"]
        dataWriteWindow = res["dataWriteWindow"]

        sliceDespeckleHeight = dataReadWindow["BOTTOM"] - dataReadWindow["TOP"]
        sliceDespeckleWidth = dataReadWindow["RIGHT"] - dataReadWindow["LEFT"]

        meanReader = SingleBandTiffFile(self.meanFile)
        sliceMeanArr = meanReader.ReadForPixelLims(xLims=(dataReadWindow["LEFT"], dataReadWindow["RIGHT"]),
                                              yLims=(dataReadWindow["TOP"], dataReadWindow["BOTTOM"]))
        sdReader = SingleBandTiffFile(self.stdFile)
        sliceSDArr = sdReader.ReadForPixelLims(xLims=(dataReadWindow["LEFT"], dataReadWindow["RIGHT"]),
                                              yLims=(dataReadWindow["TOP"], dataReadWindow["BOTTOM"]))

        coastReader = SingleBandTiffFile(self.coastFile)
        sliceCoastArr = coastReader.ReadForPixelLims(xLims=(dataReadWindow["LEFT"], dataReadWindow["RIGHT"]),
                                              yLims=(dataReadWindow["TOP"], dataReadWindow["BOTTOM"]))

        for calendarDay, fileList in self.inputFileDict.iteritems():
            sliceDataStack = np.empty(shape=(len(fileList), sliceDespeckleHeight, sliceDespeckleWidth),
                                 dtype='Float32')
            for y in range(len(fileList)):
                f = SingleBandTiffFile(fileList[y])
                sliceDataStack[y] = f.ReadForPixelLims(xLims=(dataReadWindow["LEFT"], dataReadWindow["RIGHT"]),
                                              yLims=(dataReadWindow["TOP"], dataReadWindow["BOTTOM"]))
                assert f.GetNdv() == DataSpecificConfig["NODATA_VALUE"]
            assert sliceDataStack.flags.c_contiguous
            sliceFlagStack = np.zeros(shape=(len(fileList), sliceDespeckleHeight, sliceDespeckleWidth),
                                    dtype=np.uint8)

            dataDict = {
                "Data": sliceDataStack,
                "Flags": sliceFlagStack,
                "Means": sliceMeanArr,
                "Stds": sliceSDArr,
                "LandMask": sliceCoastArr,
                "DistTemplate": None,
                "KnownUnfillable": None
            }

            # Run the despeckle
            despeckleResult = setSpeckleFlags(dataDict, despeckleMargins)

            # Replace the data with the despeckled data, for A1
            dataDict["Data"] = despeckleResult[0]
            dataDict["Flags"] = despeckleResult[1]

            # todo enable a different start position
            a1Result = a1_core(DataStacks=dataDict, Margins=a1Margins, RunFillFromPos=0)
            # result is tuple of (output, dists, new flags, log info)
            del dataDict
            del sliceDataStack
            del sliceFlagStack

            if not ((DataSpecificConfig["DATA_UPPER_LIMIT"] == DataSpecificConfig["NODATA_VALUE"]) and
                DataSpecificConfig["DATA_LOWER_LIMIT"] == DataSpecificConfig["NODATA_VALUE"]):
                sliceMeanArr = np.copy(sliceMeanArr[a1Margins["TOP"]:a1Margins["TOP"]+a1Margins["FILLHEIGHT"],
                                     a1Margins["LEFT"]:a1Margins["LEFT"]+a1Margins["FILLWIDTH"]])
                sliceSDArr = np.copy(sliceSDArr[a1Margins["TOP"]:a1Margins["TOP"] + a1Margins["FILLHEIGHT"],
                                       a1Margins["LEFT"]:a1Margins["LEFT"] + a1Margins["FILLWIDTH"]])
                MinMaxClip3D( dataImage=a1Result[0], flagsImage=a1Result[2],
                              meansImage=sliceMeanArr, stdImage=sliceSDArr,
                              flagToCheck=FlagValues["A1_FILLED"],
                              flagToSet=FlagValues["CLIPPED"],
                              floor_ceiling_value=DataSpecificConfig["FLOOR_CEILING_ZSCORE"],
                              _NDV=DataSpecificConfig["NODATA_VALUE"],
                              upperHardLimit=DataSpecificConfig["DATA_UPPER_LIMIT"],
                              lowerHardLimit=DataSpecificConfig["DATA_LOWER_LIMIT"])
            if len(self.slices) > 1:
                # use intermediate files.
                # in this code we're not saving separate tile files, we're saving global tiffs and just writing
                # part of them at once
                for y in range(len(fileList)):
                    if y < stackFillFromPos:
                        continue
                    inputFilename = os.path.basename(fileList[y])
                    outNameData = os.path.join(tempDir, "A1IntermediateData_" + inputFilename)
                    outNameDist = os.path.join(tempDir, "A1IntermediateDists_" + inputFilename)
                    outNameFlag = os.path.join(tempDir, "A1IntermediateFlags_" + inputFilename)
                    d = SingleBandTiffFile(outNameData)
                    d.SavePart(np.asarray(a1Result[0][y], dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                    d = SingleBandTiffFile(outNameDist)
                    d.SavePart(np.asarray(a1Result[1][y], dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                    d = SingleBandTiffFile(outNameFlag)
                    d.SavePart(np.asarray(a1Result[2][y], dataWriteWindow["LEFT"], dataWriteWindow["TOP"]))
                    d = None

    def A2Caller(self, filename, meanFileObject):

