from typing import NamedTuple

# The flags output is an 8 bit raster which represents 8 separate flag conditions as a bitmask, defined here

#FlagItems = namedtuple("FlagItems", ["OCEAN", "FAILURE", "EXTREME", "SPECKLE", "A1_FILLED", "A1_FULL",
#                                       "A2_FILLED", "CLIPPED"])
class FlagItems(NamedTuple):
    OCEAN: int
    FAILURE: int
    EXTREME: int
    SPECKLE: int
    A1_FILLED: int
    A1_FULL: int
    A2_FILLED: int
    CLIPPED: int

    @classmethod
    def from_yaml_config(cls, runParams):
        fP = runParams['Flags']
        return cls(int(fP['ocean']), int(fP['failure']), int(fP['extreme']), int(fP['speckle']),
                   int(fP['a1_filled']), int(fP['a1_full']), int(fP['a2_filled']), int(fP['clipped']))

class GapfillJobConfig(NamedTuple):
    XMin: float
    YMin: float
    XMax: float
    YMax: float
    CalendarDaysToFill: list
    StartYear: int
    ClipMinMax: bool
    RunA2: bool

    @classmethod
    def from_yaml_config(cls, runConfig):
        jC = runConfig['FillJob']
        jobConfigParsed = cls(
            XMin=float(jC['XMin']),
            XMax=float(jC['XMax']),
            YMin=float(jC['YMin']),
            YMax=float(jC['YMax']),
            CalendarDaysToFill=float(jC['CalendarDaysToFill']),
            StartYear=float(jC['StartYear']),
            ClipMinMax=float(jC['ClipMinMax']),
            RunA2=float(jC['RunA2'])
        )
        return jobConfigParsed


class SpiralSearchConfig(NamedTuple):
    # Number of cells to search; the radius in pixel distance terms is approx sqrt(value/pi)
    MAX_NBRS_TO_SEARCH: int
    # Search is successful if we find at least this many within the max number
    MIN_REQUIRED_NBRS: int
    # Stop searching early after this many even if we haven't gone to the full radius
    MAX_USED_NBRS: int

    @classmethod
    def from_yaml_config(cls, singleFillConfig):
        sC = singleFillConfig['spiral']
        spiralParsed = cls(
            MAX_NBRS_TO_SEARCH=float(sC['max_nbrs_to_search']),
            MIN_REQUIRED_NBRS=float(sC['min_required_nbrs']),
            MAX_USED_NBRS=float(sC['max_used_nbrs'])
        )
        return spiralParsed


class DespeckleConfig(NamedTuple):
    SPIRAL_CONFIG: SpiralSearchConfig
    EXTREME_BEYOND_SD: float
    SPECKLE_BEYOND_SD: float
    SPECKLE_NBR_Z_THRESH: float

    @classmethod
    def from_yaml_config(cls, runParams):
        dC = runParams['Despeckle']
        spiralParsed = SpiralSearchConfig.from_yaml_config(dC)
        configParsed = cls(
            SPIRAL_CONFIG=spiralParsed,
            EXTREME_BEYOND_SD=float(dC['extreme_sd_threshold']),
            SPECKLE_BEYOND_SD=float(dC['speckle_sd_threshold']),
            SPECKLE_NBR_Z_THRESH=float(dC['speckle_neighbour_zscore_threshold'])
        )
        return configParsed


class A1SearchConfig(NamedTuple):
    SPIRAL_CONFIG: SpiralSearchConfig
    # Generate fill values from neighbour values by comparing "RATIO" or "DIFFERENCE" between them?
    FILL_GENERATION_METHOD: str
    # If the fill generation method is "RATIO" then what should be the maximum allowable ratio, to allow for
    # near-zero divisors?
    MAX_ALLOWABLE_RATIO: float
    # If true then the highest and lowest single partial fill values from neighbours will be dropped as a
    # further level of protection against outliers
    TRIM_FULL_OUTLIERS: bool

    @classmethod
    def from_yaml_config(cls, runParams):
        a1C = runParams['A1']
        spiralParsed = SpiralSearchConfig.from_yaml_config(a1C)
        configParsed = cls(
            SPIRAL_CONFIG=spiralParsed,
            FILL_GENERATION_METHOD="RATIO" if bool(a1C['use_ratio_fill']) else "DIFFERENCE",
            MAX_ALLOWABLE_RATIO=float(a1C['max_allowable_ratio']),
            TRIM_FULL_OUTLIERS=bool(a1C['trim_full_outliers'])
        )
        return configParsed


class A2SearchConfig(NamedTuple):
    # n neighbours should normally just be 8 for A2 as the search uses previously-generated
    # values on each step, i.e. "smears" values across in the direction of the search
    SPIRAL_CONFIG: SpiralSearchConfig
    # Generate fill values from neighbour values by comparing "RATIO" or "DIFFERENCE" between them?
    FILL_GENERATION_METHOD: str
    # If the fill generation method is "RATIO" then what should be the maximum allowable ratio, to allow for
    # near-zero divisors?
    MAX_ALLOWABLE_RATIO: float
    # If true then the highest and lowest single partial fill values from neighbours will be dropped as a
    # further level of protection against outliers
    PASS_AVERAGE_TYPE: str

    @classmethod
    def from_yaml_config(cls, runParams):
        a2C = runParams['A2']
        spiralParsed = SpiralSearchConfig.from_yaml_config(a2C)
        configParsed = cls(
            SPIRAL_CONFIG=spiralParsed,
            FILL_GENERATION_METHOD="RATIO" if bool(a2C['use_ratio_fill']) else "DIFFERENCE",
            MAX_ALLOWABLE_RATIO=float(a2C['max_allowable_ratio']),
            PASS_AVERAGE_TYPE=a1C['pass_averaging_method']
        )
        return configParsed


DataSpecificConfig = namedtuple("DataSpecificConfig", [
    "CEILING_VALUE", # Hard upper limit to clip fill values to. Use NODATA_VALUE to not clip in this direction
    "FLOOR_VALUE", # Hard lower limit to clip fill values to. Use NODATA_VALUE to not clip in this direction
    "FLOOR_CEILING_ZSCORE", # Clip fill values more than this number of SD from the mean, only if ceiling/floor also set
    "CORRECTION_OFFSET", # Value to add to all data before further processing, e.g. 0.15 if celsius files were made a
    # bit wrong with 273 degree offset instead of 271.15 (totally random example, asking for a friend)
    "NODATA_VALUE", # Should generally be the same as whatever nodata is stored as in the tiff files
    "ABSOLUTE_ZERO_OFFSET" # Offset to convert data into an absolute scale if appropriate, e.g. -273.15 to convert celsius
    # into kelvin
])


DespeckleDiagnostics = namedtuple("DespeckleDiagnostics", [
    "SpeckleCellCount",
    "ExtremeCellCount",
    "GoodCellCount",
    "OceanCellCount",
    "ClearedSpeckleCellCount",
    "TimeSeconds"
])

A1Diagnostics = namedtuple("A1Diagnostics", [
    "TotalCellsSearched",
    "CellsWithGoodData",
    "OceanCells",
    "NeverDataLocations",
    "GapCellsTotal",
    "GapCellsTooBig",
    "PermanentGapCells",
    "GapCellsFullyFilled",
    "GapCellsPartFilled",
    "FailedInsufficientPairs",
    "FailedNoPairs",
    "TotalAlternateYearsScanned",
    "TotalNbrsChecked",
    "TotalNbrsUsed",
    "TimeSeconds"
])

A2Diagnostics = namedtuple("A2Diagnostics", [
    "GapCellsTotal",
    "GapCellsFilled",
    "TimeSeconds"
])

GapfillFilePaths = namedtuple("GapfillFilePaths", [
    "DATA_FILES_GLOB_PATTERN",
    "SYNOPTIC_MEAN_FILE",
    "SYNOPTIC_SD_FILE",
    "COASTLINE_FILE",
    "OUTPUT_FOLDER",
    "TEMP_FOLDER"
])

