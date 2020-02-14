from collections import namedtuple
# The flags output is an 8 bit raster which represents 8 separate flag conditions as a bitmask, defined here

FlagItems = namedtuple("FlagItems", ["OCEAN", "FAILURE", "EXTREME", "SPECKLE", "A1_FILLED", "A1_FULL",
                                       "A2_FILLED", "CLIPPED"])

SpiralSearchConfig = namedtuple("SpiralSearchConfig", [
    # Number of cells to search; the radius in pixel distance terms is approx sqrt(value/pi)
    "MAX_NBRS_TO_SEARCH",
    # Search is successful if we find at least this many within the max number
    "MIN_REQUIRED_NBRS",
    # Stop searching early after this many even if we haven't gone to the full radius
    "MAX_USED_NBRS"
])

DespeckleConfig = namedtuple("DespeckleConfig", [
    "SPIRAL_CONFIG",
    "EXTREME_BEYOND_SD",
    "SPECKLE_BEYOND_SD",
    "SPECKLE_NBR_Z_THRESH"
])

DataCharacteristicsConfig = namedtuple("DataSpecificConfig", [
    "CEILING_VALUE", # Hard upper limit to clip fill values to. Use NODATA_VALUE to not clip in this direction
    "FLOOR_VALUE", # Hard lower limit to clip fill values to. Use NODATA_VALUE to not clip in this direction
    "FLOOR_CEILING_ZSCORE", # Clip fill values more than this number of SD from the mean, only if ceiling/floor also set
    "CORRECTION_OFFSET", # Value to add to all data before further processing, e.g. 0.15 if celsius files were made a
    # bit wrong with 273 degree offset instead of 271.15 (totally random example, asking for a friend)
    "NODATA_VALUE", # Should generally be the same as whatever nodata is stored as in the tiff files
    "ABSOLUTE_ZERO_OFFSET" # Offset to convert data into an absolute scale if appropriate, e.g. -273.15 to convert celsius
    # into kelvin
])

A1SearchConfig = namedtuple("A1SearchConfig", [
    # an instance of SpiralSearchConfig
    "SPIRAL",
    # Generate fill values from neighbour values by comparing "RATIO" or "DIFFERENCE" between them?
    "FILL_GENERATION_METHOD",
    # If the fill generation method is "RATIO" then what should be the maximum allowable ratio, to allow for
    # near-zero divisors?
    "MAX_ALLOWABLE_RATIO",
    # If true then the highest and lowest single partial fill values from neighbours will be dropped as a
    # further level of protection against outliers
    "TRIM_FULL_OUTLIERS"
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

A2SearchConfig = namedtuple("A2SearchConfig", [
    # An instance of SpiralSearchConfig
    # n neighbours should normally just be 8 for A2 as the search uses previously-generated
    # values on each step, i.e. "smears" values across in the direction of the search
    "SPIRAL",
    # "MEAN" or "MEDIAN" to choose how to create a single fill value from the 8 directional passes at each cell
    "PASS_AVERAGE_TYPE",
    # Generate fill values from neighbour values by comparing "RATIO" or "DIFFERENCE" between them?
    "FILL_GENERATION_METHOD",
    # If the fill generation method is "RATIO" then what should be the maximum allowable ratio, to allow for
    # near-zero divisors?
    "MAX_ALLOWABLE_RATIO"
])

GapfillFilePaths = namedtuple("GapfillFilePaths", [
    "DATA_FILES_GLOB_PATTERN",
    "SYNOPTIC_MEAN_FILE",
    "SYNOPTIC_SD_FILE",
    "COASTLINE_FILE",
    "OUTPUT_FOLDER",
    "TEMP_FOLDER"
])

