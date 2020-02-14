import gapfill_config

# Temporary file to run the code prior to developing a read-at-runtime config system

# The flags are set / checked using bitwise operators so should be powers of two
DefaultFlagValues = gapfill_config.FlagItems(
    OCEAN=1,
    FAILURE=2,
    EXTREME=4,
    SPECKLE=8,
    A1_FILLED=16,
    A1_FULL=32,
    A2_FILLED=64,
    CLIPPED=128)

DefaultDataConfig_LST = gapfill_config.DataCharacteristicsConfig(
    CEILING_VALUE = 100.0,
    FLOOR_VALUE = -100.0,
    FLOOR_CEILING_ZSCORE = 2.58,
    CORRECTION_OFFSET = 0,
    NODATA_VALUE = -9999,
    ABSOLUTE_ZERO_OFFSET = -273.15)

DefaultDataConfig_EVI = gapfill_config.DataCharacteristicsConfig(
    CEILING_VALUE = 1.0,
    FLOOR_VALUE = 0.0,
    FLOOR_CEILING_ZSCORE = 2.58,
    CORRECTION_OFFSET = 0,
    NODATA_VALUE = -9999,
    ABSOLUTE_ZERO_OFFSET = 0)

DefaultDataConfig_TCW = gapfill_config.DataCharacteristicsConfig(
    CEILING_VALUE = 2.0,
    FLOOR_VALUE = -1.0,
    FLOOR_CEILING_ZSCORE = 2.58,
    CORRECTION_OFFSET = 0,
    NODATA_VALUE = -9999,
    ABSOLUTE_ZERO_OFFSET = 0)

DefaultDataConfig_TCB = gapfill_config.DataCharacteristicsConfig(
    CEILING_VALUE = 2.0,
    FLOOR_VALUE = -1.0,
    FLOOR_CEILING_ZSCORE = 2.58,
    CORRECTION_OFFSET = 0,
    NODATA_VALUE = -9999,
    ABSOLUTE_ZERO_OFFSET = 0)

# The despeckle and gapfill algorithms search for neighbour cells in an outward spiral i.e. an
# increasing search radius until enough neighbours are found for the checks or until the maximum radius is reached.
# Define the size of the search and the number of neighbours required. For A1 this determines the fill radius i.e.
# the size of gap that can be filled.
# All MODIS gapfilling done in MAP since 2014 has used the following parameters which give an effective search
# radius of approx 31 pixels (i.e. sqrt(3142/pi)) and
DefaultDespeckleSpiralConfig = gapfill_config.SpiralSearchConfig(
    MAX_NBRS_TO_SEARCH = 3142,
    MIN_REQUIRED_NBRS = 320,
    MAX_USED_NBRS = 640)

DefaultA1SpiralConfig = gapfill_config.SpiralSearchConfig(
    MAX_NBRS_TO_SEARCH = 3142,
    MIN_REQUIRED_NBRS = 480,
    MAX_USED_NBRS = 960)

# for A2 we only normally search the immediate neighbours (as the subsequent steps use previous step outputs,
# so values are smeared across anyway), but this does not need to be the case
DefaultA2SpiralConfig = gapfill_config.SpiralSearchConfig(
    MAX_NBRS_TO_SEARCH = 8,
    MIN_REQUIRED_NBRS = 1,
    MAX_USED_NBRS = 8)

DefaultDespeckleConfig = gapfill_config.DespeckleConfig(
    SPIRAL = DefaultDespeckleSpiralConfig,
    EXTREME_BEYOND_SD = 2.58,
    SPECKLE_BEYOND_SD = 1.64,
    SPECKLE_NBR_Z_THRESH = 0.2)

DefaultA1Config = gapfill_config.A1SearchConfig(
    SPIRAL = DefaultA1SpiralConfig,
    FILL_GENERATION_METHOD = "DIFFERENCE",
    MAX_ALLOWABLE_RATIO = 0,
    TRIM_FULL_OUTLIERS = True)

DefaultA2Config = gapfill_config.A2SearchConfig(
    SPIRAL = DefaultA2SpiralConfig,
    PASS_AVERAGE_TYPE = "MEAN",
    FILL_GENERATION_METHOD = "DIFFERENCE",
    MAX_ALLOWABLE_RATIO = 0)
