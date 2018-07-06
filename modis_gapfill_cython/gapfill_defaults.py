from gapfill_config import *

# Temporary file to run the code prior to developing a read-at-runtime config system

# The flags are set / checked using bitwise operators so should be powers of two
DefaultFlagValues = FlagItems(
    OCEAN=1,
    FAILURE=2,
    EXTREME=4,
    SPECKLE=8,
    A1_FILLED=16,
    A1_FULL=32,
    A2_FILLED=64,
    CLIPPED=128)



DataSpecificDefaultConfig = DataCharacteristicsConfig(
    CEILING_VALUE=2.0,
    FLOOR_VALUE=-1.0,
    FLOOR_CEILING_ZSCORE=2.58,
    CORRECTION_OFFSET=0,
    NODATA_VALUE=-9999,
    ABSOLUTE_ZERO_OFFSET=0)

# The despackle and gapfill algorithms search for neighbour cells in an outward spiral i.e. an
# increasing search radius until enough neighbours are found for the checks or until the maximum radius is reached.
# Define the size of the search and the number of neighbours required. For A1 this determines the fill radius i.e.
# the size of gap that can be filled.
DespeckleSpiralDefaultConfig = SpiralSearchConfig(MAX_NBRS_TO_SEARCH=3142, MIN_REQUIRED_NBRS=320, MAX_USED_NBRS=640)
A1SpiralDefaultConfig = SpiralSearchConfig(MAX_NBRS_TO_SEARCH=3142, MIN_REQUIRED_NBRS=480, MAX_USED_NBRS=960)
A2SpiralDefaultConfig = SpiralSearchConfig(MAX_NBRS_TO_SEARCH=8, MIN_REQUIRED_NBRS=1, MAX_USED_NBRS=8)

DespeckleThresholdDefaultConfig = SpeckleSelectionConfig(SPIRAL=DespeckleSpiralDefaultConfig,
                                                         EXTREME_BEYOND_SD=2.58, SPECKLE_BEYOND_SD=1.64,
                                                         SPECKLE_NBR_Z_THRESH=0.2)

A1DefaultParams = A1SearchConfig(SPIRAL=A1SpiralDefaultConfig,
                                 FILL_GENERATION_METHOD="DIFFERENCE", MAX_ALLOWABLE_RATIO=5.0,
                                 TRIM_FULL_OUTLIERS=True)

A2DefaultParams = A2SearchConfig(SPIRAL=A2SpiralDefaultConfig,
                                 PASS_AVERAGE_TYPE="MEAN", FILL_GENERATION_METHOD="DIFFERENCE",
                                 MAX_ALLOWABLE_RATIO=5.0)
