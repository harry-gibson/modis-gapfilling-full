from typing import NamedTuple


class FlagItems(NamedTuple):
    """The flags output is an 8 bit raster holding 8 separate flag conditions as a bitmask, values defined here"""
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
        """Factory method to create instance of this NamedTuple from YAML config, which must contain "Flags" key"""
        fP = runParams['Flags']
        return cls(int(fP['ocean']), int(fP['failure']), int(fP['extreme']), int(fP['speckle']),
                   int(fP['a1_filled']), int(fP['a1_full']), int(fP['a2_filled']), int(fP['clipped']))


class GapfillJobConfig(NamedTuple):
    """The gapfill job config defines the bounding box, the calendar days, and years to fill, as well as whether or
     not the A2 algorithm should be run in addition to A1, and whether or not to clamp the output to the floor/ceiling
      values"""
    XMin_Deg: float = None
    YMin_Deg: float = None
    XMax_Deg: float = None
    YMax_Deg: float = None
    CalendarDaysToFill: list = None
    StartYear: int = None
    ClipMinMax: bool = True
    RunA2: bool = True
    NCores: int = None
    MemTargetBytes: int = None

    @classmethod
    def from_yaml_config(cls, runConfig):
        """Factory method to create instance of this NamedTuple from YAML config, which must contain "FillJob" key"""
        jC = runConfig['FillJob']
        jobConfigParsed = cls(
            XMin_Deg=float(jC['XMin']),
            XMax_Deg=float(jC['XMax']),
            YMin_Deg=float(jC['YMin']),
            YMax_Deg=float(jC['YMax']),
            CalendarDaysToFill=list(jC['CalendarDaysToFill']),
            StartYear=float(jC['StartYear']),
            ClipMinMax=float(jC['ClipMinMax']),
            RunA2=float(jC['RunA2']),
            NCores=int(jC['NCores']),
            MemTargetBytes=int(jC['MemTargetBytes'])
        )
        return jobConfigParsed


class SpiralSearchConfig(NamedTuple):
    """ The spiral search config determines the number of neighbour cells to search and thus the gapfill radius
    Attributes:
        MAX_NBRS_TO_SEARCH: Number of cells to search; the radius in pixel distance terms is approx sqrt(value/pi)
        MIN_REQUIRED_NBRS   Search is successful if we find at least this many within the max number
        MAX_USED_NBRS       Stop searching early after this many even if we haven't gone to the full radius
    """
    MAX_NBRS_TO_SEARCH: int
    MIN_REQUIRED_NBRS: int
    MAX_USED_NBRS: int

    @classmethod
    def from_yaml_config(cls, singleFillConfig):
        """Factory method to create instance of this NamedTuple from YAML config, which must contain "Spiral" key"""
        sC = singleFillConfig['spiral']
        spiralParsed = cls(
            MAX_NBRS_TO_SEARCH=int(sC['max_nbrs_to_search']),
            MIN_REQUIRED_NBRS=int(sC['min_required_nbrs']),
            MAX_USED_NBRS=int(sC['max_used_nbrs'])
        )
        return spiralParsed


class DespeckleConfig(NamedTuple):
    """ The Despeckle config configures how far the despeckle algorithm will search and what counts as as speckle.
    Attributes:
        SPIRAL_CONFIG:      Instance of SpiralSearchConfig, gives radius of search and required density of neighbours
        EXTREME_BEYOND_SD:  Number of SDs from the mean beyond which a cell is definitely speckle
        SPECKLE_BEYOND_SD:  Number of SDs from the mean beyond which a cell may be speckle, depending on neighbours
        SPECKLE_NBR_Z_THRESH: Number of SDs a 'maybe speckle' can differ from neighbours to count as 'similar'
    """
    SPIRAL_CONFIG: SpiralSearchConfig
    EXTREME_BEYOND_SD: float
    SPECKLE_BEYOND_SD: float
    SPECKLE_NBR_Z_THRESH: float

    @classmethod
    def from_yaml_config(cls, runParams):
        """Factory method to create instance of this NamedTuple from YAML config, which must contain "Despeckle" key"""
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
    """ The A1 config configures how far the A1 algorithm will search and  how fill values are calculated.
        Attributes:
            SPIRAL_CONFIG:      Instance of SpiralSearchConfig, gives radius of search and required density of neighbours
            FILL_GENERATION_METHOD:  "RATIO" or "DIFFERENCE": how a cell's value is compared to a neighbour
            MAX_ALLOWABLE_RATIO:  For "RATIO", max allowed ratio, to guard against infinity with near-zero divisors
            TRIM_FULL_OUTLIERS: Drop most extreme partial fill values from neighbours, to guard against outliers?
    """
    SPIRAL_CONFIG: SpiralSearchConfig
    FILL_GENERATION_METHOD: str
    MAX_ALLOWABLE_RATIO: float
    TRIM_FULL_OUTLIERS: bool

    @classmethod
    def from_yaml_config(cls, runParams):
        """Factory method to create instance of this NamedTuple from YAML config, which must contain "A1" key"""
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
    """ The A2 config configures how far the A1 algorithm will search and  how fill values are calculated.
            Attributes:
                SPIRAL_CONFIG:      Instance of SpiralSearchConfig, gives radius of search and required density of
                                    neighbours. Expected to be configured to search immediate neighbours only.
                FILL_GENERATION_METHOD:  "RATIO" or "DIFFERENCE": how a cell's value is compared to a neighbour
                MAX_ALLOWABLE_RATIO:  For "RATIO", max allowed ratio, to guard against infinity with near-zero divisors
                PASS_AVERAGE_TYPE:  "MEAN" / "MEDIAN": how to select fill value from the 8 directional passes
    """

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
        """Factory method to create instance of this NamedTuple from YAML config, which must contain "A2" key"""
        a2C = runParams['A2']
        spiralParsed = SpiralSearchConfig.from_yaml_config(a2C)
        configParsed = cls(
            SPIRAL_CONFIG=spiralParsed,
            FILL_GENERATION_METHOD="RATIO" if bool(a2C['use_ratio_fill']) else "DIFFERENCE",
            MAX_ALLOWABLE_RATIO=float(a2C['max_allowable_ratio']),
            PASS_AVERAGE_TYPE=a2C['pass_averaging_method']
        )
        return configParsed


class DataLimitsConfig(NamedTuple):
    """ Configures numerical parameters specific to this dataset.
    Attributes:
        CEILING_VALUE:          Hard upper limit to clip fill values to. Use NODATA_VALUE to not clip in this direction
        FLOOR_VALUE:            Hard lower limit to clip fill values to. Use NODATA_VALUE to not clip in this direction
        FLOOR_CEILING_ZSCORE:   Clip fill values more than this number of SD from the mean, iif ceiling/floor also set
        CORRECTION_OFFSET:      Value to add to all data before further processing if required
        ABSOLUTE_ZERO_OFFSET:   Offset to convert data into an absolute scale if appropriate, e.g. -273.15 to convert
                                celsius into kelvin so that ratio-based comparisons can validly be made
        NODATA_VALUE:           Value to treat as nodata. Will not currently be read from files directly.
                                Must currently be the same across all data files.
    """
    CEILING_VALUE: float
    FLOOR_VALUE: float
    FLOOR_CEILING_ZSCORE: float
    CORRECTION_OFFSET: float
    ABSOLUTE_ZERO_OFFSET: float
    NODATA_VALUE: float

    @classmethod
    def from_yaml_config(cls, runConfig):
        """Factory method to create instance from YAML config, which must contain "DataLimitParams" key"""
        dL = runConfig['DataLimitParams']
        configParsed = cls(
            CEILING_VALUE=float(dL['ceiling_value']),
            FLOOR_VALUE=float(dL['floor_value']),
            CORRECTION_OFFSET=float(dL['correction_offset']),
            FLOOR_CEILING_ZSCORE=float(dL['floor_ceiling_zscore']),
            ABSOLUTE_ZERO_OFFSET=float(dL['absolute_zero_for_ratio']),
            NODATA_VALUE=float(dL['nodata_value'])
        )
        return configParsed

    def GetSummaryMessage(self):
        print("Despeckle: Rejecting data beyond {0!s}s.d. of mean. Nbr search on data beyond {1!s} s.d. of mean.".
              format(stDevValidityThreshold, speckleDevThreshold))
        print("Nbr searching for {0!s} - {1!s} nbrs within {2!s} spiral steps for z-score tolerance of {3!s}".
              format(_SPECKLE_NBR_MIN_THRESHOLD, _SPECKLE_NBR_MAX_THRESHOLD, _MAX_NEIGHBOURS_TO_CHECK,
                     _SPECKLE_ZSCORE_THRESHOLD))

class GapfillFilePaths(NamedTuple):
    DATA_FILES_GLOB_PATTERN: str
    SYNOPTIC_MEAN_FILE: str
    SYNOPTIC_SD_FILE: str
    COASTLINE_FILE: str
    OUTPUT_FOLDER: str
    TEMP_FOLDER: str

    @classmethod
    def from_yaml_config(cls, runConfig):
        """Factory method to create instance of this NamedTuple from YAML config, which must contain "FilePaths" key"""
        fP = runConfig['FilePaths']
        configParsed = cls(
            DATA_FILES_GLOB_PATTERN=fP['UnfilledFilesGlobPattern'],
            SYNOPTIC_MEAN_FILE=fP['UnfilledSynopticMean'],
            SYNOPTIC_SD_FILE=fP['UnfilledSynopticSD'],
            COASTLINE_FILE=fP['CoastlineTemplate'],
            OUTPUT_FOLDER=fP['OutputFolder'],
            TEMP_FOLDER=fP['TemporaryFolder']
        )
        # todo implement path validity / file existence checks here
        return configParsed


class DespeckleDiagnostics(NamedTuple):
    SpeckleCellCount: int
    ExtremeCellCount: int
    GoodCellCount: int
    OceanCellCount: int
    ClearedSpeckleCellCount: int
    TimeSeconds: float

    def GetSummaryMessage(self):
        message = f"""
    Despeckle report:
        Speckle cell count:             {self.SpeckleCellCount}
        Extreme cell count:             {self.ExtremeCellCount}
        Good cell count:                {self.GoodCellCount}
        Ocean cell count:               {self.OceanCellCount}
        Cleared speckle count:          {self.ClearedSpeckleCellCount}
        Total time for despeckle:       {round(self.TimeSeconds, 2)}s
        """
        return message


class A1Diagnostics(NamedTuple):
    TotalCellsSearched: int
    CellsWithGoodData: int
    OceanCells: int
    NeverDataLocations: int
    GapCellsTotal: int
    GapCellsTooBig: int
    PermanentGapCells: int
    GapCellsFullyFilled: int
    GapCellsPartFilled: int
    FailedInsufficientPairs: int
    FailedNoPairs: int
    TotalAlternateYearsScanned: int
    TotalNbrsChecked: int
    TotalNbrsUsed: int
    TimeSeconds: float

    def GetSummaryMessage(self):
        message = f"""
    A1 report:
        Total cells scanned:            {self.TotalCellsSearched}
        Total cells with good data:     {self.CellsWithGoodData}
        Total cells ocean:              {self.OceanCells}
        Total cells never-data:         {self.NeverDataLocations}
        Total processed gaps:           {self.GapCellsTotal}
        Total gaps too large too fill:  {self.GapCellsTooBig}
        Total gaps at unfillable locs   {self.PermanentGapCells}
            (i.e. never-data * years)
        Total cells filled fully:       {self.GapCellsFullyFilled}
        Total cells filled partially:   {self.GapCellsPartFilled}
        Total w/insufficient nbr pairs: {self.FailedInsufficientPairs}
        Total with no nbr pairs:        {self.FailedNoPairs}
        Total yrs scanned f/b in stack: {self.TotalAlternateYearsScanned}
        Total nbr cells scanned:        {self.TotalNbrsChecked}
        Total used nbr pairs:           {self.TotalNbrsUsed}
        Total time for A1:              {round(self.TimeSeconds, 2)}s
        """
        return message


class A2Diagnostics(NamedTuple):
    GapCellsTotal: int
    GapCellsFilled: int
    GapCellsFilledByPass: list
    TimeSeconds: float
    TimeSecondsByPass: list

    def GetSummaryMessage(self):
        tByPass = " / ".join([f"{p}: {round(t, 2)}s" for p, t in enumerate(self.TimeSecondsByPass)])
        nByPass = " / ".join([f"{p}: {n}" for p, n in enumerate(self.GapCellsFilledByPass)])
        message = f"""
    A2 report:
        Total processed gaps:           {self.GapCellsTotal}
        Total cells filled:             {self.GapCellsFilled}
            Filled by pass:             {nByPass}
        Total time for A2:              {self.TimeSeconds}s
            Time by pass:               {tByPass}
        """
        return message


