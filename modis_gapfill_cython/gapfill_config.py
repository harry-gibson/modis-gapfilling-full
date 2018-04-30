# The flags output is an 8 bit raster which represents 8 separate flag conditions as a bitmask, defined here
FlagValues = {
    "OCEAN" : 1,
    "FAILURE" : 2,
    "EXTREME": 4,
    "SPECKLE" : 8,
    "A1_FILLED" : 16,
    "A1_FULL" : 32,
    "A2_FILLED" : 64,
    "CLIPPED" : 128
}

DataSpecificConfig = {
    # Hard upper / lower limits that the fill values will be clipped to.
    # If one or the other is set to the same as the NODATA_VALUE then values won't be clipped in that
    # direction. If both are set to this then no clipping will be done.
    "DATA_UPPER_LIMIT": 2.0,
    "DATA_LOWER_LIMIT": -1.0,
    # An offset may be specified which will be applied (added) to all data before any further processing
    # for example to convert from kelvin to celsius specify 273.15
    "DATA_CORRECTION_OFFSET": 0,
    # A value may be specified which will be ignored; this should generally be the nodata value in the
    # GeoTIFF files
    "NODATA_VALUE": -9999,
    # If the fill generation method is ratio, this makes no sense if the data are not ratiometric/absolute
    # for example temperature in celsius. It may be possible to convert the data into an absolute scale to make
    # a ratio-based fill appropriate, for example by putting a celsius temperature into kelvin
    "DATA_ABSOLUTE_ZERO_OFFSET": 0
}

# The de-speckle search runs in an outward spiral i.e. an increasing search radius until enough neighbours are
# found for the checks or until the maximum radius is reached. Define the size of the search and the number of
# neighbours required here:
# todo define a namedtuple type for this
DespeckleSearchConfig = {
    # Number of cells to search; the radius in pixel distance terms is approx sqrt(value/pi)
    "MAX_NBRS_TO_SEARCH": 3142,
    # Search is successful if we find at least this many within the max number
    "MIN_REQ_NBRS": 320,
    # Stop searching early after this many even if we haven't gone to the full radius
    "MAX_USED_NBRS": 640
}

DespeckleThresholdConfig = {
    # Number of stds beyond the mean beyond which a value will be unconditionally discarded in despeckle
    "EXTREME_BEYOND_SD": 2.58,
    # Number of stds beyond the mean beyond which a value MAY be discarded as being a speckle, depending
    # on neighborhood similarity
    "SPECKLE_BEYOND_SD": 1.64,
    # Max difference in Z-score between a potential speckle and the average of its neighbours, for it
    # to be accepted as not being a speckle
    "SPECKLE_NBR_Z_THRESH": 0.2 #
}

# The A1 search runs in an outward spiral i.e. an increasing search radius until enough neighbours are
# found for the checks or until the maximum radius is reached. Define the size of the search and the number of
# neighbours required here:
# todo define a namedtuple type for this
A1SearchConfig = {
    # Number of cells to search; the radius in pixel distance terms is approx sqrt(value/pi)
    "MAX_NBRS_TO_SEARCH": 3142,
    # Search is successful if we find at least this many within the max number
    "MIN_REQ_NBRS": 480,
    # Stop searching early after this many even if we haven't gone to the full radius
    "MAX_USED_NBRS": 960,
    # Generate fill values from neighbour values by comparing "RATIO" or "DIFFERENCE" between them?
    "FILL_GENERATION_METHOD": "DIFFERENCE",
    # If the fill generation method is "RATIO" then what should be the maximum allowable ratio, to allow for
    # near-zero divisors?
    "MAX_ALLOWABLE_RATIO": 5.0,
    # If true then the highest and lowest single partial fill values from neighbours will be dropped as a
    # further level of protection against outliers
    "TRIM_FILL_OUTLIERS": True
}

A2SearchConfig = {
    # Number of cells to search; the radius in pixel distance terms is approx sqrt(value/pi)
    # This should normally just be immediate neighbours (8) for A2 as the search uses previously-generated
    # values on each step, i.e. "smears" values across in the direction of the search
    "MAX_NBRS_TO_SEARCH": 8,
    # "MEAN" or "MEDIAN" to choose how to create a single fill value from the 8 directional passes at each cell
    "PASS_SELECTOR_MEAN_OR_MEDIAN": "MEAN",
    # Generate fill values from neighbour values by comparing "RATIO" or "DIFFERENCE" between them?
    "FILL_GENERATION_METHOD": "DIFFERENCE",
    # If the fill generation method is "RATIO" then what should be the maximum allowable ratio, to allow for
    # near-zero divisors?
    "MAX_ALLOWABLE_RATIO": 5.0
}