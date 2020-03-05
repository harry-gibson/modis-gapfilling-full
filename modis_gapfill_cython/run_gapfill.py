from gapfiller import Gapfiller

from gapfill_defaults import DefaultFlagValues, DefaultDespeckleConfig, \
    DefaultA1Config, DefaultA2Config, \
    DefaultDataConfig_LST, DefaultDataConfig_EVI, DefaultDataConfig_TCB, DefaultDataConfig_TCW

from modis_gapfill_cython.gapfill_config_types import GapfillFilePaths, DataSpecificConfig, GapfillJobConfig
import argparse
import os, sys
import ruamel.yaml


def main():
    gapfiller = Gapfiller(gapfillFilePaths="",
                          despeckleConfig=None, a1Config=None, a2Config=None,
                          dataSpecificConfig=None, flagValues=None,
                          fillForLatLims=None, fillForLonLims=None
                          )
    gapfiller.RunFill(onlyForDays=None, startYear=None)


def is_valid_file(parser, arg):
    if not (arg == 'SAMPLE' or os.path.exists(arg)):
        parser.error("The file %s does not exist!" % arg)
        return False
    else:
        return arg


def is_valid_directory(parser, arg):
    if os.path.isdir(arg):
        return arg
    else:
        parser.error("The directory %s does not exist!" % arg)
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This accepts a glob pattern to find data files in the MAP "
                                                 "mastergrid naming format, and paths to temporary and output "
                                                 "folders, and runs gapfilling on them using default MAP "
                                                 "parameters")
    parser.add_argument("-c", "--config", dest="runConfig", required=True,
                        help="path to the gapfill_run_config_sample_master.yml file which points to "
                             "the files that need filling and specifies which area and times "
                             "they should be filled for, and what limit values the fill should be "
                             "constrained to. Run with 'SAMPLE' for this "
                             "argument to write a sample file to the current folder which "
                             "you can edit.",
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-p" "--params", dest="fillParams", required=False,
                        help="path to the gapfill_params_defaults.yml file which configures the search "
                             "parameters for gapfilling (radius of search, number of neighbours "
                             "required, etc) for despeckle, A1, A2. If this parameter is not provided then MAP default "
                             "values will be used. Run with 'SAMPLE' for this argument to write "
                             "a sample file to the current folder which you can edit.",
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    args = parser.parse_args()

    runConfigYml = args.runConfig
    fillParamsYml = args.fillParams
    canRun = True

    yamlParser = ruamel.yaml.YAML()
    if runConfigYml == 'SAMPLE':
        defaultYamlConfigFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'gapfill_run_config_sample_master.yml')
        if not is_valid_file(defaultYamlConfigFile):
            raise RuntimeError("Sample file not found!")
        with open(defaultYamlConfigFile) as stream:
            sampleConfig = yamlParser.load(stream)
        outputYamlSampleConfigFile = os.path.join(os.getcwd(), 'gapfill_run_config_sample.yml')
        yamlParser.dump(sampleConfig, outputYamlSampleConfigFile)
        print("sample run config file written to {}. Please edit, and rerun using this as the -c parameter.")
        canRun = False
    else:
        with open(runConfigYml) as stream:
            runConfig = yamlParser.load(stream)

    # load gapfill_params_defaults always (in case not all parts provided in input), and write as sample if requested
    defaultYamlParamsFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'gapfill_params_defaults.yml')
    if not is_valid_file(defaultYamlParamsFile):
        raise RuntimeError("Default parameters file not found!")
    with open(defaultYamlParamsFile) as stream:
        defaultParams = yamlParser.load(stream)
    if fillParamsYml is None:
        runParams = defaultParams
    elif fillParamsYml == 'SAMPLE':
        outputYamlParamsFile = os.path.join(os.getcwd(), 'gapfill_run_params.yml')
        yamlParser.dump(sampleConfig, outputYamlSampleConfigFile)
        print("sample parameters file written to {}. Please edit, and rerun using this as the -p parameter.")
        canRun = False
    else:
        with open(fillParamsYml) as stream:
            runParams = yamlParser.load(stream)
        # todo: grab missing bits from default

    if not canRun:
        sys.exit(1)

def parseFilePathsConfig(runConfig):
    """ essentially maps the orderedDict from yaml parsing into the namedtuple we want to use.

    Yes, we could just use orderedDict throughout but i prefer the security of the class-like namedtuple.
    Returns GapfillFilePaths namedtuple object."""
    fP = runConfig['FilePaths']
    filePathsParsed = GapfillFilePaths(
        DATA_FILES_GLOB_PATTERN=fP['UnfilledFilesGlobPattern'],
        SYNOPTIC_MEAN_FILE=fP['UnfilledSynopticMean'],
        SYNOPTIC_SD_FILE=fP['UnfilledSynopticSD'],
        COASTLINE_FILE=fP['CoastlineTemplate'],
        OUTPUT_FOLDER=fP['OutputFolder'],
        TEMP_FOLDER=fP['TemporaryFolder'])
    return filePathsParsed

def parseDataLimitsConfig(runConfig):
    """ essentially maps the orderedDict from yaml parsing into the namedtuple we want to use.

    Yes, we could just use orderedDict throughout but i prefer the security of the class-like namedtuple.
    Returns DataSpecificConfig namedtuple object."""
    dP = runConfig['DataLimitParams']
    dataLimitsParsed = DataSpecificConfig(
        CEILING_VALUE=dP['ceiling_value'],
        FLOOR_VALUE=dP['floor_value'],
        FLOOR_CEILING_ZSCORE=dP['floor_ceiling_zscore'],
        CORRECTION_OFFSET=dP['correction_offset'],
        NODATA_VALUE=None, # not implementing this at present: read from files
        ABSOLUTE_ZERO_OFFSET=dP['absolute_zero_for_ratio']
    )
    return dataLimitsParsed


def parseJobConfig(runConfig):
    """ essentially maps the orderedDict from yaml parsing into the namedtuple we want to use.

    Yes, we could just use orderedDict throughout but i prefer the security of the class-like namedtuple.
    Returns GapfillJobConfig namedtuple object."""
    jC = runConfig['FillJob']
    jobConfigParsed = GapfillJobConfig(
        XMin=jC['XMin'],
        XMax=jC['XMax'],
        YMin=jC['YMin'],
        YMax=jC['YMax'],
        CalendarDaysToFill=jC['CalendarDaysToFill'],
        StartYear=jC['StartYear'],
        ClipMinMax=jC['ClipMinMax'],
        RunA2=jC['RunA2']
    )
    return jobConfigParsed

def