from modis_gapfill_cython.gapfiller import GapFiller

from modis_gapfill_cython.gapfill_config_types import GapfillFilePaths, GapfillJobConfig, DataLimitsConfig, \
    DespeckleConfig, A1SearchConfig, A2SearchConfig, FlagItems

import argparse
import os, sys
import ruamel.yaml


def main(gapfillFilePaths: GapfillFilePaths,
         despeckleConfig: DespeckleConfig,
         a1Config: A1SearchConfig,
         a2Config: A2SearchConfig,
         dataSpecificConfig: DataLimitsConfig,
         flagValues: FlagItems,
         jobDetails: GapfillJobConfig):
    gapfiller = GapFiller(gapfillFilePaths=gapfillFilePaths,
                          despeckleConfig=despeckleConfig, a1Config=a1Config, a2Config=a2Config,
                          dataSpecificConfig=dataSpecificConfig, flagValues=flagValues,
                          jobDetails=jobDetails
                          )
    gapfiller.RunFill()


def is_valid_file(arg):
    if not (arg == 'SAMPLE' or os.path.exists(arg)):
        parser.error("The file %s does not exist!" % arg)
        return False
    else:
        return arg


def is_valid_directory(arg):
    if os.path.isdir(arg):
        return arg
    else:
        parser.error("The directory %s does not exist!" % arg)
        return False


def try_load_default_limits(filePaths: GapfillFilePaths):
    aFileName = filePaths.DATA_FILES_GLOB_PATTERN
    basename = os.path.basename(aFileName)
    defaultYamlLimitsFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                         'config_defaults',
                                         'MAP_limits_defaults.yml')
    yamlParser = ruamel.yaml.YAML()
    with open(defaultYamlLimitsFile) as stream:
        allDefaultLimits = yamlParser.load(stream)
    matched = [o for o in allDefaultLimits.keys() if o.upper() in basename.upper()]
    if len(matched) == 1:
        print("Found default limit configuration for variable {} based on filename pattern {}".format(
            matched[0], aFileName))
        return {'DataLimitParams': allDefaultLimits[matched[0]]}
    else:
        raise KeyError("No default limit configuration found based on filename pattern {}".format(aFileName))


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
                        type=is_valid_file)
    parser.add_argument("-p" "--params", dest="fillParams", required=False,
                        help="path to the gapfill_params_defaults.yml file which configures the search "
                             "parameters for gapfilling (radius of search, number of neighbours "
                             "required, etc) for despeckle, A1, A2. If this parameter is not provided then MAP default "
                             "values will be used. Run with 'SAMPLE' for this argument to write "
                             "a sample file to the current folder which you can edit.",
                        metavar="FILE",
                        type=is_valid_file)
    args = parser.parse_args()

    runConfigYml = args.runConfig
    fillParamsYml = args.fillParams
    canRun = True

    yamlParser = ruamel.yaml.YAML()
    if runConfigYml == 'SAMPLE':
        defaultYamlConfigFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'config_defaults',
                                             'gapfill_run_config_sample_master.yml')
        if not is_valid_file(defaultYamlConfigFile):
            raise RuntimeError("Sample file not found!")
        with open(defaultYamlConfigFile) as stream:
            sampleConfig = yamlParser.load(stream)
        outputYamlSampleConfigFile = os.path.join(os.getcwd(),
                                                  'gapfill_run_config_sample.yml')
        with open(outputYamlSampleConfigFile, 'w') as stream:
            yamlParser.dump(sampleConfig, stream)
        print("sample run config file written to {}. Please edit, and rerun using this as the -c parameter."
              .format(outputYamlSampleConfigFile))
        canRun = False
    else:
        with open(runConfigYml) as stream:
            runConfig = yamlParser.load(stream)

    # load gapfill_params_defaults always (in case not all parts provided in input), and write as sample if requested
    defaultYamlParamsFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                         'config_defaults',
                                        'gapfill_params_defaults.yml')
    if not is_valid_file(defaultYamlParamsFile):
        raise RuntimeError("Default parameters file not found!")
    with open(defaultYamlParamsFile) as stream:
        defaultParams = yamlParser.load(stream)
    if fillParamsYml is None:
        runParams = defaultParams
    elif fillParamsYml == 'SAMPLE':
        outputYamlParamsFile = os.path.join(os.getcwd(),
                                            'gapfill_run_params_sample.yml')
        with open(outputYamlParamsFile, 'w') as stream:
            yamlParser.dump(sampleConfig, stream)
        print("sample parameters file written to {}. Please edit, and rerun using this as the -p parameter."
              .format(outputYamlParamsFile))
        canRun = False
    else:
        with open(fillParamsYml) as stream:
            runParams = yamlParser.load(stream)
        # todo: grab missing bits from default

    if not canRun:
        sys.exit(1)
    try:
        filePaths = GapfillFilePaths.from_yaml_config(runConfig)
    except KeyError:
        raise ValueError("Could not parse file paths from yaml config")
    try:
        dataLimits = DataLimitsConfig.from_yaml_config(runConfig)
    except KeyError:
        yamlishBlock = try_load_default_limits(filePaths)
        dataLimits = DataLimitsConfig.from_yaml_config(yamlishBlock)
    try:
        jobConfig = GapfillJobConfig.from_yaml_config(runConfig)
    except KeyError:
        raise ValueError("Could not parse job request from yaml config")
    try:
        despeckleConfig = DespeckleConfig.from_yaml_config(runParams)
    except KeyError:
        try:
            print("Loading default despeckle config as none was specified")
            despeckleConfig = DespeckleConfig.from_yaml_config(defaultParams)
        except KeyError:
            raise ValueError("Could not parse despeckle configuration from default yaml config")
    try:
        a1Config = A1SearchConfig.from_yaml_config(runParams)
    except KeyError:
        try:
            print("Loading default A1 config as none was specified")
            a1Config = A1SearchConfig.from_yaml_config(defaultParams)
        except:
            raise ValueError("Could not parse A1 configuration from default yaml config")
    try:
        a2Config = A2SearchConfig.from_yaml_config(runParams)
    except KeyError:
        if jobConfig.RunA2:
            try:
                print("Loading default A2 config as none was specified")
                a2Config = A2SearchConfig.from_yaml_config(defaultParams)
            except:
                raise ValueError("A2 fill was requested but could not parse A2 configuration from default yaml config")
    try:
        flagItems = FlagItems.from_yaml_config(runParams)
    except KeyError:
        try:
            print("Loading default flag values")
            flagItems = FlagItems.from_yaml_config(defaultParams)
        except:
            raise ValueError("Could not parse flag values from default yaml config")

    main(gapfillFilePaths=filePaths, despeckleConfig=despeckleConfig, a1Config=a1Config,
         a2Config=a2Config, dataSpecificConfig=dataLimits, flagValues=flagItems,
         jobDetails=jobConfig)
