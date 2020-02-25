from gapfiller import Gapfiller

from gapfill_defaults import DefaultFlagValues, DefaultDespeckleConfig, \
    DefaultA1Config, DefaultA2Config, \
    DefaultDataConfig_LST, DefaultDataConfig_EVI, DefaultDataConfig_TCB, DefaultDataConfig_TCW

import argparse

def main()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This accepts a glob pattern to find data files in the MAP "
                                                 "mastergrid naming format, and paths to temporary and output "
                                                 "folders, and runs gapfilling on them using default MAP "
                                                 "parameters")
    parser.add_argument("-f", "--filepattern", dest="globPattern", required=True,
                        help="glob pattern to match the mastergrid-formatted filenames for gapfilling")
    parser.add_argument("-m", "--meanfile", dest="meanFile", required=True,
                        help="path to the synoptic mean file for use in gapfilling")
    parser.add_argument("-s", "--sdfile", dest="sdFile", required=True,
                        help="path to the synoptic SD file for use in gapfilling")
    parser.add_argument("-c", "--coastlinefile", dest="coastlineFile", required=True,
                        help="path to the coastline template file to identify land/sea areas")
    parser.add_argument("-a1", "--a1config", dest="a1Config", required=False,
                        help="A1 gapfill configuration")
    parser.add_argument("-a2", "--a2config", dest="a2Config", required=False,
                        help="A2 gapfill configuration")
