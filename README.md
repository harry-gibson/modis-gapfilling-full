This repository contains gapfilling code used in MAP to fill gaps in MODIS LST, EVI, TCB and TCW imagery.

The repository presently contains two versions of the code:

Firstly, the "original" implementation contained within an IPython notebook `Gapfill_BySlice_Cython_v04.ipynb`. 
This notebook contains all necessary code to run the gapfilling but is hard to do in a systematic way due to the need to run cells manually and babysit it.

Secondly the main implementation, which is a commandline script contained in `run_gapfill.py`. This requires one or more yaml files for configuration and can be run unattended with the appropriate parameters set in these.

In either case a python environment is required which contains Cython and a working compiler for the platform - on Windows this means installing the Windows 10 SDK, through Visual Studio (other less reliable and constantly breaking workarounds are available).
 
At present there is a dependency on another MAP software package which is also available in our github / gitlab repositories, but it is intended to remove this dependency asap so further details are not provided here (correct as of 2020-03-24).

To set up the commandline version: clone this repository, change to the `modis_gapfill_cython` subdirectory, and run `python setup.py build_ext --inplace` (using an appropriate python environment). Assuming the build proceeds without errors, in the parent folder (repository root) run `python run_gapfill.py` with appropriate parameters to run the gapfill (use `-h` to get help on generating a configuration file)
  
Reference for the original gapfilling algorithms re-implemented here: 
Weiss, D.J., Atkinson, P.M., Bhatt, S., Mappin, B., Hay, S.I. & Gething, P.W. (2014) An effective approach for gap-filling continental scale remotely sensed time-series. ISPRS Journal of Photogrammetry and Remote Sensing, 98, 106-118
