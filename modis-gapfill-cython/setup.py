# Set this to False if you don't have cython installed and want to use the pregenerated .c files
USE_CYTHON = True


from distutils.core import setup
from distutils.extension import Extension
import sys
import _version

if USE_CYTHON:
    try:
        from Cython.Distutils import build_ext
    except ImportError:
        if USE_CYTHON == 'auto':
            USE_CYTHON = False
        else:
            raise

gapfill_cython_exts = []
cmdclass = {}

if sys.version_info[0] == 2:
    base_dir = 'python2'
elif sys.version_info[0] == 3:
    # Still build from python2 code, but use build_py_2to3 to translate.
    base_dir = 'python2'
    from distutils.command.build_py import build_py_2to3 as build_py
    cmdclass.update({ 'build_py': build_py })

if USE_CYTHON:
    gapfill_cython_exts += [
        Extension(
            "gapfill_cython_core",
            ["gapfill_core_a1.pyx",
             "gapfill_core_a2.pyx",
             "gapfill_core_clamp.pyx",
             "gapfill_core_despeckle_and_flag.pyx"],
            extra_compile_args=['/openmp'],
            extra_link_args=['/openmp']
        )
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    gapfill_cython_exts += [
        Extension( "gapfill_cython_core",
            ["gapfill_core_a1.c",
             "gapfill_core_a2.c",
             "gapfill_core_clamp.c",
             "gapfill_core_despeckle_and_flag.c"],
            extra_compile_args=['/openmp'],
            extra_link_args=['/openmp'])
    ]

setup(
    name="MODIS Gapfilling algorithms",
    version = _version.__version__,
    description='Cython version of the Weiss MODIS gapfilling models',
    author = 'Harry Gibson',
    author_email='harry.s.gibson@gmail.com',
    cmdclass = cmdclass,

    ext_modules = gapfill_cython_exts
)