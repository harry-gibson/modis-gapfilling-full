# Set this to False if you don't have cython installed and want to use the pregenerated .c files
USE_CYTHON = True


from distutils.core import setup
from distutils.extension import Extension
import sys
import os
import platform
#import _version

if USE_CYTHON:
    try:
        from Cython.Distutils import build_ext
    except ImportError:
        if USE_CYTHON == 'auto':
            USE_CYTHON = False
        else:
            raise

def scandir(dir, files=[]):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if USE_CYTHON:
            ext = ".pyx"
        else:
            ext = ".c"
        if os.path.isfile(path) and path.endswith(ext):
            if dir == ".":
                files.append(file[:-4])
            else:
                files.append((path.replace(os.path.sep, ".")[:-4]))
        elif os.path.isdir(path):
            # recursive scan
            scandir(path, files)
    return files

def makeExtension(extName):
    if USE_CYTHON:
        ext = ".pyx"
    else:
        ext = ".c"
    extPath = extName.replace(".", os.path.sep) + ext
    # there's probably a more robust way to do this check; we actually care whether it's GCC or intel compiler
    # rather than the OS. The compilers need a different switch syntax for triggering openmp
    if platform.uname()[0] == 'Windows':
        compArgs = ['/openmp', '-O3']
        linkArgs = ['/openmp']
    else:
        compArgs = ['-fopenmp', '-O3']
        linkArgs = ['-fopenmp']
    return Extension(
        extName,
        [extPath],
        include_dirs=["."],
        extra_compile_args=compArgs,
        extra_link_args=linkArgs
    )

extensionNames = scandir(".")
extensions = [makeExtension(name) for name in extensionNames]

#x
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
    cmdclass.update({'build_ext': build_ext})

setup(
    name="MODIS Gapfilling algorithms",
    #version = _version.__version__,
    description='Cython version of the Weiss MODIS gapfilling models',
    author = 'Harry Gibson',
    author_email='harry.s.gibson@gmail.com',
    cmdclass = cmdclass,
    ext_modules = extensions
)