import numpy as np
cimport cython
from libc.math cimport sqrt
#cdef int _SEARCH_RADIUS = 10
from cython.parallel import prange, parallel


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef MinMaxClip(float[:,::1] dataImage,
                 unsigned char[:,::1] flagsImage,
                 float[:,::1] meansImage,
                 float[:,::1] stdImage,
                 unsigned char flagToCheck,
                 unsigned char flagToSet,
                 float floor_ceiling_value,
                 float _NDV,
                 float upperHardLimit,
                 float lowerHardLimit
                 ):
    '''
    Clips (clamps) an image to not exceed +- n std from the mean image, or a hard upper / lower limit

    '''
    cdef:
        Py_ssize_t x, y, xShape, yShape
        float value, minAllowed, maxAllowed
        #char localFlag
        float negInf, posInf

    yShape = dataImage.shape[0]
    xShape = dataImage.shape[1]
    assert xShape == meansImage.shape[1]
    assert xShape == stdImage.shape[1]
    assert xShape == flagsImage.shape[1]
    assert yShape == meansImage.shape[0]
    assert yShape == stdImage.shape[0]
    assert yShape == flagsImage.shape[0]

    posInf = np.inf
    negInf = -np.inf

    #with nogil, parallel(num_threads=20):
    if 1:
        for y in range(yShape):
            value = -1
            maxAllowed = posInf
            minAllowed = negInf
            for x in range(xShape):
                if (flagsImage[y, x] & flagToCheck) != flagToCheck:
                    continue
                value = dataImage[y,x]
                if value == _NDV:
                    continue
                maxAllowed = meansImage[y, x] + (floor_ceiling_value * stdImage[y, x])
                minAllowed = meansImage[y, x] - (floor_ceiling_value * stdImage[y, x])
                if maxAllowed>=200.0:
                    print ("Whoops! Location {0!s},{1!s} had value {2!s}. Mean={3!s} and std={4!s} giving maxallowed of {5!s}"
                           .format(x,y,value,meansImage[y,x],stdImage[y,x],maxAllowed)
                    )
                    # crash
                    assert False
                if minAllowed<=-200.0:
                    print ("Whoops! Location {0!s},{1!s} had value {2!s}. Mean={3!s} and std={4!s} giving minallowed of {5!s}"
                           .format(x,y,value,meansImage[y,x],stdImage[y,x],minAllowed)
                    )
                    # crash
                    assert False
                if maxAllowed > upperHardLimit:
                    maxAllowed = upperHardLimit
                if minAllowed < lowerHardLimit:
                    minAllowed = lowerHardLimit

                if value > maxAllowed:
                    dataImage[y, x] = maxAllowed
                    flagsImage[y, x] = flagsImage[y, x] | flagToSet
                    continue
                if value < minAllowed:
                    dataImage[y, x] = minAllowed
                    flagsImage[y, x] = flagsImage[y, x] | flagToSet
                    continue

cpdef MinMaxClip3D(float[:,:,::1] dataImage,
                 unsigned char[:,:,::1] flagsImage,
                 float[:,::1] meansImage,
                 float[:,::1] stdImage,
                 unsigned char flagToCheck,
                 unsigned char flagToSet,
                 float floor_ceiling_value,
                 float _NDV,
                 float upperHardLimit,
                 float lowerHardLimit
                 ):
    '''
    Clips / clamps a stack of images to not exceed +- n stds from a mean image or a hard upper/lower limit
    '''
    cdef:
        Py_ssize_t zSize, z
    zSize = dataImage.shape[0]
    assert zSize == flagsImage.shape[0]
    for z in range(zSize):
        MinMaxClip(dataImage[z], flagsImage[z],
                   meansImage, stdImage,
                   flagToCheck, flagToSet, floor_ceiling_value, _NDV, upperHardLimit, lowerHardLimit)