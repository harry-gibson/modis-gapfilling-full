cimport cython
import numpy as np

cdef class A1DataStack:
    """ Holds the various data stacks needed to run the a1_core algorithm """
    cdef public float[:,:,::1] DataArray3D
    cdef public unsigned char [:,:,::1] FlagsArray3D
    cdef public float[:,:,::1] DistanceTemplate3D
    cdef public float[:,::1] MeansArray2D
    cdef public float[:,::1] SDArray2D
    cdef public char[:,::1] KnownUnfillable2D
    cdef public char[:,::1] Coastline2D
    cdef public Py_ssize_t FillFromZPosition

    def __init__(self):
        pass

cdef class A2DataStack:
    cdef public float[:,::1] DataArray2D
    cdef public unsigned char[:,::1] FlagsArray2D
    cdef public float[:,::1] DistanceArray2D
    cdef public float[:,::1] MeansArray2D

    def __init__(self):
        pass


cdef class A2PassData:
    """Holds the images necessary to run the a2_core algorithm, in the correct direction
    The core algorithm iterates through the data in c-native order for efficiency. To implement the 8 directional
    passes, we transpose the data we pass to it, such that this ordering over the transposed pixels is equivalent to
    the required direction over the source pixels.
    e.g. if we want to iterate right to left, we flip the data first, then iterating left to right as normal is
    equivalent to iterating right to left in the original data.
    In the original Jupyter Notebook implementation, we did this by just re-striding the arrays. This has the desired
    behaviour regarding data order, but because the data are not "really" flipped in memory, access is slow, and those
    passes of the A2 algorithm ran ~ 9* slower. So here we physically copy the arrays into the new order,
    then copy them back again for output."""

    # all members are invisible to python, we want python code to only see things that are untransposed and can be
    # iterated in c-normal order, so use the getter methods to access the outputs
    cdef readonly float[:,::1] TransformedDataArray2D
    cdef readonly unsigned char[:,::1] TransformedFlagsArray2D
    cdef readonly float[:,::1] TransformedDistanceArray2D
    cdef readonly float[:,::1] TransformedMeanArray2D
    cdef public float[:,::1] TransformedDataArrayOutput2D
    cdef readonly float[:,::1] TransformedSumOfPassDistancesArray2D

    cdef char passNumber
    cdef float[:,::1] outputData

    def __cinit__(self, passNumber, dataArray, flagsArray, distanceArray, meanArray, sumDistArray):

        self.passNumber = passNumber

        self.TransformedDataArray2D = self.__transformArray__(dataArray)
        self.TransformedFlagsArray2D = self.__transformArray__(flagsArray)
        self.TransformedDistanceArray2D = self.__transformArray__(distanceArray)
        self.TransformedMeanArray2D = self.__transformArray__(meanArray)
        self.TransformedSumOfPassDistancesArray2D = self.__transformArray__(sumDistArray)

        self.TransformedDataArrayOutput2D = np.empty_like(self.TransformedDataArray2D)

    def getOutputData(self):
        """Returns the filled data array held in this object, transformed back to c-normal order"""
        return self.__transformArray__(self.TransformedDataArray2D)

    def getOutputDists(self):
        """Returns the sum-of-fill-distances array held in this object transformed back to c-normal order"""
        return self.__transformArray__(self.TransformedDistanceArray2D)

    def __transformArray__(self, data):
        """Transforms the input data arrays such that iterating over them in C-normal order will effectively
        pass over the data in one of the 8 cardinal directions, according to which pass number this is.
        This implementation actually copies the data into a new C-normal array, so that all A2 pass runs
        should take the same amount of time. Not copying the data i.e. not using np.copy results in the
        transposed passed 1-4 being ~ 9* slower due to the inefficient memory access.
        All transforms are reversible so can use the same call for input and output"""
        if self.passNumber == 0:
            return np.copy(data.T, order='C')
        elif self.passNumber == 1:
            return np.copy(data[:,::-1].T, order='C')
        elif self.passNumber == 2:
            return np.copy(data[::-1,:].T, order='C')
        elif self.passNumber == 3:
            return np.copy(data[::-1,::-1].T, order='C')
        elif self.passNumber == 4:
            return np.copy(data, order='C')
        elif self.passNumber == 5:
            return np.copy(data[:,::-1], order='C')
        elif self.passNumber == 6:
            return np.copy(data[::-1,:], order='C')
        elif self.passNumber ==7:
            return np.copy(data[::-1,::-1], order='C')
        else:
            raise ValueError()


cdef class PixelMargins:
    """ Simple object to hold 4 values representing margins for processing data """
    cdef readonly Py_ssize_t Top, Bottom, Left, Right

    def __init__(self, top=0, bottom=0, left=0, right=0):
        self.Top = top
        self.Bottom = bottom
        self.Left = left
        self.Right = right