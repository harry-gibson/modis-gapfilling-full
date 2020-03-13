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

        self.TransformedDataArray2D = dataArray
        self.TransformedFlagsArray2D = flagsArray
        self.TransformedDistanceArray2D = distanceArray
        self.TransformedMeanArray2D = meanArray
        self.TransformedSumOfPassDistancesArray2D = sumDistArray

        self.passNumber = passNumber
        #self.dataBits=[self.DataArray2D, self.FlagsArray2D, self.DistanceArray2D, self.MeanArray2D, self.SumOfPassDistancesArray2D]

        self.__TransformData()

        self.TransformedDataArrayOutput2D = np.empty_like(dataArray)

    def getOutputData(self):
        """Returns the filled data array held in this object, transformed back to c-normal order"""
        if self.passNumber == 0:
            return np.copy(self.TransformedDataArrayOutput2D.T)
        elif self.passNumber == 1:
            return np.copy(self.TransformedDataArrayOutput2D[:,::-1].T)
        elif self.passNumber == 2:
            return np.copy(self.TransformedDataArrayOutput2D[::-1,:].T)
        elif self.passNumber == 3:
            return np.copy(self.TransformedDataArrayOutput2D[::-1,::-1].T)
        elif self.passNumber == 4:
            return np.copy(self.TransformedDataArrayOutput2D)
        elif self.passNumber == 5:
            return np.copy(self.TransformedDataArrayOutput2D[:,::-1])
        elif self.passNumber == 6:
            return np.copy(self.TransformedDataArrayOutput2D[::-1,:])
        elif self.passNumber ==7:
            return np.copy(self.TransformedDataArrayOutput2D[::-1,::-1])
        else:
            raise ValueError()

    def getOutputDists(self):
        """Returns the sum-of-fill-distances array held in this object transformed back to c-normal order"""
        if self.passNumber == 0:
            return np.copy(self.TransformedSumOfPassDistancesArray2D.T)
        elif self.passNumber == 1:
            return np.copy(self.TransformedSumOfPassDistancesArray2D[:,::-1].T)
        elif self.passNumber == 2:
            return np.copy(self.TransformedSumOfPassDistancesArray2D[::-1,:].T)
        elif self.passNumber == 3:
            return np.copy(self.TransformedSumOfPassDistancesArray2D[::-1,::-1].T)
        elif self.passNumber == 4:
            return np.copy(self.TransformedSumOfPassDistancesArray2D)
        elif self.passNumber == 5:
            return np.copy(self.TransformedSumOfPassDistancesArray2D[:,::-1])
        elif self.passNumber == 6:
            return np.copy(self.TransformedSumOfPassDistancesArray2D[::-1,:])
        elif self.passNumber ==7:
            return np.copy(self.TransformedSumOfPassDistancesArray2D[::-1,::-1])
        else:
            raise ValueError()

    def __TransformData(self):
        """Transforms the input data arrays such that iterating over them in C-normal order will effectively
        pass over the data in one of the 8 cardinal directions. This implementation actually copies the data into
        a new C-normal array, so that all A2 pass runs should take the same amount of time. Not copying the data
        i.e. not using np.copy results in the transposed passed 1-4 being ~ 9* slower due to the inefficient
        memory access"""
        dataBits = [self.TransformedDataArray2D, self.TransformedFlagsArray2D, self.TransformedDistanceArray2D,
                    self.TransformedMeanArray2D, self.TransformedSumOfPassDistancesArray2D]
        passNumber = self.passNumber
        for item in(self.TransformedDataArray2D, self.TransformedFlagsArray2D, self.TransformedDistanceArray2D,
                    self.TransformedMeanArray2D, self.TransformedSumOfPassDistancesArray2D):
            if passNumber == 0:
                item = np.copy(item.T)
            elif passNumber == 1:
                item = np.copy(item.T[:,::-1])
            elif passNumber == 2:
                item = np.copy(item.T[::-1,:])
            elif passNumber == 3:
                item = np.copy(item.T[::-1,::-1])
            elif passNumber == 4:
                pass
            elif passNumber == 5:
                item = np.copy(item[:,::-1])
            elif passNumber == 6:
                item = np.copy(item[::-1,:])
            elif passNumber == 7:
                item = np.copy(item[::-1,::-1])
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