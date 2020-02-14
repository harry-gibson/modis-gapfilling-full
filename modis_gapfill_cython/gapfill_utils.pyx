cimport cython
cimport numpy as np

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

cdef class A2PassData:
    # all members are invisible to python, we want python code to only see things that are untransposed and can be
    # iterated in c-normal order, so use the getter methods to access the outputs
    cdef float[:,::1] DataArray2D
    cdef unsigned char[:,::1] FlagsArray2D
    cdef float[:,::1] DistanceArray2D
    cdef float[:,::1] MeanArray2D
    cdef float[:,::1] DataArrayOutput2D
    cdef char passNumber
    cdef float[:,::1] outputData
    cdef float[:,::1] SumOfPassDistancesArray2D

    def __cinit__(self, passNumber, dataArray, flagsArray, distanceArray, meanArray, sumDistArray):

        self.DataArray2D = dataArray
        self.FlagsArray2D = flagsArray
        self.DistanceArray2D = distanceArray
        self.MeanArray2D = meanArray
        self.SumOfPassDistancesArray2D = sumDistArray

        self.passNumber = passNumber
        self.dataBits=[self.DataArray2D, self.FlagsArray2D, self.DistanceArray2D, self.MeanArray2D, self.SumOfPassDistancesArray2D]

        self.__TransformData()

        self.outputData = np.empty_like(dataArray)

    def getOutputData(self):
        """Returns the filled data array held in this object, transformed back to c-normal order"""
        if self.passNumber == 0:
            return np.copy(self.outputData.T)
        elif self.passNumber == 1:
            return np.copy(self.outputData[:,::-1].T)
        elif self.passNumber == 2:
            return np.copy(self.outputData[::-1,:].T)
        elif self.passNumber == 3:
            return np.copy(self.outputData[::-1,::-1].T)
        elif self.passNumber == 4:
            return np.copy(self.outputData)
        elif self.passNumber == 5:
            return np.copy(self.outputData[:,::-1])
        elif self.passNumber == 6:
            return np.copy(self.outputData[::-1,:])
        elif self.passNumber ==7:
            return np.copy(self.outputData[::-1,::-1])
        else:
            raise ValueError()

    def getOutputDists(self):
        """Returns the sum-of-fill-distances array held in this object transformed back to c-normal order"""
        if self.passNumber == 0:
            return np.copy(self.SumOfPassDistancesArray2D.T)
        elif self.passNumber == 1:
            return np.copy(self.SumOfPassDistancesArray2D[:,::-1].T)
        elif self.passNumber == 2:
            return np.copy(self.SumOfPassDistancesArray2D[::-1,:].T)
        elif self.passNumber == 3:
            return np.copy(self.SumOfPassDistancesArray2D[::-1,::-1].T)
        elif self.passNumber == 4:
            return np.copy(self.SumOfPassDistancesArray2D)
        elif self.passNumber == 5:
            return np.copy(self.SumOfPassDistancesArray2D[:,::-1])
        elif self.passNumber == 6:
            return np.copy(self.SumOfPassDistancesArray2D[::-1,:])
        elif self.passNumber ==7:
            return np.copy(self.SumOfPassDistancesArray2D[::-1,::-1])
        else:
            raise ValueError()

    def __TransformData(self, passNumber):
        """Transforms the input data arrays such that iterating over them in C-normal order will effectively
        pass over the data in one of the 8 cardinal directions. This implementation actually copies the data into
        a new C-normal array, so that all A2 pass runs should take the same amount of time. Not copying the data
        i.e. not using np.copy results in the transposed passed 1-4 being ~ 9* slower due to the inefficient
        memory access"""
        if passNumber == 0:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i].T)
        elif passNumber == 1:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i].T[:,::-1])
        elif passNumber == 1:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i].T[:,::-1])
        elif passNumber == 2:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i].T[::-1,:])
        elif passNumber == 3:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i].T[::-1,::-1])
        elif passNumber == 4:
            pass
        elif passNumber == 5:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i][:,::-1])

        elif passNumber == 6:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i][::-1,:])

        elif passNumber == 7:
            for i in range(len(self.dataBits)):
                self.dataBits[i] = np.copy(self.dataBits[i][::-1,::-1])
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