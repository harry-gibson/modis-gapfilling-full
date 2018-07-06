cimport cython

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

cdef class A2Data:
    cdef float[:,::1] DataArray2D
    cdef unsigned char[:,::1] FlagsArray2D
    cdef float[:,::1] DistanceArray2D

    def __init__(self, dataArray, flagsArray, distanceArray):
        self.DataArray2D = dataArray
        self.FlagsArray2D = flagsArray
        self.DistanceArray2D = distanceArray



cdef class PixelMargins:
    """ Simple object to hold 4 values representing margins for processing data """
    cdef readonly Py_ssize_t Top, Bottom, Left, Right

    def __init__(self, top=0, bottom=0, left=0, right=0):
        self.Top = top
        self.Bottom = bottom
        self.Left = left
        self.Right = right