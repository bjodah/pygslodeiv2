# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.vector cimport vector

cdef extern from "gslodeiv2_numpy.hpp" namespace "gslodeiv2":
    cdef cppclass PyGslOdeiv2:
        size_t ny
        vector[double] xout
        vector[double] yout

        PyGslOdeiv2(PyObject*, PyObject*, size_t)
        size_t integrate_adaptive(PyObject*, double, double, double, double, double, int)
