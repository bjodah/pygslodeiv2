# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from gsl_odeiv2_cxx cimport StepType

cdef extern from "gsl_odeiv2_numpy.hpp" namespace "gsl_odeiv2_numpy":
    cdef cppclass PyErrorHandler:
        pass

    cdef cppclass PyGslOdeiv2:
        size_t ny, nfev, njev
        double time_cpu, time_wall

        PyGslOdeiv2(PyObject*, PyObject*, size_t)
        pair[vector[double], vector[double]] adaptive(PyObject*, double, double, double, double,
                                                      StepType, double, double, double, int) except +
        void predefined(PyObject*, PyObject*, double, double, StepType, double, double, double, int) except +
