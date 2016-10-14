# -*- coding: utf-8; mode: cython -*-

from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "gsl_odeiv2_cxx.hpp" namespace "gsl_odeiv2_cxx":
    cdef cppclass StepType:
        pass

    cdef StepType styp_from_name(string) nogil except +
    cdef bool requires_jacobian(StepType) nogil


cdef extern from "gsl_odeiv2_cxx.hpp" namespace "gsl_odeiv2_cxx::StepType":
    cdef StepType RK2
    cdef StepType RK4
    cdef StepType RKF45
    cdef StepType RKCK
    cdef StepType RK8PD
    cdef StepType RK1IMP
    cdef StepType RK2IMP
    cdef StepType RK4IMP
    cdef StepType BSIMP
    cdef StepType MSADAMS
    cdef StepType MSBDF
