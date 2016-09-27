# -*- coding: utf-8; mode: cython -*-

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp cimport bool

cdef extern from "gsl_odeiv2_cxx.hpp" namespace "gsl_odeiv2_cxx":
    cdef cppclass StepType:
        pass

    cdef void simple_predefined[U](
        U * const, const double, const double, const StepType, const double * const,
        const size_t, const double * const, double * const, const double,
        const double, const double, const long int
    ) nogil except +

    cdef pair[vector[double], vector[double]] simple_adaptive[U](
        U * const,
        const double,
        const double,
        const StepType,
        const double * const,
        const double,
        const double,
        const double,
        const double,
        const double,
        const long int
    ) nogil except +

    cdef StepType styp_from_name(string) except +
    cdef bool requires_jacobian(StepType)


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
