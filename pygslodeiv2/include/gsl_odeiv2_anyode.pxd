# -*- coding: utf-8; mode: cython -*-

from libcpp.vector cimport vector
from libcpp.utility cimport pair

from gsl_odeiv2_cxx cimport StepType

cdef extern from "gsl_odeiv2_anyode.hpp" namespace "gsl_odeiv2_anyode":

    cdef pair[vector[double], vector[double]] simple_adaptive[U](
        U * const,
        const double,
        const double,
        const StepType,
        const double * const,
        const double,
        const double,
        const long int,
        const double,
        const double,
        const double
    ) except +

    cdef void simple_predefined[U](
        U * const,
        const double,
        const double,
        const StepType,
        const double * const,
        const size_t,
        const double * const,
        double * const,
        const long int,
        const double,
        const double,
        const double
    ) except +
