# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from gsl_odeiv2_cxx cimport StepType

cdef extern from "gsl_odeiv2_anyode_parallel.hpp" namespace "gsl_odeiv2_anyode_parallel":
    cdef vector[pair[vector[double], vector[double]]] multi_adaptive[U](
        vector[U*],
        double,
        double,
        StepType,
        const double * const,
        const double *,
        const double *,
        long int,
        double *,
        double *,
        double *,
        int,
        bool
    ) except + nogil

    cdef vector[int] multi_predefined[U](
        vector[U*],
        double,
        double,
        StepType,
        const double * const,
        size_t,
        const double * const,
        double * const,
        long int,
        double *,
        double *,
        double *,
        int,
        bool
    ) except + nogil
