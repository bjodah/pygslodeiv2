# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++

from cpython.object cimport PyObject
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport numpy as cnp

from anyode_numpy cimport PyOdeSys
from gsl_odeiv2_cxx cimport styp_from_name
from gsl_odeiv2_anyode cimport simple_adaptive, simple_predefined

import numpy as np

cnp.import_array()  # Numpy C-API initialization

requires_jac = ('rk1imp', 'rk2imp', 'rk4imp', 'bsimp', 'msbdf')
steppers = requires_jac + ('rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd', 'msadams')

cdef dict get_last_info(PyOdeSys * odesys, success=True):
    info = {str(k.decode('utf-8')): v for k, v in dict(odesys.last_integration_info).items()}
    info.update({str(k.decode('utf-8')): v for k, v in dict(odesys.last_integration_info_dbl).items()})
    info['nfev'] = odesys.nfev
    info['njev'] = odesys.njev
    info['success'] = success
    return info

def adaptive(rhs, jac, cnp.ndarray[cnp.float64_t, mode='c'] y0, double x0, double xend, double atol,
             double rtol, str method='bsimp', long int nsteps=500, double dx0=0.0, double dx_min=0.0,
             double dx_max=0.0, int autorestart=0, bool return_on_error=False, cb_kwargs=None):
    cdef:
        int ny = y0.shape[y0.ndim - 1]
        PyOdeSys * odesys

    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    if np.isnan(y0).any():
        raise ValueError("NaN found in y0")

    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, NULL,
                          <PyObject *>cb_kwargs, -1, -1, 0)
    try:
        xout, yout = map(np.asarray, simple_adaptive[PyOdeSys](
            odesys, atol, rtol, styp_from_name(method.lower().encode('UTF-8')),
            &y0[0], x0, xend, nsteps, dx0, dx_min, dx_max, autorestart, return_on_error))
        info = get_last_info(odesys, False if return_on_error and xout[-1] != xend else True)
        return xout, yout.reshape((xout.size, ny)), info
    finally:
        del odesys


def predefined(rhs, jac,
               cnp.ndarray[cnp.float64_t, mode='c'] y0,
               cnp.ndarray[cnp.float64_t, ndim=1] xout,
               double atol, double rtol, str method='bsimp', int nsteps=500, double dx0=0.0,
               double dx_min=0.0, double dx_max=0.0, cb_kwargs=None):
    cdef:
        int ny = y0.shape[y0.ndim - 1]
        cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty((xout.size, ny))
        PyOdeSys * odesys

    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    if np.isnan(y0).any():
        raise ValueError("NaN found in y0")
    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, NULL, <PyObject *>cb_kwargs, -1, -1, 0)
    try:
        simple_predefined[PyOdeSys](odesys, atol, rtol, styp_from_name(method.lower().encode('UTF-8')), &y0[0],
                                    xout.size, &xout[0], <double *>yout.data, nsteps,
                                    dx0, dx_min, dx_max)
        return yout.reshape((xout.size, ny)), get_last_info(odesys)
    finally:
        del odesys
