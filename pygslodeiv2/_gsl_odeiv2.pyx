# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++

from cpython.object cimport PyObject
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport numpy as cnp

from anyode_numpy cimport PyOdeSys
from gsl_odeiv2_cxx cimport styp_from_name, fpes as _fpes
from gsl_odeiv2_anyode cimport simple_adaptive, simple_predefined

import numpy as np

cnp.import_array()  # Numpy C-API initialization

requires_jac = ('rk1imp', 'rk2imp', 'rk4imp', 'bsimp', 'msbdf')
steppers = requires_jac + ('rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd', 'msadams')
fpes = {str(k.decode('utf-8')): v for k, v in dict(_fpes).items()}

cdef dict get_last_info(PyOdeSys * odesys, success=True):
    info = {str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_int).items()}
    info.update({str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_dbl).items()})
    info.update({str(k.decode('utf-8')): np.array(v, dtype=np.float64) for k, v in dict(odesys.current_info.nfo_vecdbl).items()})
    info.update({str(k.decode('utf-8')): np.array(v, dtype=int) for k, v in dict(odesys.current_info.nfo_vecint).items()})
    info['nfev'] = odesys.nfev
    info['njev'] = odesys.njev
    info['success'] = success
    return info

def adaptive(rhs, jac, cnp.ndarray[cnp.float64_t, mode='c'] y0, double x0, double xend, double atol,
             double rtol, str method='bsimp', long int nsteps=500, double dx0=0.0, double dx_min=0.0,
             double dx_max=0.0, int autorestart=0, bool return_on_error=False, cb_kwargs=None,
             bool record_rhs_xvals=False, bool record_jac_xvals=False, bool record_order=False,
             bool record_fpe=False, dx0cb=None, dx_max_cb=None):
    cdef:
        int ny = y0.shape[y0.ndim - 1]
        PyOdeSys * odesys

    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    if np.isnan(y0).any():
        raise ValueError("NaN found in y0")

    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, NULL, NULL, NULL,
                          <PyObject *>cb_kwargs, -1, -1, 0, 0, <PyObject *>dx0cb, <PyObject *>dx_max_cb)
    odesys.record_rhs_xvals = record_rhs_xvals
    odesys.record_jac_xvals = record_jac_xvals
    odesys.record_order = record_order
    odesys.record_fpe = record_fpe

    try:
        xout, yout = map(np.asarray, simple_adaptive[PyOdeSys](
            odesys, atol, rtol, styp_from_name(method.lower().encode('UTF-8')),
            &y0[0], x0, xend, nsteps, dx0, dx_min, dx_max, autorestart, return_on_error))
        info = get_last_info(odesys, False if return_on_error and xout[-1] != xend else True)
        info['atol'], info['rtol'] = atol, rtol
        return xout, yout.reshape((xout.size, ny)), info
    finally:
        del odesys


def predefined(rhs, jac,
               cnp.ndarray[cnp.float64_t, mode='c'] y0,
               cnp.ndarray[cnp.float64_t, ndim=1] xout,
               double atol, double rtol, str method='bsimp', int nsteps=500, double dx0=0.0,
               double dx_min=0.0, double dx_max=0.0, int autorestart=0,
               bool return_on_error=False, cb_kwargs=None, bool record_rhs_xvals=False,
               bool record_jac_xvals=False, bool record_order=False, bool record_fpe=False,
               dx0cb=None, dx_max_cb=None):
    cdef:
        int ny = y0.shape[y0.ndim - 1]
        cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty((xout.size, ny))
        int nreached
        PyOdeSys * odesys

    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    if np.isnan(y0).any():
        raise ValueError("NaN found in y0")
    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, NULL, NULL, NULL, <PyObject *>cb_kwargs,
                          -1, -1, 0, 0, <PyObject *>dx0cb, <PyObject *>dx_max_cb)
    odesys.record_rhs_xvals = record_rhs_xvals
    odesys.record_jac_xvals = record_jac_xvals
    odesys.record_order = record_order
    odesys.record_fpe = record_fpe
    try:
        nreached = simple_predefined[PyOdeSys](odesys, atol, rtol, styp_from_name(method.lower().encode('UTF-8')), &y0[0],
                                               xout.size, &xout[0], <double *>yout.data, nsteps,
                                               dx0, dx_min, dx_max, autorestart, return_on_error)
        info = get_last_info(odesys, success=False if return_on_error and nreached < xout.size else True)
        info['nreached'] = nreached
        info['atol'], info['rtol'] = atol, rtol
        return yout.reshape((xout.size, ny)), info
    finally:
        del odesys
