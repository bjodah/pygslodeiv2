# -*- coding: utf-8 -*-
"""
Python binding for odeiv2 in GNU Scientific Library (GSL).
"""

from __future__ import division, absolute_import

import numpy as np

from ._gsl_odeiv2 import adaptive, predefined, requires_jac, steppers, fpes
from ._util import _check_callable, _check_indexing, _ensure_5args
from ._release import __version__


def get_include():
    from pkg_resources import resource_filename, Requirement
    return resource_filename(Requirement.parse(__name__),
                             '%s/include' % __name__)


def integrate_adaptive(rhs, jac, y0, x0, xend, atol, rtol, dx0=.0,
                       dx_min=.0, dx_max=.0, method='bsimp', nsteps=500,
                       check_callable=False, check_indexing=False,
                       autorestart=0, return_on_error=False, cb_kwargs=None, **kwargs):
    """ Integrates a system of ordinary differential equations (solver chosen output).

    Parameters
    ----------
    rhs : callable
        Function with signature f(t, y, fout) which modifies fout *inplace*.
    jac : callable
        Function with signature j(t, y, jmat_out, dfdx_out) which modifies
        jmat_out and dfdx_out *inplace*.
    y0 : array_like
        initial values of the dependent variables
    x0 : float
        initial value of the independent variable
    xend : float
        stopping value for the independent variable
    atol : float
        absolute tolerance
    rtol : float
        relative tolerance
    dx0 : float
        initial step-size
    dx_min : float
        minimum step (default: 0.0)
    dx_max : float
        maximum step (default: 0.0)
    method : str
        One of: 'rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd', 'rk1imp',
            'rk2imp', 'rk4imp', 'bsimp', 'msadams', 'msbdf'
    nsteps : int
        maximum number of steps (default: 500)
    check_callable : bool (default: False)
        perform signature sanity checks on ``rhs`` and ``jac``
    check_indexing : bool (default: False)
        perform item setting sanity checks on ``rhs`` and ``jac``.
    autorestart : int
        Autorestarts on error (requires autonomous system).
    return_on_error : bool
        Instead of raising an exception return silently (see info['success']).
    record_rhs_xvals : bool
        When True: will return x values for rhs calls in ``info['rhs_xvals']``.
    record_jac_xvals : bool
        When True will return x values for jac calls in ``info['jac_xvals']``.
    record_order : bool
        When True will return used time stepper order in ``info['orders']``.
    record_fpe : bool
        When True will return observed floating point errors in ``info['fpes']``. (see ``fpes``)
    cb_kwargs: dict
        Extra keyword arguments passed to ``rhs``, ``jac`` and possibly ``dx0cb``.
    dx0cb : callable
        Callback for calculating dx0 (make sure to pass dx0==0.0) to enable.
        Signature: ``f(x, y[:]) -> float``.
    dx_max_cb: callable
        Callback for calculating dx_max. Signature: ``f(x, y[:]) -> float``.

    Returns
    -------
    (xout, yout, info):
        xout: 1-dimensional array of values for the independent variable
        yout: 2-dimensional array of the dependent variables (axis 1) for
            values corresponding to xout (axis 0)
        info: dictionary with information about the integration

    """
    # Sanity checks to reduce risk of having a segfault:
    jac = _ensure_5args(jac)
    if check_callable:
        _check_callable(rhs, jac, x0, y0)

    if check_indexing:
        _check_indexing(rhs, jac, x0, y0)

    return adaptive(rhs, jac, np.ascontiguousarray(y0, dtype=np.float64), x0,
                    xend, atol, rtol, method, nsteps, dx0, dx_min, dx_max,
                    autorestart, return_on_error, cb_kwargs, **kwargs)


def integrate_predefined(rhs, jac, y0, xout, atol, rtol, dx0=.0,
                         dx_min=.0, dx_max=.0, method='bsimp', nsteps=500,
                         check_callable=False, check_indexing=False,
                         autorestart=0, return_on_error=False, cb_kwargs=None, **kwargs):
    """ Integrates a system of ordinary differential equations (user chosen output).

    Parameters
    ----------
    rhs : callable
        Function with signature f(t, y, fout) which modifies fout *inplace*.
    jac : callable
        Function with signature j(t, y, jmat_out, dfdx_out) which modifies
        jmat_out and dfdx_out *inplace*.
    y0 : array_like
        initial values of the dependent variables
    xout : array_like
        values of the independent variable
    atol : float
        absolute tolerance
    rtol : float
        relative tolerance
    dx0 : float
        initial step-size
    dx_min : float
        minimum step (default: 0.0)
    dx_max : float
        maximum step (default: 0.0)
    method : str
        One of: 'rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd', 'rk1imp',
        'rk2imp', 'rk4imp', 'bsimp', 'msadams', 'msbdf'.
    nsteps : int
        maximum number of steps (default: 500).
    check_callable : bool
        perform signature sanity checks on ``rhs`` and ``jac``.
    check_indexing : bool
        perform item setting sanity checks on ``rhs`` and ``jac``.
    autorestart : int
        Autorestarts on error (requires autonomous system).
    return_on_error : bool
        Instead of raising an exception return silently (see info['success'] & info['nreached']).
    record_rhs_xvals : bool
        When True: will return x values for rhs calls in ``info['rhs_xvals']``.
    record_jac_xvals : bool
        When True will return x values for jac calls in ``info['jac_xvals']``.
    record_order : bool
        When True will return used time stepper order in ``info['orders']``.
    record_fpe : bool
        When True will return observed floating point errors in ``info['fpes']``. (see ``fpes``)
    cb_kwargs : dict
        Extra keyword arguments passed to ``rhs`` and ``jac``.
    dx0cb : callable
        Callback for calculating dx0 (make sure to pass dx0==0.0) to enable.
        Signature: ``f(x, y[:]) -> float``.
    dx_max_cb: callable
        Callback for calculating dx_max. Signature: ``f(x, y[:]) -> float``.

    Returns
    -------
    (result, info):
        result: 2-dimensional array of the dependent variables (axis 1) for
            values corresponding to xout (axis 0)
        info: dictionary with information about the integration

    """
    # Sanity checks to reduce risk of having a segfault:
    jac = _ensure_5args(jac)
    if check_callable:
        _check_callable(rhs, jac, xout[0], y0)

    if check_indexing:
        _check_indexing(rhs, jac, xout[0], y0)

    return predefined(rhs, jac,
                      np.ascontiguousarray(y0, dtype=np.float64),
                      np.ascontiguousarray(xout, dtype=np.float64), atol, rtol,
                      method, nsteps, dx0, dx_min, dx_max,
                      autorestart, return_on_error, cb_kwargs, **kwargs)
