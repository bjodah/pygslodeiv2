# -*- coding: utf-8 -*-
"""
Python binding for odeiv2 in GNU Scientific Library (GSL).
"""

from __future__ import division, absolute_import

import numpy as np

from ._gsl_odeiv2 import adaptive, predefined, requires_jac, steppers
from ._util import _check_callable, _check_indexing, _ensure_5args
from ._release import __version__


def get_include():
    from pkg_resources import resource_filename, Requirement
    return resource_filename(Requirement.parse(__name__),
                             '%s/include' % __name__)


def integrate_adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol,
                       dx_min=.0, dx_max=.0, method='bsimp', nsteps=500,
                       check_callable=False, check_indexing=False,
                       autorestart=0, return_on_error=False, cb_kwargs=None):
    """
    Integrates a system of ordinary differential equations.

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
    dx0 : float
        initial step-size
    atol : float
        absolute tolerance
    rtol : float
        relative tolerance
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
    cb_kwargs: dict
        Extra keyword arguments passed to ``rhs`` and ``jac``.
    autorestart : int
        Autorestarts on error (requires autonomous system).
    return_on_error : bool
        Instead of raising an exception return silently (see info['success']).
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
                    autorestart, return_on_error, cb_kwargs)


def integrate_predefined(rhs, jac, y0, xout, dx0, atol, rtol,
                         dx_min=.0, dx_max=.0, method='bsimp', nsteps=500,
                         check_callable=False, check_indexing=False,
                         cb_kwargs=None):
    """
    Integrates a system of ordinary differential equations.

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
    dx0 : float
        initial step-size
    atol : float
        absolute tolerance
    rtol : float
        relative tolerance
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
    cb_kwargs : dict
        Extra keyword arguments passed to ``rhs`` and ``jac``.

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
                      method, nsteps, dx0, dx_min, dx_max, cb_kwargs)
