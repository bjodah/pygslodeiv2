# -*- coding: utf-8 -*-

import numpy as np
import pytest

from pygslodeiv2 import (
    integrate_adaptive, integrate_predefined, requires_jac
)

decay_defaults = dict(y0=[0.7, 0.3, 0.5],
                      xout=np.linspace(0, 3, 31),
                      dx0=1e-10, atol=1e-8, rtol=1e-8)
decay_default_k = [2.0, 3.0, 4.0]

decay_analytic = {
    0: lambda y0, k, t: (
        y0[0] * np.exp(-k[0]*t)),
    1: lambda y0, k, t: (
        y0[1] * np.exp(-k[1] * t) + y0[0] * k[0] / (k[1] - k[0]) *
        (np.exp(-k[0]*t) - np.exp(-k[1]*t))),
    2: lambda y0, k, t: (
        y0[2] * np.exp(-k[2] * t) + y0[1] * k[1] / (k[2] - k[1]) *
        (np.exp(-k[1]*t) - np.exp(-k[2]*t)) +
        k[1] * k[0] * y0[0] / (k[1] - k[0]) *
        (1 / (k[2] - k[0]) * (np.exp(-k[0]*t) - np.exp(-k[2]*t)) -
         1 / (k[2] - k[1]) * (np.exp(-k[1]*t) - np.exp(-k[2]*t))))
}


def decay_get_Cref(k, y0, tout):
    coeffs = list(k) + [0]*(3-len(k))
    return np.column_stack([
        decay_analytic[i](y0, coeffs, tout) for i in range(
            min(3, len(k)+1))])


def _get_f_j(k):
    k0, k1, k2 = k

    def f(t, y, fout):
        fout[0] = -k0*y[0]
        fout[1] = k0*y[0] - k1*y[1]
        fout[2] = k1*y[1] - k2*y[2]

    def j(t, y, jmat_out, dfdx_out):
        jmat_out[0, 0] = -k0
        jmat_out[0, 1] = 0
        jmat_out[0, 2] = 0
        jmat_out[1, 0] = k0
        jmat_out[1, 1] = -k1
        jmat_out[1, 2] = 0
        jmat_out[2, 0] = 0
        jmat_out[2, 1] = k1
        jmat_out[2, 2] = -k2
        dfdx_out[0] = 0
        dfdx_out[1] = 0
        dfdx_out[2] = 0
    return f, j

methods = [('bsimp', 3e-4), ('msadams', 5), ('rkf45', 0.5),
           ('rkck', 0.3), ('rk8pd', 0.04), ('rk4imp', 0.8),
           ('msbdf', 23)]
# ['rk2', 'rk4', 'rk1imp', 'rk2imp']


@pytest.mark.parametrize("method,forgiveness", methods)
def test_integrate_adaptive(method, forgiveness):
    use_jac = method in requires_jac
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    f, j = _get_f_j(k)
    if not use_jac:
        j = None
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(x0=0, xend=3, dx0=1e-10, atol=atol, rtol=rtol,
                  method=method)
    # Run twice to catch possible side-effects:
    xout, yout, info = integrate_adaptive(f, j, y0, **kwargs)
    xout, yout, info = integrate_adaptive(f, j, y0, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref,
                       rtol=forgiveness*rtol,
                       atol=forgiveness*atol)
    assert info['nfev'] > 0
    if method in requires_jac:
        assert info['njev'] > 0

    integrate_adaptive(f, j, y0, nsteps=1, **kwargs)


@pytest.mark.parametrize("method,forgiveness", methods)
def test_integrate_predefined(method, forgiveness):
    use_jac = method in requires_jac
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    f, j = _get_f_j(k)
    if not use_jac:
        j = None
    xout = np.linspace(0, 3, 31)
    dx0 = 1e-10
    atol, rtol = 1e-8, 1e-8
    # Run twice to catch possible side-effects:
    yout, info = integrate_predefined(f, j, y0, xout, dx0, atol,
                                      rtol, method=method)
    yout, info = integrate_predefined(f, j, y0, xout, dx0, atol,
                                      rtol, method=method)
    yref = decay_get_Cref(k, y0, xout)
    print(yout)
    print(yref)
    assert np.allclose(yout, yref,
                       rtol=forgiveness*rtol,
                       atol=forgiveness*atol)
    assert info['nfev'] > 0
    if method in requires_jac:
        assert info['njev'] > 0


def test_bad_f():
    k0, k1, k2 = decay_default_k

    def f(t, y, fout):
        y[0] = -1  # read-only! should raise ValueError
        fout[0] = -k0*y[0]
        fout[1] = k0*y[0] - k1*y[1]
        fout[2] = k1*y[1] - k2*y[2]
    with pytest.raises(ValueError):
        yout, info = integrate_predefined(f, None, method='rkck',
                                          check_callable=False,
                                          check_indexing=False,
                                          **decay_defaults)
        assert yout  # silence pyflakes


def test_bad_j():
    k0, k1, k2 = decay_default_k

    def f(t, y, fout):
        fout[0] = -k0*y[0]
        fout[1] = k0*y[0] - k1*y[1]
        fout[2] = k1*y[1] - k2*y[2]

    def j(t, y, jmat_out, dfdx_out):
        y[0] = -1  # read-only! should raise ValueError
        jmat_out[0, 0] = -k0
        jmat_out[0, 1] = 0
        jmat_out[0, 2] = 0
        jmat_out[1, 0] = k0
        jmat_out[1, 1] = -k1
        jmat_out[1, 2] = 0
        jmat_out[2, 0] = 0
        jmat_out[2, 1] = k1
        jmat_out[2, 2] = -k2
        dfdx_out[0] = 0
        dfdx_out[1] = 0
        dfdx_out[2] = 0
    with pytest.raises(ValueError):
        yout, info = integrate_predefined(f, j, method='bsimp',
                                          check_callable=False,
                                          check_indexing=False,
                                          **decay_defaults)
        assert yout  # silence pyflakes
