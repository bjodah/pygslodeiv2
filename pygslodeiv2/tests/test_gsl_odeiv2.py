# -*- coding: utf-8 -*-

import os
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
    assert info['success']
    assert info['atol'] == atol and info['rtol'] == rtol
    assert np.allclose(yout, yref,
                       rtol=forgiveness*rtol,
                       atol=forgiveness*atol)
    assert info['nfev'] > 0
    if method in requires_jac:
        assert info['njev'] > 0

    with pytest.raises(RuntimeError) as excinfo:
        integrate_adaptive(f, j, y0, nsteps=7, **kwargs)
    assert 'steps' in str(excinfo.value).lower()
    assert '7' in str(excinfo.value).lower()


@pytest.mark.parametrize("method,forgiveness", methods)
def test_integrate_predefined(method, forgiveness):
    use_jac = method in requires_jac
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    f, j = _get_f_j(k)
    if not use_jac:
        j = None
    xout = np.linspace(0, 3, 31)
    kwargs = dict(atol=1e-8, rtol=1e-8, dx0=1e-10, method=method)
    atol, rtol = 1e-8, 1e-8
    # Run twice to catch possible side-effects:
    yout, info = integrate_predefined(f, j, y0, xout, **kwargs)
    yout, info = integrate_predefined(f, j, y0, xout, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert info['success']
    assert info['atol'] == atol and info['rtol'] == rtol
    assert np.allclose(yout, yref,
                       rtol=forgiveness*rtol,
                       atol=forgiveness*atol)
    assert info['nfev'] > 0
    if method in requires_jac:
        assert info['njev'] > 0
    if os.name == 'posix':
        assert info['time_cpu'] >= 0
        assert info['time_wall'] >= 0


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


def test_adaptive_return_on_error():
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(x0=0, xend=3, dx0=1e-10, atol=atol, rtol=rtol,
                  method='bsimp')
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, nsteps=7, return_on_error=True, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)
    assert xout.size == 8
    assert xout[-1] > 1e-6
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] is False
    assert xout[-1] < kwargs['xend']  # obviously not strict


def test_predefined_autorestart():
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    x0, xend = 0, 3
    kwargs = dict(dx0=1e-10, atol=atol, rtol=rtol,
                  method='bsimp', nsteps=62,
                  autorestart=10)
    f, j = _get_f_j(k)
    xout = np.linspace(x0, xend)
    yout, info = integrate_predefined(f, j, y0, xout, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref,
                       rtol=10*rtol,
                       atol=10*atol)
    assert xout[-1] > 1e-6
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success']
    assert xout[-1] == xend


def test_predefined_return_on_error():
    k = k0, k1, k2 = 2.0, 300.0, 4.0
    y0 = [0.7, 0.0, 0.2]
    atol, rtol = 1e-12, 1e-12
    kwargs = dict(dx0=1e-11, atol=atol, rtol=rtol, method='bsimp',
                  return_on_error=True, nsteps=10)
    f, j = _get_f_j(k)
    xout = np.logspace(-10, 1, 5)
    yout, info = integrate_predefined(f, j, y0, xout, **kwargs)
    yref = decay_get_Cref(k, y0, xout - xout[0])
    assert np.allclose(yout[:info['nreached'], :], yref[:info['nreached'], :],
                       rtol=10*rtol, atol=10*atol)
    assert 2 < info['nreached'] < 5  # not strict
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] is False


def test_dx0cb():
    k = 1e23, 3.0, 4.0
    y0 = [.7, .0, .0]
    x0, xend = 0, 5
    kwargs = dict(atol=1e-8, rtol=1e-8, method='bsimp', dx0cb=lambda x, y: y[0]*1e-30)
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, x0, xend, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref, atol=40*kwargs['atol'], rtol=40*kwargs['rtol'])
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] is True
    assert xout[-1] == xend


def test_dx_max_cb():
    k = 1e23, 3.0, 4.0
    y0 = [.7, .0, .0]
    x0, xend = 0, 5
    kwargs = dict(atol=1e-8, rtol=1e-8, method='bsimp', dx_max_cb=lambda x, y: 1e-3, nsteps=xend*1050)
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, x0, xend, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref, atol=40*kwargs['atol'], rtol=40*kwargs['rtol'])
    assert info['n_steps'] > 5000
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] is True
    assert xout[-1] == xend
