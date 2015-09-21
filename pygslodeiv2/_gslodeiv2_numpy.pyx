# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
cimport numpy as cnp
import numpy as np

from gslodeiv2_numpy cimport PyGslOdeiv2

cnp.import_array()  # Numpy C-API initialization

cdef class GslOdeiv2:

    cdef PyGslOdeiv2 *thisptr

    def __cinit__(self, object f, object j, size_t ny):
        self.thisptr = new PyGslOdeiv2(<PyObject *>f, <PyObject *>j, ny)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                 double t0, double tend,
                 double atol, double rtol,
                 double hstart=0.0, int step_type_idx=8):
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        return self.thisptr.adaptive(<PyObject*>y0, t0, tend, atol,
                                     rtol, hstart, step_type_idx)

    def predefined(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                   cnp.ndarray[cnp.float64_t, ndim=1] xout,
                   double dx0, double atol, double rtol,
                   int step_type_idx=8, double dx_max=0, double dx_min=0):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty((xout.size, y0.size),
                                                                dtype=np.float64)
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        yout[0, :] = y0
        self.thisptr.predefined(<PyObject*>y0, <PyObject*>xout, <PyObject*> yout,
                                dx0, atol, rtol, step_type_idx, dx_max, dx_min)
        return yout

    def get_xout(self, size_t nsteps):
        cdef cnp.ndarray[cnp.float64_t, ndim=1] xout = np.empty(nsteps, dtype=np.float64)
        cdef int i
        for i in range(nsteps):
            xout[i] = self.thisptr.xout[i]
        return xout

    def get_yout(self, size_t nsteps):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty((nsteps, self.thisptr.ny),
                                                                dtype=np.float64)
        cdef int i
        cdef int ny = self.thisptr.ny
        for i in range(nsteps):
            for j in range(ny):
                yout[i, j] = self.thisptr.yout[i*ny + j]
        return yout


step_type_indices = ['rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd', 'rk1imp',
                     'rk2imp', 'rk4imp', 'bsimp', 'msadams', 'msbdf']
requires_jac = ('rk1imp', 'rk2imp', 'rk4imp', 'bsimp', 'msbdf')

def adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol, method='bsimp'):
    cdef size_t nsteps
    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = GslOdeiv2(rhs, jac, len(y0))
    nsteps = integr.adaptive(np.asarray(y0, dtype=np.float64),
                             x0, xend, dx0, atol, rtol,
                             step_type_indices.index(method))
    return integr.get_xout(nsteps), integr.get_yout(nsteps)


def predefined(rhs, jac, y0, xout, dx0, atol, rtol, method='bsimp'):
    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = GslOdeiv2(rhs, jac, len(y0))
    return integr.predefined(np.asarray(y0, dtype=np.float64),
                             np.asarray(xout, dtype=np.float64),
                             dx0, atol, rtol, step_type_indices.index(method))
