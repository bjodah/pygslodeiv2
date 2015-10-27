# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
cimport numpy as cnp
import numpy as np

from gslodeiv2_numpy cimport PyGslOdeiv2

cnp.import_array()  # Numpy C-API initialization

cdef class GslOdeiv2:

    cdef PyGslOdeiv2 *thisptr
    cdef object rhs

    def __cinit__(self, object f, object j, size_t ny):
        self.rhs = f
        self.thisptr = new PyGslOdeiv2(<PyObject *>f, <PyObject *>j, ny)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                 double t0, double tend, double atol, double rtol,
                 int step_type_idx=8, double dx0=.0, double dx_min=.0,
                 double dx_max=.0):
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        if dx0 <= 0.0:
            raise ValueError("dx must be > 0")
        return self.thisptr.adaptive(<PyObject*>y0, t0, tend, atol,
                                     rtol, step_type_idx, dx0, dx_min, dx_max)

    def predefined(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                   cnp.ndarray[cnp.float64_t, ndim=1] xout,
                   double atol, double rtol, int step_type_idx=8, double dx0=.0,
                   double dx_max=0, double dx_min=0):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty(
            (xout.size, y0.size), dtype=np.float64)
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        if dx0 <= 0.0:
            raise ValueError("dx must be > 0")
        yout[0, :] = y0
        self.thisptr.predefined(<PyObject*>xout, <PyObject*> yout, atol,
                                rtol, step_type_idx, dx0, dx_max, dx_min)
        return yout

    def get_xout(self, size_t nsteps):
        cdef cnp.ndarray[cnp.float64_t, ndim=1] xout = np.empty(
            nsteps, dtype=np.float64)
        cdef size_t i
        for i in range(nsteps):
            xout[i] = self.thisptr.xout[i]
        return xout

    def get_yout(self, size_t nsteps):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty((
            nsteps, self.thisptr.ny), dtype=np.float64)
        cdef size_t i
        cdef size_t ny = self.thisptr.ny
        for i in range(nsteps):
            for j in range(ny):
                yout[i, j] = self.thisptr.yout[i*ny + j]
        return yout

    def post_process(self, cnp.ndarray[cnp.float64_t, ndim=1] xout,
                     cnp.ndarray[cnp.float64_t, ndim=2] yout, int nderiv=0):
        cdef size_t nx = yout.shape[0]
        cdef size_t ny = yout.shape[1]
        cdef size_t idx_x, idx_y
        cdef cnp.ndarray[cnp.float64_t, ndim=1] fout
        cdef cnp.ndarray[cnp.float64_t, ndim=3] yout2
        if nderiv == 0:
            return yout
        elif nderiv == 1:
            yout2 = np.empty((nx, 2, ny))
            yout2[:, 0, :] = yout
            fout = np.empty(ny)
            for idx_x in range(nx):
                self.rhs(xout[idx_x], yout[idx_x, :], fout)
                yout2[idx_x, 1, :] = fout
            return yout2
        else:
            raise NotImplementedError("Higher derivatives not available")

    def get_info(self):
        return {'nrhs': self.thisptr.nrhs, 'njac': self.thisptr.njac}


steppers = (
    'rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd', 'rk1imp',
    'rk2imp', 'rk4imp', 'bsimp', 'msadams', 'msbdf'
)
requires_jac = (
    'rk1imp', 'rk2imp', 'rk4imp', 'bsimp', 'msbdf'
)


def adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol, dx_min=.0, dx_max=.0,
             nderiv=0, method='bsimp'):
    cdef size_t nsteps
    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = GslOdeiv2(rhs, jac, len(y0))
    nsteps = integr.adaptive(np.array(y0, dtype=np.float64),
                             x0, xend, atol, rtol,
                             steppers.index(method),
                             dx0, dx_min, dx_max)
    xout = integr.get_xout(nsteps)
    return xout, integr.post_process(xout, integr.get_yout(nsteps), nderiv), integr.get_info()


def predefined(rhs, jac, y0, xout, dx0, atol, rtol, dx_min=.0, dx_max=.0,
               nderiv=0, method='bsimp'):
    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = GslOdeiv2(rhs, jac, len(y0))
    xout = np.asarray(xout, dtype=np.float64)
    yout = integr.predefined(np.asarray(y0, dtype=np.float64),
                             xout,
                             atol, rtol, steppers.index(method),
                             dx0, dx_min, dx_max)
    return integr.post_process(xout, yout, nderiv), integr.get_info()
