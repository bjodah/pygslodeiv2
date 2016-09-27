# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.string cimport string

cimport numpy as cnp
import numpy as np

from gsl_odeiv2_cxx cimport styp_from_name
from gsl_odeiv2_numpy cimport PyGslOdeiv2

cnp.import_array()  # Numpy C-API initialization

requires_jac = ('rk1imp', 'rk2imp', 'rk4imp', 'bsimp', 'msbdf')
steppers = requires_jac + ('rk2', 'rk4', 'rkf45', 'rkck', 'rk8pd', 'msadams')

cdef class GslOdeiv2:

    cdef PyGslOdeiv2 *thisptr
    cdef object rhs
    cdef bint success

    def __cinit__(self, object f, object j, size_t ny):
        self.rhs = f
        self.thisptr = new PyGslOdeiv2(<PyObject *>f, <PyObject *>j, ny)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                 double t0, double tend, double atol, double rtol,
                 str method='bsimp', double dx0=.0, double dx_min=.0,
                 double dx_max=.0, int nsteps=0):
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        if dx0 <= 0.0:
            raise ValueError("dx must be > 0")
        self.success = False
        try:
            xout, yout = self.thisptr.adaptive(<PyObject*>y0, t0, tend, atol, rtol,
                                               styp_from_name(method.upper().encode('UTF-8')),
                                               dx0, dx_min, dx_max, nsteps)
        except:
            raise
        else:
            self.success = True
        return np.asarray(xout), np.asarray(yout).reshape((len(xout), self.thisptr.ny))

    def predefined(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                   cnp.ndarray[cnp.float64_t, ndim=1] xout,
                   double atol, double rtol, str method='bsimp', double dx0=.0,
                   double dx_max=0, double dx_min=0, int nsteps=0):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty(
            (xout.size, y0.size), dtype=np.float64)
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        if dx0 <= 0.0:
            raise ValueError("dx must be > 0")
        yout[0, :] = y0
        self.success = False
        try:
            self.thisptr.predefined(<PyObject*>xout, <PyObject*> yout, atol,
                                    rtol, styp_from_name(method.upper().encode('UTF-8')),
                                    dx0, dx_min, dx_max, nsteps)
        except:
            raise
        else:
            self.success = True
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
        return {
            'success': self.success,
            'nfev': self.thisptr.nfev,
            'njev': self.thisptr.njev,
            'time_cpu': self.thisptr.time_cpu,
            'time_wall': self.thisptr.time_wall
        }


def adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol, dx_min=.0, dx_max=.0,
             nsteps=0, nderiv=0, method='bsimp'):
    cdef size_t steps_taken
    if method.lower() in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = GslOdeiv2(rhs, jac, len(y0))
    xout, yout = integr.adaptive(np.array(y0, dtype=np.float64),
                             x0, xend, atol, rtol, method,
                             dx0, dx_min, dx_max, nsteps)
    return xout, integr.post_process(xout, yout, nderiv), integr.get_info()


def predefined(rhs, jac, y0, xout, dx0, atol, rtol, dx_min=.0, dx_max=.0,
               nsteps=0, nderiv=0, method='bsimp'):
    if method.lower() in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = GslOdeiv2(rhs, jac, len(y0))
    xout = np.asarray(xout, dtype=np.float64)
    yout = integr.predefined(np.asarray(y0, dtype=np.float64),
                             xout,
                             atol, rtol, method,
                             dx0, dx_min, dx_max, nsteps)
    return integr.post_process(xout, yout, nderiv), integr.get_info()
