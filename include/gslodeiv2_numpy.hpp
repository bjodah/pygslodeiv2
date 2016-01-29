#ifndef GSLODEIV2_H_RVEFMVJH25BGNEJDX4WHGA73YI
#define GSLODEIV2_H_RVEFMVJH25BGNEJDX4WHGA73YI
#include <Python.h>
#include <numpy/arrayobject.h>

#include <chrono> // std::clock
#include <utility> // std::pair
#include <vector> // std::vector

// based on GSL v1.16
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>


namespace gslodeiv2{
    const gsl_odeiv2_step_type * get_step_type(int index){
        switch(index){
        case 0:
            return gsl_odeiv2_step_rk2;
        case 1:
            return gsl_odeiv2_step_rk4;
        case 2:
            return gsl_odeiv2_step_rkf45;
        case 3:
            return gsl_odeiv2_step_rkck;
        case 4:
            return gsl_odeiv2_step_rk8pd;
        case 5:
            return gsl_odeiv2_step_rk1imp;
        case 6:
            return gsl_odeiv2_step_rk2imp;
        case 7:
            return gsl_odeiv2_step_rk4imp;
        case 8:
            return gsl_odeiv2_step_bsimp;
        case 9:
            return gsl_odeiv2_step_msadams;
        case 10:
            return gsl_odeiv2_step_msbdf;
        default:
            throw std::logic_error("Unknown steptype index");
        }
    }
    // typedef int (*RhsFn)(double t, const double y[], double dydt[], void *params);
    // typedef int (*JacFn)(double t, const double y[], double *dfdy, double dfdt[], void *params);
    int rhs(double t, const double y[], double f[], void * params);
    int jac(double x, const double y[], double *dfdy, double dfdt[], void *params);

    struct System : gsl_odeiv2_system {
        gsl_odeiv2_system m_sys;
        System(size_t ny, void *params) : m_sys({rhs, jac, ny, params}) {}
    };
    struct Driver {
        gsl_odeiv2_driver *m_driver;
        Driver(const int step_type_idx, const size_t dim,
               const System& sys, const double hstart,
               const double epsabs, const double epsrel //,
               ) : m_driver(gsl_odeiv2_driver_alloc_y_new(&sys.m_sys, get_step_type(step_type_idx),
                                                          hstart, epsabs, epsrel))
               //               const double a_y=1, const double a_dydt=0,
               //const double scale_abs[]=nullptr) :
            // m_driver((scale_abs == nullptr) ?
            //          gsl_odeiv2_driver_alloc_standard_new(&sys.m_sys, get_step_type(step_type_idx),
            //                                               hstart, epsabs, epsrel, a_y, a_dydt) :
            //          gsl_odeiv2_driver_alloc_scaled_new(&sys.m_sys, get_step_type(step_type_idx),
            //                                             hstart, epsabs, epsrel, a_y, a_dydt,
            //                                             scale_abs))
               {}
        ~Driver() { gsl_odeiv2_driver_free(this->m_driver); }
        int set_hmin(double hmin) { return gsl_odeiv2_driver_set_hmin(m_driver, hmin); }
        int set_hmax(double hmax) { return gsl_odeiv2_driver_set_hmax(m_driver, hmax); }
        int set_nmax(const unsigned long int nmax) { return gsl_odeiv2_driver_set_nmax(m_driver, nmax); }
        int apply(double * t, const double t1, double y[]){
            return gsl_odeiv2_driver_apply(m_driver, t, t1, y);
        }
        int reset() { return gsl_odeiv2_driver_reset(m_driver); }
        int reset_hstart(const double hstart) { return gsl_odeiv2_driver_reset_hstart(m_driver, hstart); }
    };
    struct Step {
        gsl_odeiv2_step *m_step;
        Step(const int step_type_idx, const size_t dim) :
            m_step(gsl_odeiv2_step_alloc(get_step_type(step_type_idx), dim)) {}
        ~Step() { gsl_odeiv2_step_free(m_step); }
        int reset() { return gsl_odeiv2_step_reset(m_step); }
        unsigned int order() { return gsl_odeiv2_step_order(m_step); }
        int set_driver(Driver& d) { return gsl_odeiv2_step_set_driver(m_step, d.m_driver); }
        int apply(double t, double h, double y[], double yerr[], const double dydt_in[],
                  double dydt_out[], const System& sys) {
            return gsl_odeiv2_step_apply(m_step, t, h, y, yerr, dydt_in, dydt_out, &sys.m_sys);
        }
    };
    struct Control {
        gsl_odeiv2_control *m_control;
        Control(double eps_abs, double eps_rel, double a_y=1, double a_dydt=0,
                const double scale_abs[]=nullptr, size_t dim=0) :
            m_control((scale_abs == nullptr) ?
                      gsl_odeiv2_control_standard_new(eps_abs, eps_rel, a_y, a_dydt) :
                      gsl_odeiv2_control_scaled_new(eps_abs, eps_rel, a_y,
                                                    a_dydt, scale_abs, dim)) {}
        ~Control() { gsl_odeiv2_control_free(m_control); }
        int set_driver(Driver& d) { return gsl_odeiv2_control_set_driver(m_control, d.m_driver); }
        int init(Control& c, double eps_abs, double eps_rel, double a_y, double a_dydt){
            return gsl_odeiv2_control_init(m_control, eps_abs, eps_rel, a_y, a_dydt);
        }
        int hadjust(Step& s, const double y[], const double yerr[], const double dydt[], double *h){
            return gsl_odeiv2_control_hadjust(m_control, s.m_step, y, yerr, dydt, h);
        }
    };
    struct Evolve {
        gsl_odeiv2_evolve *m_evolve;
        Evolve(size_t dim) : m_evolve(gsl_odeiv2_evolve_alloc(dim)) {}
        ~Evolve() { gsl_odeiv2_evolve_free(this->m_evolve); }
        int set_driver(Driver& d) { return gsl_odeiv2_evolve_set_driver(m_evolve, d.m_driver); }
        int apply(Control& con, Step& step, System& sys, double *t,
                  double t1, double *h, double y[]) {
            return gsl_odeiv2_evolve_apply(m_evolve, con.m_control, step.m_step,
                                           &sys.m_sys, t, t1, h, y);
        }
        int reset() { return gsl_odeiv2_evolve_reset(m_evolve); }
    };
    void py_handle_error(const char * reason,
                         const char * file,
                         int line,
                         int gsl_errno) {
        PyErr_SetString(PyExc_RuntimeError, reason);
    }

    struct PyErrorHandler {
        gsl_error_handler_t *m_ori_handler;
        PyErrorHandler() {
            m_ori_handler = gsl_set_error_handler(&py_handle_error);
        }
        ~PyErrorHandler() {
            if (m_ori_handler)
                this->release();
        }
        void release() {
            gsl_set_error_handler(m_ori_handler);
        }
    };

    class PyGslOdeiv2 {
    public:
        PyObject *py_rhs, *py_jac;
        size_t ny;
        size_t nrhs, njac;
        double time_cpu;
        std::vector<double> xout;
        std::vector<double> yout;

        PyGslOdeiv2(PyObject * py_rhs, PyObject * py_jac, size_t ny) :
            py_rhs(py_rhs), py_jac(py_jac), ny(ny) {}
        size_t adaptive(PyObject *py_y0, double x0, double xend,
                        double atol, double rtol, int step_type_idx,
                        double dx0, double dx_min=.0, double dx_max=.0, int nsteps=0){
            std::clock_t cputime0 = std::clock();
            auto handler = PyErrorHandler();
            System sys(this->ny, static_cast<void*>(this));
            Driver drv(step_type_idx, this->ny, sys, dx0, atol, rtol);
            Step step(step_type_idx, this->ny);
            step.set_driver(drv);
            Control control(atol, rtol);
            control.set_driver(drv);
            Evolve evolve(this->ny);
            evolve.set_driver(drv);
            int info;
            size_t steps_taken = 0;
            if (dx_min > 0) drv.set_hmin(dx_min);
            if (dx_max > 0) drv.set_hmax(dx_max);
            drv.set_nmax(nsteps);
            nrhs = 0; njac = 0;
            while (x0 < xend){
                info = evolve.apply(control, step, sys, &x0, xend, &dx0,
                                    (double*)PyArray_GETPTR1(py_y0, 0));
                if (info != GSL_SUCCESS)
                    throw std::runtime_error("evolve.apply failed");
                this->xout.push_back(x0);
                for (size_t i=0; i<this->ny; ++i)
                    this->yout.push_back(*(double*)PyArray_GETPTR1(py_y0, i));
                steps_taken++;
            }
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
            return steps_taken;
        }
        void predefined(PyObject *py_xout, PyObject *py_yout,
                        double atol, double rtol, int step_type_idx,
                        double dx0, double dx_min=0.0, double dx_max=0.0, int nsteps=0)
        {
            std::clock_t cputime0 = std::clock();
            auto handler = PyErrorHandler();
            System sys(this->ny, static_cast<void*>(this));
            Driver drv(step_type_idx, this->ny, sys, dx0, atol, rtol);
            int info;
            double xval = *(double*)PyArray_GETPTR1(py_xout, 0);
            npy_intp n_steps = PyArray_DIMS(py_xout)[0] - 1;
            if (dx_min > 0) drv.set_hmin(dx_min);
            if (dx_max > 0) drv.set_hmax(dx_max);
            drv.set_nmax(nsteps);
            nrhs = 0; njac = 0;
            for (npy_intp ix=0; ix<n_steps; ++ix){
                double * prev_ydata = (double *)PyArray_GETPTR2(py_yout, ix, 0);
                std::copy(prev_ydata, prev_ydata+this->ny,
                          (double*)PyArray_GETPTR2(py_yout, ix+1, 0));
                info = drv.apply(&xval, *(double*)PyArray_GETPTR1(py_xout, ix+1),
                                 (double*)PyArray_GETPTR2(py_yout, ix+1, 0));
                if (info != GSL_SUCCESS)
                    throw std::runtime_error("driver.apply failed");
            }
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
        }
    };

} // namespace gslodeiv2
int gslodeiv2::rhs(double xval, const double y[], double dydx[], void * params)
{
    auto obj = static_cast<gslodeiv2::PyGslOdeiv2*>(params);
    npy_intp dims[1] { static_cast<npy_intp>(obj->ny) } ;
    npy_intp y_strides[1] { sizeof(double) };
    PyObject* py_yarr = PyArray_New(
                &PyArray_Type, 1, dims, NPY_DOUBLE, y_strides,
                static_cast<void *>(const_cast<double *>(y)), sizeof(double),
                NPY_ARRAY_C_CONTIGUOUS, nullptr);
    PyObject * py_dydx = PyArray_SimpleNewFromData(
        1, dims, NPY_DOUBLE, static_cast<void*>(dydx));
    PyObject * py_arglist = Py_BuildValue("(dOO)", xval, py_yarr, py_dydx);
    PyObject * py_result = PyEval_CallObject(obj->py_rhs, py_arglist);
    Py_DECREF(py_arglist);
    Py_DECREF(py_dydx);
    Py_DECREF(py_yarr);
    obj->nrhs++;
    if (py_result == nullptr){
        return GSL_EBADFUNC;
    } else if (py_result != Py_None){
        // py_result is not None
        PyErr_SetString(PyExc_RuntimeError, "rhs() did not return None");
        Py_DECREF(py_result);
        return GSL_EBADFUNC;
    }
    Py_DECREF(py_result);
    return GSL_SUCCESS;
}

int gslodeiv2::jac(double xval, const double y[], double *dfdy, double dfdx[], void *params)
{
    auto obj = static_cast<gslodeiv2::PyGslOdeiv2*>(params);
    npy_intp ydims[1] { static_cast<npy_intp>(obj->ny) };
    npy_intp y_strides[1] { sizeof(double) };
    npy_intp Jdims[2] { static_cast<npy_intp>(obj->ny), static_cast<npy_intp>(obj->ny) };
    PyObject* py_yarr = PyArray_New(
                &PyArray_Type, 1, ydims, NPY_DOUBLE, y_strides,
                static_cast<void *>(const_cast<double *>(y)), sizeof(double),
                NPY_ARRAY_C_CONTIGUOUS, nullptr);
    PyObject * py_jmat = PyArray_SimpleNewFromData(2, Jdims, NPY_DOUBLE, const_cast<double *>(dfdy));
    PyObject * py_dfdx = PyArray_SimpleNewFromData(1, ydims, NPY_DOUBLE, const_cast<double *>(dfdx));
    PyObject * py_arglist = Py_BuildValue("(dOOO)", xval, py_yarr, py_jmat, py_dfdx);
    PyObject * py_result = PyEval_CallObject(obj->py_jac, py_arglist);
    Py_DECREF(py_arglist);
    Py_DECREF(py_dfdx);
    Py_DECREF(py_jmat);
    Py_DECREF(py_yarr);
    obj->njac++;
    if (py_result == nullptr){
        return GSL_EBADFUNC;
    } else if (py_result != Py_None){
        // py_result is not None
        PyErr_SetString(PyExc_RuntimeError, "jac() did not return None");
        Py_DECREF(py_result);
        return GSL_EBADFUNC;
    }
    Py_DECREF(py_result);
    return GSL_SUCCESS;
}

#endif /* GSLODEIV2_H_RVEFMVJH25BGNEJDX4WHGA73YI */
