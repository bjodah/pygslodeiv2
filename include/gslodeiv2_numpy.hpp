#ifndef GSLODEIV2_H_RVEFMVJH25BGNEJDX4WHGA73YI
#define GSLODEIV2_H_RVEFMVJH25BGNEJDX4WHGA73YI
#include <Python.h>
#include <numpy/arrayobject.h>

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
    // struct System : gsl_odeiv2_system {
    //     System() : gsl_odeiv2_system() {}
    // };
    using System = gsl_odeiv2_system;
    struct Driver {
        gsl_odeiv2_driver *m_driver;
        Driver(const int step_type_idx, const size_t dim,
               const System sys, const double hstart,
               const double epsabs, const double epsrel,
               const double a_y=1, const double a_dydt=0,
               const double scale_abs[]=nullptr) :
            m_driver((scale_abs == nullptr) ?
                     gsl_odeiv2_driver_alloc_standard_new(&sys, get_step_type(step_type_idx),
                                                          hstart, epsabs, epsrel, a_y, a_dydt) :
                     gsl_odeiv2_driver_alloc_scaled_new(&sys, get_step_type(step_type_idx), hstart,
                                                        epsabs, epsrel, a_y, a_dydt, scale_abs)) {}
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
        int set_driver(Driver& d) { gsl_odeiv2_step_set_driver(m_step, d.m_driver); }
        int apply(double t, double h, double y[], double yerr[], const double dydt_in[],
                  double dydt_out[], const System sys) {
            return gsl_odeiv2_step_apply(m_step, t, h, y, yerr, dydt_in, dydt_out, &sys);
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
        int init(Control& c, double eps_abs, double eps_rel, double a_y, double a_dydt){
            gsl_odeiv2_control_init(m_control, eps_abs, eps_rel, a_y, a_dydt);
        }
        int hadjust(Step& s, const double y[], const double yerr[], const double dydt[], double *h){
            return gsl_odeiv2_control_hadjust(m_control, s.m_step, y, yerr, dydt, h);
        }
    };
    struct Evolve {
        gsl_odeiv2_evolve *m_evolve;
        Evolve(size_t dim) : m_evolve(gsl_odeiv2_evolve_alloc(dim)) {}
        ~Evolve() { gsl_odeiv2_evolve_free(this->m_evolve); }
        int apply(Control& con, Step& step, System sys, double *t,
                  double t1, double *h, double y[]) {
            gsl_odeiv2_evolve_apply(m_evolve, con.m_control, step.m_step, &sys, t, t1, h, y);
        }
        int reset() { return gsl_odeiv2_evolve_reset(m_evolve); }
    };

    int func(double t, const double y[], double f[], void * params);
    int jac(double x, const double y[], double *dfdy, double dfdt[], void *params);

    class PyGslOdeiv2 {
    public:
        PyObject *py_f, *py_j;
        size_t ny;
        std::vector<double> xout;
        std::vector<double> yout;

        PyGslOdeiv2(PyObject * py_f, PyObject * py_j, size_t ny) :
            py_f(py_f), py_j(py_j), ny(ny) {}
        size_t integrate_adaptive(PyObject *py_y0, double x0, double xend,
                                  double dx0, double atol, double rtol,
                                  int step_type_idx){
            System system = {func, jac, this->ny, static_cast<void*>(this)};
            Step step(step_type_idx, this->ny);
            Control control(atol, rtol);
            Evolve evolve(this->ny);
            std::vector<double> tout;
            std::vector<double> yout;
            int info;
            size_t nsteps = 0;
            while (x0 < xend){
                info = evolve.apply(control, step, system, &x0, xend, &dx0,
                                    (double*)PyArray_GETPTR1(py_y0, 0));
                if (info != GSL_SUCCESS)
                    throw std::runtime_error("evolve.apply failed");
                this->xout.push_back(x0);
                for (size_t i=0; i<this->ny; ++i)
                    this->yout.push_back(*(double*)PyArray_GETPTR1(py_y0, i));
                nsteps++;
            }
            return nsteps;
        }
        int integrate_fixed_step (double t, double t1, double *y, size_t n_steps,
                                  double hstart, double epsabs,
                                  double epsrel,
                                  int step_type_idx, double * tout, double * yout,
                                  double hmin=0, double hmax=0)
        {
            System sys = {func, jac, this->ny, static_cast<void*>(this)};
            Driver d(step_type_idx, this->ny, sys, hstart, epsabs, epsrel);
            double ti, dt = (t1-t)/n_steps;
            int status;
            if (hmax > 0) d.set_hmax(hmax);
            if (hmin > 0) d.set_hmin(hmin);
            for (unsigned int i=0; i<n_steps; ++i){
                ti = t + dt;
                status = d.apply(&t, ti, y);
                if (status != GSL_SUCCESS)
                    break;
                tout[i] = t;
                for (unsigned int j=0; i<(this->ny); ++j)
                    yout[i*this->ny + j] = y[j];
            }
            return status;
        }
    };

} // namespace gslodeiv2

int gslodeiv2::func(double xval, const double y[], double dydx[], void * params)
{
    auto obj = static_cast<gslodeiv2::PyGslOdeiv2*>(params);
    npy_intp dims[1] { obj->ny } ;
    PyObject * py_yarr = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, const_cast<double *>(y));
    PyObject * py_dydx = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, const_cast<double *>(dydx));
    PyObject * py_arglist = Py_BuildValue("(dOO)", xval, py_yarr, py_dydx);
    PyObject * py_result = PyEval_CallObject(obj->py_f, py_arglist);
    Py_DECREF(py_arglist);
    Py_DECREF(py_dydx);
    Py_DECREF(py_yarr);
    if (py_result == nullptr){
        PyErr_SetString(PyExc_RuntimeError, "f() failed");
        return 1;
    } else if (py_result != Py_None){
        // py_result is not None
        PyErr_SetString(PyExc_RuntimeError, "f() did not return None");
    }
    Py_DECREF(py_result);
    return GSL_SUCCESS;
}

int gslodeiv2::jac(double xval, const double y[], double *dfdy, double dfdx[], void *params)
{
    auto obj = static_cast<gslodeiv2::PyGslOdeiv2*>(params);
    npy_intp ydims[1] { obj->ny };
    npy_intp Jdims[2] { obj->ny, obj->ny };
    PyObject * py_yarr = PyArray_SimpleNewFromData(1, ydims, NPY_DOUBLE, const_cast<double *>(y));
    PyObject * py_jmat = PyArray_SimpleNewFromData(2, Jdims, NPY_DOUBLE, const_cast<double *>(dfdy));
    PyObject * py_dfdx = PyArray_SimpleNewFromData(1, ydims, NPY_DOUBLE, const_cast<double *>(dfdx));
    PyObject * py_arglist = Py_BuildValue("(dOOO)", xval, py_yarr, py_jmat, py_dfdx);
    PyObject * py_result = PyEval_CallObject(obj->py_j, py_arglist);
    Py_DECREF(py_arglist);
    Py_DECREF(py_dfdx);
    Py_DECREF(py_jmat);
    Py_DECREF(py_yarr);
    if (py_result == nullptr){
        PyErr_SetString(PyExc_RuntimeError, "f() failed");
        return 1;
    } else if (py_result != Py_None){
        // py_result is not None
        PyErr_SetString(PyExc_RuntimeError, "f() did not return None");
    }
    Py_DECREF(py_result);
    return GSL_SUCCESS;
}

// class Block {
//     gsl_block * m_block;
// public:
//     Block(size_t dim) : m_block(gsl_block_calloc(dim)) {}
//     ~Block() { gsl_block_free(m_block); }
// }
// class Matrix {
//     gsl_matrix * m_matrix;
// public:
//     Matrix(size_t dim1, size_t dim2) : m_matrix(gsl_matrix_calloc(dim1, dim2)) {}
//     ~Matrix() { gsl_matrix_free(m_matrix); }
// }
#endif /* GSLODEIV2_H_RVEFMVJH25BGNEJDX4WHGA73YI */