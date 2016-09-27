#pragma once

#include <Python.h>
#include <numpy/arrayobject.h>

#include <chrono> // std::clock
#include <vector> // std::vector

#include "gsl_odeiv2_cxx.hpp"

namespace gsl_odeiv2_numpy{

    using gsl_odeiv2_cxx::StepType;

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


    struct PyGslOdeiv2 : public AnyODE::OdeSysBase {
        PyObject *py_rhs, *py_jac;
        size_t ny;
        size_t nfev, njev;
        double time_cpu, time_wall;

        PyGslOdeiv2(PyObject * py_rhs, PyObject * py_jac, size_t ny) :
            py_rhs(py_rhs), py_jac(py_jac), ny(ny) {}

        virtual int get_ny() const override { return ny; }

        std::pair<std::vector<double>, std::vector<double> >
        adaptive(PyObject *py_y0, double x0, double xend,
                 double atol, double rtol, StepType styp,
                 double dx0, double dx_min=.0, double dx_max=.0, int mxsteps=0){
            std::clock_t cputime0 = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();
            nfev = 0; njev = 0;
            auto handler = PyErrorHandler();
            auto result = simple_adaptive(this, atol, rtol, styp, (double*)PyArray_GETPTR1(py_y0, 0),
                                          x0, xend, dx0, dx_min, dx_max, mxsteps);
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
            this->time_wall = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();
            return result;
        }

        void predefined(PyObject *py_xout, PyObject *py_yout,
                        double atol, double rtol, StepType styp,
                        double dx0, double dx_min=0.0, double dx_max=0.0, int mxsteps=0)
        {
            std::clock_t cputime0 = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();
            auto handler = PyErrorHandler();
            npy_intp nt = PyArray_DIMS(py_xout)[0];
            double * xout = (double*)PyArray_GETPTR1(py_xout, 0);
            double * y0 = (double *)PyArray_GETPTR2(py_yout, 0, 0);
            double * yout = (double *)PyArray_GETPTR2(py_yout, 0, 0);
            nfev = 0; njev = 0;
            simple_predefined(this, atol, rtol, styp, y0, nt, xout, yout, dx0, dx_min, dx_max, mxsteps);
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
            this->time_wall = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();
        }

        AnyODE::Status handle_py_status_(PyObject * py_result, const std::string what_arg){
            if (py_result == nullptr){
                throw std::runtime_error(what_arg + " failed");
            } else if (py_result == Py_None){
                Py_DECREF(py_result);
                return AnyODE::Status::success;
            }
            long result = PyInt_AsLong(py_result);
            Py_DECREF(py_result);
            if ((PyErr_Occurred() and result == -1) or result == static_cast<long int>(AnyODE::Status::unrecoverable_error))
                return AnyODE::Status::unrecoverable_error;
            else if (result == static_cast<long int>(AnyODE::Status::recoverable_error))
                return AnyODE::Status::recoverable_error;
            else if (result == static_cast<long int>(AnyODE::Status::success))
                return AnyODE::Status::success;
            throw std::runtime_error(what_arg + " did not return None, -1, 0 or 1");
        }

        virtual AnyODE::Status rhs(double xval,
                                   const double * const y,
                                   double * const dydx) override
        {
            npy_intp dims[1] { static_cast<npy_intp>(this->ny) } ;
            npy_intp y_strides[1] { sizeof(double) };
            PyObject* py_yarr = PyArray_New(
                                            &PyArray_Type, 1, dims, NPY_DOUBLE, y_strides,
                                            static_cast<void *>(const_cast<double *>(y)), sizeof(double),
                                            NPY_ARRAY_C_CONTIGUOUS, nullptr);
            PyObject * py_dydx = PyArray_SimpleNewFromData(
                                                           1, dims, NPY_DOUBLE, static_cast<void*>(dydx));
            PyObject * py_arglist = Py_BuildValue("(dOO)", xval, py_yarr, py_dydx);
            PyObject * py_result = PyEval_CallObject(this->py_rhs, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_dydx);
            Py_DECREF(py_yarr);
            this->nfev++;
            return handle_py_status_(py_result, "rhs");
        }

        virtual AnyODE::Status dense_jac_rmaj(double xval,
                                              const double * const y,
                                              const double * const fy, // always nullptr when using GSL
                                              double * const dfdy,
                                              long int ldim,
                                              double * const dfdx) override
        {
            npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
            npy_intp y_strides[1] { sizeof(double) };
            npy_intp Jdims[2] { static_cast<npy_intp>(this->ny), static_cast<npy_intp>(this->ny) };
            PyObject* py_yarr = PyArray_New(
                                            &PyArray_Type, 1, ydims, NPY_DOUBLE, y_strides,
                                            static_cast<void *>(const_cast<double *>(y)), sizeof(double),
                                            NPY_ARRAY_C_CONTIGUOUS, nullptr);
            PyObject * py_jmat = PyArray_SimpleNewFromData(2, Jdims, NPY_DOUBLE, const_cast<double *>(dfdy));
            PyObject * py_dfdx = PyArray_SimpleNewFromData(1, ydims, NPY_DOUBLE, const_cast<double *>(dfdx));
            PyObject * py_arglist = Py_BuildValue("(dOOO)", xval, py_yarr, py_jmat, py_dfdx);
            PyObject * py_result = PyEval_CallObject(this->py_jac, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_dfdx);
            Py_DECREF(py_jmat);
            Py_DECREF(py_yarr);
            this->njev++;
            return handle_py_status_(py_result, "rhs");
        }

    };

}
