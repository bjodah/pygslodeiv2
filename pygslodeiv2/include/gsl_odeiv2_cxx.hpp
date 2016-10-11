#pragma once

#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

namespace {
    class StreamFmt
    {
        std::stringstream m_s;
    public:
        StreamFmt() {}
        ~StreamFmt() {}

        template <typename T>
        StreamFmt& operator << (const T& v) {
            this->m_s << v;
            return *this;
        }

        std::string str() const {
            return this->m_s.str();
        }
        operator std::string() const {
            return this->m_s.str();
        }

    };
}


namespace gsl_odeiv2_cxx {

    enum class StepType : int {RK2=0, RK4=1, RKF45=2, RKCK=3, RK8PD=4,
            RK1IMP=5, RK2IMP=6, RK4IMP=7, BSIMP=8, MSADAMS=9, MSBDF=10};

    StepType styp_from_name(std::string name){
        if (name == "rk2")
            return StepType::RK2;
        else if (name == "rk4")
            return StepType::RK4;
        else if (name == "rkf45")
            return StepType::RKF45;
        else if (name == "rkck")
            return StepType::RKCK;
        else if (name == "rk8pd")
            return StepType::RK8PD;
        else if (name == "rk1imp")
            return StepType::RK1IMP;
        else if (name == "rk2imp")
            return StepType::RK2IMP;
        else if (name == "rk4imp")
            return StepType::RK4IMP;
        else if (name == "bsimp")
            return StepType::BSIMP;
        else if (name == "msadams")
            return StepType::MSADAMS;
        else if (name == "msbdf")
            return StepType::MSBDF;
        else
            throw std::runtime_error(StreamFmt() << "Unknown stepper type name: " << name);
    }
    bool requires_jacobian(StepType styp){
        if ((styp == StepType::RK1IMP) or (styp == StepType::RK2IMP) or
            (styp == StepType::RK4IMP) or (styp == StepType::BSIMP) or
            (styp == StepType::MSBDF))
            return true;
        else
            return false;
    }


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

    typedef int (*RhsFn)(double t, const double y[], double dydt[], void *params);
    typedef int (*JacFn)(double t, const double y[], double *dfdy, double dfdt[], void *params);

    struct Driver {
        gsl_odeiv2_driver *m_driver;
        Driver(const gsl_odeiv2_system * sys,
               const gsl_odeiv2_cxx::StepType styp,
               const double init_step,
               const double atol,
               const double rtol) : m_driver(gsl_odeiv2_driver_alloc_y_new(
                   sys, get_step_type(static_cast<int>(styp)), init_step, atol, rtol))
               {}
        ~Driver() { gsl_odeiv2_driver_free(this->m_driver); }
        int set_init_step(const double hstart) { return gsl_odeiv2_driver_reset_hstart(m_driver, hstart); }
        int set_min_step(double hmin) { return gsl_odeiv2_driver_set_hmin(m_driver, hmin); }
        int set_max_step(double hmax) { return gsl_odeiv2_driver_set_hmax(m_driver, hmax); }
        int set_max_num_steps(const unsigned long int nmax) { return gsl_odeiv2_driver_set_nmax(m_driver, nmax); }
        int apply(double * t, const double t1, double y[]){
            return gsl_odeiv2_driver_apply(m_driver, t, t1, y);
        }
        int reset() { return gsl_odeiv2_driver_reset(m_driver); }
    };
    struct Step {
        gsl_odeiv2_step *m_step;
        Step(const StepType styp, const size_t dim) :
            m_step(gsl_odeiv2_step_alloc(get_step_type(static_cast<int>(styp)), dim)) {}
        ~Step() { gsl_odeiv2_step_free(m_step); }
        int reset() { return gsl_odeiv2_step_reset(m_step); }
        unsigned int order() { return gsl_odeiv2_step_order(m_step); }
        int set_driver(Driver& d) { return gsl_odeiv2_step_set_driver(m_step, d.m_driver); }
        int apply(double t, double h, double y[], double yerr[], const double dydt_in[],
                  double dydt_out[], const gsl_odeiv2_system * sys) {
            return gsl_odeiv2_step_apply(m_step, t, h, y, yerr, dydt_in, dydt_out, sys);
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
        int init(double eps_abs, double eps_rel, double a_y, double a_dydt){
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
        int apply(Control& con, Step& step, gsl_odeiv2_system * sys, double *t,
                  double t1, double *h, double y[]) {
            return gsl_odeiv2_evolve_apply(m_evolve, con.m_control, step.m_step,
                                           sys, t, t1, h, y);
        }
        int reset() { return gsl_odeiv2_evolve_reset(m_evolve); }
    };

    void stderr_handle_error(const char * reason,
                             const char * file,
                             int line,
                             int gsl_errno) {
        std::fprintf(stderr, "GSL Error: %s (in %s, line %d, gsl_errno: %d)", reason, file, line, gsl_errno);
    }

    struct ErrorHandler {
        gsl_error_handler_t *m_ori_handler;
        ErrorHandler() {
            m_ori_handler = gsl_set_error_handler(&stderr_handle_error);
        }
        ~ErrorHandler() {
            if (m_ori_handler)
                this->release();
        }
        void release() {
            gsl_set_error_handler(m_ori_handler);
        }
    };

    struct GSLIntegrator{
        gsl_odeiv2_system m_sys;
        Driver m_drv;
        Step m_stp;
        Control m_ctrl;
        Evolve m_evo;
        int ny;
        GSLIntegrator(RhsFn rhs_cb, JacFn jac_cb, int ny, const StepType styp,
                      double dx0, double atol, double rtol, void * user_data) :
            m_sys(gsl_odeiv2_system({rhs_cb, jac_cb, static_cast<std::size_t>(ny), user_data})),
            m_drv(Driver(&(m_sys), styp, dx0, atol, rtol)),
            m_stp(Step(styp, ny)),
            m_ctrl(Control(atol, rtol)),
            m_evo(Evolve(ny)),
            ny(ny)
        {
            this->m_stp.set_driver(this->m_drv);
            this->m_ctrl.set_driver(this->m_drv);
            this->m_evo.set_driver(this->m_drv);
        }

        int get_n_steps() const {
            return this->m_evo.m_evolve->count;
        }
        int get_n_failed_steps() const {
            return this->m_evo.m_evolve->failed_steps;
        }

        std::pair<std::vector<double>, std::vector<double> >
        adaptive(const double x0,
                 const double xend,
                 const double * const y0){
            ErrorHandler errh;  // has side-effects;
            std::vector<double> xout;
            std::vector<double> yout;
            double curr_x = x0;
            double curr_dx = this->m_drv.m_driver->h;
            xout.push_back(curr_x);
            yout.insert(yout.end(), y0, y0 + ny);
            while (curr_x < xend){
                if (curr_dx > this->m_drv.m_driver->hmax){
                    curr_dx = this->m_drv.m_driver->hmax;
                } else if (curr_dx < this->m_drv.m_driver->hmin ) {
                    curr_dx = this->m_drv.m_driver->hmin;
                }
                yout.insert(yout.end(), yout.end() - ny, yout.end());
                int info = this->m_evo.apply(this->m_ctrl, this->m_stp, &(this->m_sys),
                                             &curr_x, xend, &curr_dx, &(*(yout.end() - ny)));
                if (info == GSL_SUCCESS)
                    ;
                else if (info == GSL_FAILURE)
                    throw std::runtime_error(StreamFmt() << std::scientific
                        << "gsl_odeiv2_evolve_apply failed at t= " << curr_x
                        << " with stepsize=" << curr_dx << " (step size too small).");
                else
                    throw std::runtime_error(StreamFmt() << std::scientific
                        << "gsl_odeiv2_evolve_apply failed at t= " << curr_x
                        << " with stepsize=" << curr_dx
                        << " (unknown error code: "<< info << ").");
                xout.push_back(curr_x);
            }
            return std::pair<std::vector<double>, std::vector<double>>(xout, yout);
        }

        void predefined(std::size_t nt,
                        const double * const tout,
                        const double * const y0,
                        double * const yout){
            ErrorHandler errh;  // has side-effects;
            double curr_t = tout[0];
            this->m_drv.reset();
            std::copy(y0, y0 + (this->ny), yout);

            for (std::size_t idx=1; idx < nt; ++idx){
                std::copy(yout + (this->ny)*(idx-1),
                          yout + (this->ny)*idx,
                          yout + (this->ny)*idx);
                int info = this->m_drv.apply(&curr_t, tout[idx], yout + idx*(this->ny));
                if (info == GSL_SUCCESS){
                    if (tout[idx] != curr_t)
                        throw std::runtime_error("Did not reach requested time.");
                } else if (info == GSL_EMAXITER){
                    throw std::runtime_error(StreamFmt() << std::scientific
                        << "gsl_odeiv2_driver_apply failed at t= " << tout[idx + 1]
                        << " with stepsize=" << this->m_drv.m_driver->e->last_step
                        << " (maximum number of iterations reached).");
                } else if (info == GSL_ENOPROG) {
                    throw std::runtime_error(StreamFmt() << std::scientific
                        << "gsl_odeiv2_driver_apply failed at t= " << tout[idx + 1]
                        << " with stepsize=" << this->m_drv.m_driver->e->last_step
                        << " (step size too small).");
                } else {
                    throw std::runtime_error(StreamFmt() << std::scientific
                        << "gsl_odeiv2_driver_apply failed at t= " << tout[idx + 1]
                        << " with stepsize=" << this->m_drv.m_driver->e->last_step
                        << " (unknown error code: "<< info <<")");
                }
            }
        }
    };

}
