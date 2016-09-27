#pragma once

#include <limits>
#include <string>
#include <unordered_map>
#include <sstream>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "anyode/anyode.hpp"


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
        if (name == "RK2")
            return StepType::RK2;
        else if (name == "RK4")
            return StepType::RK4;
        else if (name == "RKF45")
            return StepType::RKF45;
        else if (name == "RKCK")
            return StepType::RKCK;
        else if (name == "RK8PD")
            return StepType::RK8PD;
        else if (name == "RK1IMP")
            return StepType::RK1IMP;
        else if (name == "RK2IMP")
            return StepType::RK2IMP;
        else if (name == "RK4IMP")
            return StepType::RK4IMP;
        else if (name == "BSIMP")
            return StepType::BSIMP;
        else if (name == "MSADAMS")
            return StepType::MSADAMS;
        else if (name == "MSBDF")
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

    // struct System {
    //     gsl_odeiv2_system m_sys;
    //     System(RhsFn rhs, JacFn jac, std::size_t dim, void * user_data) : m_sys({rhs, jac, dim, user_data}) {}
    // };

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
                                             << "gsl_odeiv2_evolve_apply failed at t= "
                                             << curr_x << " with stepsize="
                                             << curr_dx << " (step size too small).");
                else
                    throw std::runtime_error(StreamFmt() << std::scientific
                                             << "gsl_odeiv2_evolve_apply failed at t= "
                                             << curr_x << " with stepsize="
                                             << curr_dx << " (unknown error code: "<< info << ").");
                xout.push_back(curr_x);
            }
            return std::pair<std::vector<double>, std::vector<double>>(xout, yout);
        }

        void predefined(std::size_t nt,
                        const double * const tout,
                        const double * const y0,
                        double * const yout){
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
                                             << "gsl_odeiv2_driver_apply failed at t= "
                                             << tout[idx + 1] << " with stepsize="
                                             << this->m_drv.m_driver->e->last_step
                                             << " (maximum number of iterations reached).");
                } else if (info == GSL_ENOPROG) {
                    throw std::runtime_error(StreamFmt() << std::scientific
                                             << "gsl_odeiv2_driver_apply failed at t= "
                                             << tout[idx + 1] << " with stepsize="
                                             << this->m_drv.m_driver->e->last_step
                                             << " (step size too small).");
                } else {
                    throw std::runtime_error(StreamFmt() << std::scientific
                                             << "gsl_odeiv2_driver_apply failed at t= "
                                             << tout[idx + 1] << " with stepsize="
                                             << this->m_drv.m_driver->e->last_step
                                             << " (unknown error code: "<< info <<")");
                }
            }
        }
    };

    int handle_status_(AnyODE::Status status){
        switch (status){
        case AnyODE::Status::success:
            return GSL_SUCCESS;
        case AnyODE::Status::recoverable_error:
            return 99;  // Any non-GSL specific error code
        case AnyODE::Status::unrecoverable_error:
            return GSL_EBADFUNC;
        default:
            throw std::runtime_error("impossible (this is for silencing -Wreturn-type)");
        }
    }

    template<class OdeSys>  // the *_cb functions allows us to pass C-pointer of a method.
    int rhs_cb(double t, const double y[], double ydot[], void *user_data) {
        auto& odesys = *static_cast<OdeSys*>(user_data);
        AnyODE::Status status = odesys.rhs(t, y, ydot);
        return handle_status_(status);
    }

    template <class OdeSys>
    int jac_dense_cb(double t, const double y[], double *dfdy, double dfdt[], void *user_data) {
        // callback of req. signature wrapping OdeSys method.
        auto& odesys = *static_cast<OdeSys*>(user_data);
        AnyODE::Status status = odesys.dense_jac_rmaj(t, y, nullptr, dfdy, odesys.get_ny(), dfdt);
        return handle_status_(status);
    }


    template <class OdeSys>
    GSLIntegrator get_integrator(OdeSys * odesys,
                                 const double atol,
                                 const double rtol,
                                 const StepType styp,
                                 const double dx0=0.0,
                                 const double dx_min=0.0,
                                 const double dx_max=0.0,
                                 const long int mxsteps=0
                                 )
    {
        const int ny = odesys->get_ny();
        GSLIntegrator integr {rhs_cb<OdeSys>, jac_dense_cb<OdeSys>, ny, styp, dx0, atol, rtol, static_cast<void*>(odesys)};
        if (dx0 != 0.0)
            integr.m_drv.set_init_step(dx0);
        if (dx_min != 0.0)
            integr.m_drv.set_min_step(dx_min);
        if (dx_max != 0.0)
            integr.m_drv.set_max_step(dx_max);
        if (mxsteps)
            integr.m_drv.set_max_num_steps(mxsteps);
        return integr;
    }

    void set_integration_info(std::unordered_map<std::string, int>& info,
                              const GSLIntegrator& integrator){
        info["n_steps"] = integrator.get_n_steps();
        info["n_failed_steps"] = integrator.get_n_failed_steps();
    }


    template <class OdeSys>
    std::pair<std::vector<double>, std::vector<double> >
    simple_adaptive(OdeSys * const odesys,
                    const double atol,
                    const double rtol,
                    const StepType styp,
                    const double * const y0,
                    const double x0,
                    const double xend,
                    double dx0=0.0,
                    const double dx_min=0.0,
                    const double dx_max=0.0,
                    long int mxsteps=0)
    {
        if (dx0 == 0.0){
            if (x0 == 0)
                dx0 = std::numeric_limits<double>::epsilon() * 100;
            else
                dx0 = std::numeric_limits<double>::epsilon() * 100 * x0;
        }
        if (mxsteps == 0)
            mxsteps = 500;
        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, styp, dx0, dx_min, dx_max, mxsteps);
        odesys->integrator = static_cast<void*>(&integr);
        auto result = integr.adaptive(x0, xend, y0);
        odesys->last_integration_info.clear();
        set_integration_info(odesys->last_integration_info, integr);
        return result;
    }

    template <class OdeSys>
    void simple_predefined(OdeSys * const odesys,
                           const double atol,
                           const double rtol,
                           const StepType styp,
                           const double * const y0,
                           const std::size_t nout,
                           const double * const xout,
                           double * const yout,
                           double dx0=0.0,
                           const double dx_min=0.0,
                           const double dx_max=0.0,
                           long int mxsteps=0)
    {
        if (dx0 == 0.0){
            if (xout[0] == 0)
                dx0 = std::numeric_limits<double>::epsilon() * 100;
            else
                dx0 = std::numeric_limits<double>::epsilon() * 100 * xout[0];
        }
        if (mxsteps == 0)
            mxsteps = 500;
        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, styp, dx0, dx_min, dx_max, mxsteps);
        odesys->integrator = static_cast<void*>(&integr);
        integr.predefined(nout, xout, y0, yout);
        odesys->last_integration_info.clear();
        set_integration_info(odesys->last_integration_info, integr);
    }

}
