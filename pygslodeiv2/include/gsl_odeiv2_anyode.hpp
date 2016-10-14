#pragma once

#include <chrono>
#include "anyode/anyode.hpp"
#include "gsl_odeiv2_cxx.hpp"


namespace gsl_odeiv2_anyode {

    using gsl_odeiv2_cxx::StepType;
    using gsl_odeiv2_cxx::GSLIntegrator;


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

    template <class OdeSys>
    void set_integration_info(OdeSys * odesys, const GSLIntegrator& integrator){
        odesys->last_integration_info["nfev"] = odesys->nfev;
        odesys->last_integration_info["njev"] = odesys->njev;
        odesys->last_integration_info["n_steps"] = integrator.get_n_steps();
        odesys->last_integration_info["n_failed_steps"] = integrator.get_n_failed_steps();
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
                    long int mxsteps=0,
                    double dx0=0.0,
                    const double dx_min=0.0,
                    const double dx_max=0.0)
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
        std::time_t cput0 = std::clock();
        auto t_start = std::chrono::high_resolution_clock::now();

        auto result = integr.adaptive(x0, xend, y0);

        odesys->last_integration_info_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
        odesys->last_integration_info_dbl["time_wall"] = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();
        odesys->last_integration_info.clear();
        set_integration_info<OdeSys>(odesys, integr);
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
                           long int mxsteps=0,
                           double dx0=0.0,
                           const double dx_min=0.0,
                           const double dx_max=0.0
                           )
    {
        if (dx0 == 0.0){
            if (xout[0] == 0)
                dx0 = std::numeric_limits<double>::epsilon() * 1000;
            else
                dx0 = std::numeric_limits<double>::epsilon() * 1000 * xout[0];
        }
        if (mxsteps == 0)
            mxsteps = 500;
        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, styp, mxsteps, dx0, dx_min, dx_max);
        odesys->integrator = static_cast<void*>(&integr);
        std::time_t cput0 = std::clock();
        auto t_start = std::chrono::high_resolution_clock::now();

        integr.predefined(nout, xout, y0, yout);

        odesys->last_integration_info_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
        odesys->last_integration_info_dbl["time_wall"] = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();
        odesys->last_integration_info.clear();
        set_integration_info<OdeSys>(odesys, integr);
    }

}
