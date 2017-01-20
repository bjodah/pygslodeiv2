#pragma once

#include <chrono>
#include <functional>
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
        if (odesys.record_rhs_xvals)
            odesys.last_integration_info_vecdbl["rhs_xvals"].push_back(t);
        AnyODE::Status status = odesys.rhs(t, y, ydot);
        return handle_status_(status);
    }

    template <class OdeSys>
    int jac_dense_cb(double t, const double y[], double *dfdy, double dfdt[], void *user_data) {
        // callback of req. signature wrapping OdeSys method.
        auto& odesys = *static_cast<OdeSys*>(user_data);
        if (odesys.record_jac_xvals)
            odesys.last_integration_info_vecdbl["jac_xvals"].push_back(t);
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
                                 const long int mxsteps=0,
                                 const bool record_order=false,
                                 const bool record_fpe=false
                                 )
    {
        const int ny = odesys->get_ny();
        GSLIntegrator integr {rhs_cb<OdeSys>, jac_dense_cb<OdeSys>, ny, styp, dx0, atol, rtol, static_cast<void*>(odesys)};
        integr.record_order = record_order;
        integr.record_fpe = record_fpe;
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
                    const double dx_max=0.0,
                    int autorestart=0,
                    bool return_on_error=false)
    {
        if (dx0 == 0.0)
            dx0 = odesys->get_dx0(x0, y0);
        if (dx0 == 0.0){
            if (x0 == 0)
                dx0 = std::numeric_limits<double>::epsilon() * 100;
            else
                dx0 = std::numeric_limits<double>::epsilon() * 100 * x0;
        }
        if (mxsteps == 0)
            mxsteps = 500;
        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, styp, dx0, dx_min, dx_max, mxsteps,
                                             odesys->record_order, odesys->record_fpe);
        odesys->integrator = static_cast<void*>(&integr);

        odesys->last_integration_info.clear();
        odesys->last_integration_info_dbl.clear();
        odesys->last_integration_info_vecdbl.clear();
        odesys->last_integration_info_vecint.clear();
        if (odesys->record_rhs_xvals)
            odesys->last_integration_info_vecdbl["rhs_xvals"] = {};
        if (odesys->record_jac_xvals)
            odesys->last_integration_info_vecdbl["jac_xvals"] = {};

        std::time_t cput0 = std::clock();
        auto t_start = std::chrono::high_resolution_clock::now();

        auto result = integr.adaptive(x0, xend, y0, autorestart, return_on_error,
            ((odesys->use_get_dx_max) ? static_cast<gsl_odeiv2_cxx::get_dx_max_fn>(
                std::bind(&OdeSys::get_dx_max, odesys, std::placeholders::_1 , std::placeholders::_2))
	     : gsl_odeiv2_cxx::get_dx_max_fn()));

        odesys->last_integration_info_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
        odesys->last_integration_info_dbl["time_wall"] = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();

        if (odesys->record_order)
            odesys->last_integration_info_vecint["orders"] = integr.orders_seen;
        if (odesys->record_fpe)
            odesys->last_integration_info_vecint["fpes"] = integr.fpes_seen;

        set_integration_info<OdeSys>(odesys, integr);
        return result;
    }

    template <class OdeSys>
    int simple_predefined(OdeSys * const odesys,
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
                           const double dx_max=0.0,
                           int autorestart=0,
                           bool return_on_error=false
                           )
    {
        if (dx0 == 0.0)
            dx0 = odesys->get_dx0(xout[0], y0);
        if (dx0 == 0.0){
            if (xout[0] == 0)
                dx0 = std::numeric_limits<double>::epsilon() * 1000;
            else
                dx0 = std::numeric_limits<double>::epsilon() * 1000 * xout[0];
        }
        if (mxsteps == 0)
            mxsteps = 500;
        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, styp, dx0, dx_min, dx_max, mxsteps,
                                             odesys->record_order, odesys->record_fpe);
        odesys->integrator = static_cast<void*>(&integr);

        odesys->last_integration_info.clear();
        odesys->last_integration_info_dbl.clear();
        odesys->last_integration_info_vecdbl.clear();
        if (odesys->record_rhs_xvals)
            odesys->last_integration_info_vecdbl["rhs_xvals"] = {};
        if (odesys->record_jac_xvals)
            odesys->last_integration_info_vecdbl["jac_xvals"] = {};

        std::time_t cput0 = std::clock();
        auto t_start = std::chrono::high_resolution_clock::now();

        int nreached = integr.predefined(nout, xout, y0, yout, autorestart, return_on_error,
            ((odesys->use_get_dx_max) ? static_cast<gsl_odeiv2_cxx::get_dx_max_fn>(
                std::bind(&OdeSys::get_dx_max, odesys, std::placeholders::_1 , std::placeholders::_2))
	     : gsl_odeiv2_cxx::get_dx_max_fn()));

        odesys->last_integration_info_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
        odesys->last_integration_info_dbl["time_wall"] = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();

        if (odesys->record_order)
            odesys->last_integration_info_vecint["orders"] = integr.orders_seen;
        if (odesys->record_fpe)
            odesys->last_integration_info_vecint["fpes"] = integr.fpes_seen;

        set_integration_info<OdeSys>(odesys, integr);
        return nreached;
    }

}
