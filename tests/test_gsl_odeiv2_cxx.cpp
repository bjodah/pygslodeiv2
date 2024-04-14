#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "gsl_odeiv2_cxx.hpp"
#include "testing_utils.hpp"


TEST_CASE( "methods" ) {
    const int ny = 1;
    void * user_data = nullptr;
    double dx0 = 1e-12, atol=1e-8, rtol=1e-8;
    auto intgr = gsl_odeiv2_cxx::GSLIntegrator(
	    rhs_cb, nullptr, ny, gsl_odeiv2_cxx::StepType::MSADAMS,
	    dx0, atol, rtol, user_data);
    std::vector<double> y0(1, 1.0);
    std::vector<double> tout {{0.0, 1.0}};
    std::vector<double> yout(2);
    intgr.predefined(tout.size(), &tout[0], &y0[0], &yout[0]);
    double yref = std::exp(-tout[1]);
    REQUIRE( std::abs(yout[1] - yref) < 1e-7 );
    // REQUIRE( intgr.get_n_steps() > 3 );  // m_drv, m_evo don't share n_steps for some reason
    // REQUIRE( intgr.get_n_steps() < 100 );  // m_drv, m_evo don't share n_steps for some reason
}

double get_dx_max(double /* x */, const double * const /* y */){
    return 1e-3;
}

TEST_CASE( "adaptive" ) {
    const int ny = 1;
    bool return_on_error = false;
    int autorestart = 0;
    void * user_data = nullptr;
    double dx0 = 1e-12, atol=1e-9, rtol=1e-9;
    auto intgr = gsl_odeiv2_cxx::GSLIntegrator(
	    rhs_cb, nullptr, ny, gsl_odeiv2_cxx::StepType::MSADAMS,
	    dx0, atol, rtol, user_data);
    intgr.m_drv.set_max_num_steps(1019);
    std::vector<double> y0(1, 1.0);
    double xend = 1.0;
    auto xout_yout = intgr.adaptive(0.0, xend, &y0[0], autorestart,
				    return_on_error, get_dx_max);
    auto& xout = xout_yout.first;
    auto& yout = xout_yout.second;
    REQUIRE( xout[0] == 0.0 );
    REQUIRE( yout[0] == 1.0 );
    int nt = xout.size();
    for (int idx=1; idx<nt; ++idx){
        const double yref = std::exp(-xout[idx]);
        REQUIRE( xout[idx] > 0 );
        REQUIRE( xout[idx] <= 1 );
        REQUIRE( std::abs(yout[idx] - yref) < 1e-8 );
    }
    REQUIRE( intgr.get_n_steps() > 1000 );
}

double get_dx_max2(double /* x */, const double * const /* y */){
    return 1e-3;
}

TEST_CASE( "predefined" ) {
    const int ny = 1;
    bool return_on_error = false;
    int autorestart = 0;
    void * user_data = nullptr;
    double dx0 = 1e-12, atol=1e-9, rtol=1e-9;
    auto intgr = gsl_odeiv2_cxx::GSLIntegrator(
	    rhs_cb, nullptr, ny, gsl_odeiv2_cxx::StepType::MSADAMS,
	    dx0, atol, rtol, user_data);
    intgr.m_drv.set_max_num_steps(1019);
    int nt = 702;
    double tend = 2.0;
    std::vector<double> tout(702);
    for (int idx=0; idx<nt; ++idx)
        tout[idx] = idx*tend/(nt - 1);
    std::vector<double> yout(nt*ny);
    std::vector<double> y0(1, 1.0);
    auto nout = intgr.predefined(nt, &tout[0], &y0[0], &yout[0], autorestart,
				      return_on_error, get_dx_max2);
    for (int idx=0; idx<nt; ++idx){
        const double yref = std::exp(-tout[idx]);
        REQUIRE( std::abs(yout[idx] - yref) < 1e-8 );
    }
    REQUIRE( nout == nt );
    // REQUIRE( intgr.get_n_steps() > 1000 );  // m_drv, m_evo don't share n_steps for some reason
}
