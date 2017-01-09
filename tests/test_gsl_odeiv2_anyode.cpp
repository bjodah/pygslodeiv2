// C++11 source code.
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "gsl_odeiv2_anyode.hpp"
#include "testing_utils.hpp"


TEST_CASE( "decay_adaptive", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    double y0 = 1.0;
    double dx0 = 1e-9;
    long int mxsteps = 500;
    auto tout_yout = gsl_odeiv2_anyode::simple_adaptive(&odesys, 1e-10, 1e-10,
                                                        gsl_odeiv2_cxx::StepType::MSADAMS,
                                                        &y0, 0.0, 1.0, mxsteps, dx0);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    REQUIRE( tout.size() == yout.size() );
    for (uint i = 0; i < tout.size(); ++i){
        REQUIRE( std::abs(std::exp(-tout[i]) - yout[i]) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["nfev"] > 1 );
    REQUIRE( odesys.last_integration_info["nfev"] < 997 );
    REQUIRE( odesys.last_integration_info["n_steps"] > 1 );
    REQUIRE( odesys.last_integration_info["n_steps"] < 997 );
}


TEST_CASE( "decay_adaptive_get_dx_max", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    double y0 = 1.0;
    double dx0 = 1e-9;
    odesys.use_get_dx_max = true;
    auto tout_yout = gsl_odeiv2_anyode::simple_adaptive(&odesys, 1e-10, 1e-10,
                                                        gsl_odeiv2_cxx::StepType::MSADAMS,
                                                        &y0, 0.0, 1.0, 2100, dx0, 0.0, 1e-3);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    REQUIRE( tout.size() == yout.size() );
    for (uint i = 0; i < tout.size(); ++i){
        REQUIRE( std::abs(std::exp(-tout[i]) - yout[i]) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 2000 );
}


TEST_CASE( "decay_adaptive_dx_max", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    double y0 = 1.0;
    double dx0 = 1e-9;
    auto tout_yout = gsl_odeiv2_anyode::simple_adaptive(&odesys, 1e-10, 1e-10,
                                                        gsl_odeiv2_cxx::StepType::MSADAMS,
                                                        &y0, 0.0, 1.0, 1100, dx0, 0.0, 1e-3);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    REQUIRE( tout.size() == yout.size() );
    for (uint i = 0; i < tout.size(); ++i){
        REQUIRE( std::abs(std::exp(-tout[i]) - yout[i]) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 998 );
}
