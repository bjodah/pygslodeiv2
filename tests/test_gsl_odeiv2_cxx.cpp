// C++11 source code.
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "gsl_odeiv2_cxx.hpp"
#include "testing_utils.hpp"


TEST_CASE( "methods", "[GSLIntegrator]" ) {
    const int ny = 1;
    void * user_data = nullptr;
    double dx0 = 1e-12, atol=1e-8, rtol=1e-8;
    auto intgr = gsl_odeiv2_cxx::GSLIntegrator(rhs_cb, nullptr, ny, gsl_odeiv2_cxx::StepType::MSADAMS,
                                               dx0, atol, rtol, user_data);
    std::vector<double> y0 {{1.0}};
    std::vector<double> tout {{0.0, 1.0}};
    std::vector<double> yout(2);
    intgr.predefined(tout.size(), &tout[0], &y0[0], &yout[0]);
    double yref = std::exp(-tout[1]);
    REQUIRE( std::abs(yout[1] - yref) < 1e-7 );
}
