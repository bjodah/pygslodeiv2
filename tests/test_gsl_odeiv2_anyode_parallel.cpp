#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "gsl_odeiv2_anyode_parallel.hpp"
#include "testing_utils.hpp"

TEST_CASE( "decay_adaptive" ) {
    std::vector<double> k {{ 2.0, 3.0}};
    Decay odesys1(k[0]);
    Decay odesys2(k[1]);
    std::vector<Decay *> systems {{ &odesys1, &odesys2 }};
    std::vector<double> y0 {{ 5.0, 7.0 }};
    std::vector<double> t0 {{ 1.0, 3.0 }};  // delta = 2
    std::vector<double> tend {{ 2.0, 5.0 }};  // delta = 3
    int mxsteps = 0;  // => default
    double atol = 1e-10;
    std::vector<double> dx0 {{ 1e-10, 1e-10 }};
    std::vector<double> dx_min {{ 1e-16, 1e-16 }};
    std::vector<double> dx_max {{ 1.0, 1.0 }};

    auto result = gsl_odeiv2_anyode_parallel::multi_adaptive(
            systems, atol, 1e-10, gsl_odeiv2_cxx::StepType::RKCK, &y0[0], &t0[0], &tend[0],
            mxsteps, &dx0[0], &dx_min[0], &dx_max[0]);
    for (int idx=0; idx<2; ++idx){
        const auto& tout = result[idx].first;
        const auto& yout = result[idx].second;
        for (unsigned j=0; j<tout.size(); ++j){
            REQUIRE( std::abs(y0[idx]*std::exp(tout[0]-tout[j]) - yout[j]) < 1e-8 );
        }
        REQUIRE( systems[idx]->current_info.nfo_int["n_steps"] > 1 );
        REQUIRE( systems[idx]->current_info.nfo_int["n_steps"] < 997 );
    }
}
