CXX ?= g++
LIBS ?=-lgsl -lgslcblas -lm
CXXFLAGS ?= -std=c++11 -Wall -Wextra -Werror -pedantic -g -ggdb -O0
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../pygslodeiv2/include -I../external/anyode/include
EXTRA_FLAGS ?=
CXXFLAGS += $(EXTRA_FLAGS)
OPENMP_FLAG ?= -fopenmp
OPENMP_LIB ?= -lgomp


.PHONY: test clean

test: test_gsl_odeiv2_anyode test_gsl_odeiv2_anyode_parallel test_gsl_odeiv2_anyode_autorestart test_gsl_odeiv2_cxx _test_gsl_odeiv2_anyode.py
	./test_gsl_odeiv2_anyode --abortx 1
ifeq ($(CXX),clang++-6.0)
	@echo "skipping openmp when using clang"
else
	./test_gsl_odeiv2_anyode_parallel --abortx 1
endif
	./test_gsl_odeiv2_anyode_autorestart --abortx 1
	./test_gsl_odeiv2_cxx --abortx 1
	python3 ./_test_gsl_odeiv2_anyode.py

clean:
	rm catch.hpp
	rm test_gsl_odeiv2_anyode
	rm test_gsl_odeiv2_anyode_parallel
	rm test_gsl_odeiv2_anyode_autorestart
	rm test_gsl_odeiv2_cxx


test_%: test_%.cpp ../pygslodeiv2/include/gsl_odeiv2_cxx.hpp catch.hpp testing_utils.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $< $(LIBS)

test_gsl_odeiv2_anyode_parallel: test_gsl_odeiv2_anyode_parallel.cpp ../pygslodeiv2/include/gsl_odeiv2_*.hpp
ifeq ($(CXX),clang++-6.0)
	@echo "skipping openmp when using clang"
else
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) $(INCLUDE) -o $@ $< $(LIBS) $(OPENMP_LIB)
endif

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<
