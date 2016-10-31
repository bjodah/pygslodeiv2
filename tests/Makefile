CXX ?= g++
LIBS ?=-lgsl -lgslcblas -lm
CXXFLAGS ?= -std=c++11 -Wall -Wextra -Werror -pedantic -g -ggdb -O0
#-D_GLIBCXX_DEBUG
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../include -I../external/anyode/include
DEFINES ?=
OPENMP_FLAG ?= -fopenmp
OPENMP_LIB ?= -lgomp


.PHONY: test clean

test: test_gsl_odeiv2_anyode test_gsl_odeiv2_anyode_parallel test_gsl_odeiv2_anyode_autorestart test_gsl_odeiv2_cxx _test_gsl_odeiv2_anyode.py
	./test_gsl_odeiv2_anyode --abortx 1
	./test_gsl_odeiv2_anyode_parallel --abortx 1
	./test_gsl_odeiv2_anyode_autorestart --abortx 1
	./test_gsl_odeiv2_cxx --abortx 1
	python3 ./_test_gsl_odeiv2_anyode.py

clean:
	rm catch.hpp
	rm test_gsl_odeiv2_cxx

test_%: test_%.cpp ../include/gsl_odeiv2_cxx.hpp catch.hpp testing_utils.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS)

test_gsl_odeiv2_anyode_parallel: test_gsl_odeiv2_anyode_parallel.cpp ../include/gsl_odeiv2_*.hpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS) $(OPENMP_LIB)

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<
