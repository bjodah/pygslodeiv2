CXX ?= g++
EXTRA_LIBS ?=-lgsl -lgslcblas -lm
WARNINGS ?= \
-Wall \
-Wextra \
-Wredundant-decls \
-Wcast-align \
-Wmissing-include-dirs \
-Wswitch-enum \
-Wswitch-default \
-Winvalid-pch \
-Wredundant-decls \
-Wformat=2 \
-Wmissing-format-attribute \
-Wformat-nonliteral \
-Wodr
CXXFLAGS ?= -std=c++11 $(WARNINGS) -Werror -pedantic -g -ggdb -O0
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../pygslodeiv2/include -I../external/anyode/include
EXTRA_FLAGS ?=
CXXFLAGS += $(EXTRA_FLAGS)
OPENMP_FLAG ?= -fopenmp
OPENMP_LIB ?= -lgomp


.PHONY: test clean

test: test_gsl_odeiv2_anyode test_gsl_odeiv2_anyode_parallel test_gsl_odeiv2_anyode_autorestart test_gsl_odeiv2_cxx _test_gsl_odeiv2_anyode.py
	./test_gsl_odeiv2_anyode --abortx 1
	./test_gsl_odeiv2_anyode_parallel --abortx 1
	./test_gsl_odeiv2_anyode_autorestart --abortx 1
	./test_gsl_odeiv2_cxx --abortx 1
	env DISTUTILS_DEBUG=1 CC=$(CXX) CFLAGS="$(EXTRA_FLAGS)" LDFLAGS="$(LDFLAGS)" LD_PRELOAD="$(PY_LD_PRELOAD)" ASAN_OPTIONS=detect_leaks=0 python3 ./_test_gsl_odeiv2_anyode.py

clean:
	rm -f doctest.h
	rm -f test_gsl_odeiv2_anyode
	rm -f test_gsl_odeiv2_anyode_parallel
	rm -f test_gsl_odeiv2_anyode_autorestart
	rm -f test_gsl_odeiv2_cxx


test_%: test_%.cpp ../pygslodeiv2/include/gsl_odeiv2_cxx.hpp doctest.h testing_utils.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -o $@ $< $(LDLIBS) $(EXTRA_LIBS)

test_gsl_odeiv2_anyode_parallel: test_gsl_odeiv2_anyode_parallel.cpp ../pygslodeiv2/include/gsl_odeiv2_*.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OPENMP_FLAG) $(INCLUDE) -o $@ $< $(LDLIBS) $(EXTRA_LIBS) $(OPENMP_LIB)

doctest.h: doctest.h.bz2
	bunzip2 -k -f $<
