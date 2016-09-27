CXX ?= g++
LIBS ?=-lgsl -lgslcblas -lm
CXXFLAGS ?= -std=c++11 -Wall -Wextra -Werror -pedantic -g -ggdb -O0
#-D_GLIBCXX_DEBUG
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../include
DEFINES ?=


.PHONY: test clean

test: test_gsl_odeiv2_cxx
	./test_gsl_odeiv2_cxx --abortx 1
	python3 ./_test_gsl_odeiv2_cxx.py

clean:
	rm catch.hpp
	rm test_gsl_odeiv2_cxx

test_%: test_%.cpp ../include/gsl_odeiv2_cxx.hpp catch.hpp testing_utils.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS)

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<
