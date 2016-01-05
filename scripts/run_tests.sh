#!/bin/bash -e
# Usage
#   $ ./scripts/run_tests.sh
# or
#   $ ./scripts/run_tests.sh --cov pycvodes --cov-report html
python2.7 setup.py build_ext -i
python2.7 -m pytest
python3 setup.py build_ext -i
python3 -m pytest --doctest-modules --pep8 --flakes $@
python3 -m doctest README.rst
