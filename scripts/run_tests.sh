#!/bin/bash
# Usage
#   $ ./scripts/run_tests.sh
# or
#   $ ./scripts/run_tests.sh --cov pycvodes --cov-report html
set -ex
${PYTHON:-python3} -m pytest --pyargs ${CI_REPO_NAME:-$(basename $(realpath $(dirname $BASH_SOURCE)/../))} --doctest-modules --flakes $@
MPLBACKEND=Agg ${PYTHON:-python3} -m doctest README.rst
