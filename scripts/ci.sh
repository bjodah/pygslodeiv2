#!/bin/bash -xe
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

cd tests/; make CXX=clang++-10 EXTRA_FLAGS=-fsanitize=address; make clean; cd -
cd tests/; make CXX=clang++-10 EXTRA_FLAGS=-fsanitize=undefined; make clean; cd -

export CC=gcc-10
export CXX=g++-10

cd tests/; make EXTRA_FLAGS=-D_GLIBCXX_DEBUG; make clean; cd -
cd tests/; make EXTRA_FLAGS=-DNDEBUG; make clean; cd -

PKG_VERSION=$(python3 setup.py --version)
python3 setup.py sdist
(cd dist/; python3 -m pip install $PKG_NAME-$PKG_VERSION.tar.gz)
(cd /; python3 -m pytest --pyargs $PKG_NAME)
python3 -m pip install --user -e .[all]
PYTHONPATH=$(pwd) PYTHON=python3 ./scripts/run_tests.sh --cov $PKG_NAME --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
