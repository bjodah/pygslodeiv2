#!/bin/bash -xe
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

cd tests/; make EXTRA_FLAGS=-D_GLIBCXX_DEBUG; make clean; cd -
cd tests/; make EXTRA_FLAGS=-DNDEBUG; make clean; cd -
cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=address; make clean; cd -
cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=undefined; make clean; cd -

python2.7 setup.py sdist
for PYTHON in python2.7 python3; do
    (cd dist/; $PYTHON -m pip install $PKG_NAME-$($PYTHON ../setup.py --version).tar.gz)
    (cd /; $PYTHON -m pytest --pyargs $PKG_NAME)
    $PYTHON -m pip install --user -e .[all]
done
PYTHONPATH=$(pwd) PYTHON=python2.7 ./scripts/run_tests.sh
PYTHONPATH=$(pwd) PYTHON=python3 ./scripts/run_tests.sh --cov $PKG_NAME --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg

# Make sure repo is pip installable from git-archive zip
git archive -o /tmp/$PKG_NAME.zip HEAD
python3 -m pip install --force-reinstall /tmp/$PKG_NAME.zip
(cd /; python3 -c "from ${PKG_NAME} import get_include as gi; import os; assert 'gsl_odeiv2_cxx.pxd' in os.listdir(gi())")
