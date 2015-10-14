#!/bin/bash -e
cd $(dirname $0)/..
python2 setup.py build_ext -i
python2 -m pytest --ignore build/ --ignore doc/ --doctest-modules --pep8 --flakes $@
python3 setup.py build_ext -i
python3 -m pytest --ignore build/ --ignore doc/
