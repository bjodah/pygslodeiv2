#!/bin/bash
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
sphinx-apidoc --full --force --doc-author="$(head -n 1 AUTHORS)" --doc-version=$(python setup.py --version) -F -o doc $PKG/
sed -i 's/Contents/.. include:: ..\/README.rst\n\nContents/g' doc/index.rst
( cd doc; make html )
