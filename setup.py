#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Tested with GSL 1.16

import os
import shutil
import sys
from distutils.core import setup
from distutils.extension import Extension
import numpy as np


pkg_name = 'pygslodeiv2'

# Cythonize .pyx file if it exists (not in source distribution)
ext_modules = []

if len(sys.argv) > 1 and '--help' not in sys.argv[1:] and sys.argv[1] not in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    USE_CYTHON = os.path.exists('pygslodeiv2/_gslodeiv2_numpy.pyx')
    ext = '.pyx' if USE_CYTHON else '.cpp'
    ext_modules = [
        Extension('pygslodeiv2._gslodeiv2_numpy',
                  ['pygslodeiv2/_gslodeiv2_numpy'+ext],
                  language='c++', extra_compile_args=['-std=c++11'],
                  libraries=['gsl', 'gslcblas', 'm'])
    ]
    if USE_CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules, include_path=['./include'],
                                gdb_debug=True)

PYGSLODEIV2_RELEASE_VERSION = os.environ.get('PYGSLODEIV2_RELEASE_VERSION', '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
CONDA_BUILD = os.environ.get('CONDA_BUILD', '0') == '1'
if CONDA_BUILD:
    try:
        PYGSLODEIV2_RELEASE_VERSION = 'v' + open(
            '__conda_version__.txt', 'rt').readline().rstrip()
    except IOError:
        pass

release_py_path = os.path.join(pkg_name, '_release.py')

if (len(PYGSLODEIV2_RELEASE_VERSION) > 1 and
   PYGSLODEIV2_RELEASE_VERSION[0] == 'v'):
    TAGGED_RELEASE = True
    __version__ = PYGSLODEIV2_RELEASE_VERSION[1:]
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(open(release_py_path).read())

classifiers = [
    "Development Status :: 3 - Alpha",
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
]

tests = [
    'pygslodeiv2.tests',
]

descr = 'Python binding for odeiv2 in GNU Scientific Library (GSL).'
setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=descr,
    classifiers=classifiers,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://github.com/bjodah/' + pkg_name,
    license='GPLv3',
    packages=[pkg_name] + tests,
    ext_modules=ext_modules,
    include_dirs=[np.get_include(), './include']
)

if __name__ == '__main__':
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist
            # depending on tagged version (set PYGSLODEIV2_RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write(
                "__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__', release_py_path)
