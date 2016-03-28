#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Tested with GSL v1.16 and v2.1

import os
import shutil
import sys
from setuptools import setup
from setuptools.extension import Extension


pkg_name = 'pygslodeiv2'


def _path_under_setup(*args):
    return os.path.join(os.path.dirname(__file__), *args)

USE_CYTHON = os.path.exists(_path_under_setup(
    'pygslodeiv2', '_gslodeiv2_numpy.pyx'))

# Cythonize .pyx file if it exists (not in source distribution)
ext_modules = []

if len(sys.argv) > 1 and '--help' not in sys.argv[1:] and sys.argv[1] not in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    import numpy as np
    ext = '.pyx' if USE_CYTHON else '.cpp'
    ext_modules = [
        Extension('pygslodeiv2._gslodeiv2_numpy',
                  ['pygslodeiv2/_gslodeiv2_numpy'+ext],
                  language='c++', extra_compile_args=['-std=c++11'],
                  libraries=['gsl', 'gslcblas', 'm'],
                  include_dirs=['./include', np.get_include()])
    ]
    if USE_CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(
            ext_modules, include_path=[_path_under_setup('include')],
            gdb_debug=True)

RELEASE_VERSION = os.environ.get('PYGSLODEIV2_RELEASE_VERSION', '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
CONDA_BUILD = os.environ.get('CONDA_BUILD', '0') == '1'
if CONDA_BUILD:
    try:
        RELEASE_VERSION = 'v' + open(
            '__conda_version__.txt', 'rt').readline().rstrip()
    except IOError:
        pass

release_py_path = _path_under_setup(pkg_name, '_release.py')

if (len(RELEASE_VERSION) > 1 and
   RELEASE_VERSION[0] == 'v'):
    TAGGED_RELEASE = True
    __version__ = RELEASE_VERSION[1:]
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

with open(_path_under_setup(pkg_name, '__init__.py'), 'rt') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
assert 10 < len(short_description) < 255
long_description = open(_path_under_setup('README.rst')).read()
assert len(long_description) > 100

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=short_description,
    long_description=long_description,
    classifiers=classifiers,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://github.com/bjodah/' + pkg_name,
    license='GPLv3',
    packages=[pkg_name] + tests,
    install_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    setup_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    ext_modules=ext_modules,
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
