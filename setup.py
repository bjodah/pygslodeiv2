#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Tested with GSL v1.16 and v2.1

import io
import os
import shutil
import sys
from setuptools import setup
from setuptools.extension import Extension


pkg_name = 'pygslodeiv2'


def _path_under_setup(*args):
    return os.path.join(os.path.dirname(__file__), *args)

USE_CYTHON = os.path.exists(_path_under_setup(
    pkg_name, '_gsl_odeiv2.pyx'))
package_include = os.path.join(pkg_name, 'include')

# Cythonize .pyx file if it exists (not in source distribution)
ext_modules = []

if len(sys.argv) > 1 and '--help' not in sys.argv[1:] and sys.argv[1] not in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    import numpy as np
    ext = '.pyx' if USE_CYTHON else '.cpp'
    sources = [os.path.join('pygslodeiv2', '_gsl_odeiv2' + ext)]
    ext_modules = [Extension('%s._gsl_odeiv2' % pkg_name, sources)]
    if USE_CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules, include_path=[
            package_include,
            os.path.join('external', 'anyode', 'cython_def')
        ])
    ext_modules[0].language = 'c++'
    ext_modules[0].extra_compile_args = ['-std=c++11']
    ext_modules[0].include_dirs = [
        np.get_include(), package_include,
        os.path.join('external', 'anyode', 'include')]
    ext_modules[0].libraries += ['gsl', 'gslcblas', 'm']

_version_env_var = '%s_RELEASE_VERSION' % pkg_name.upper()
RELEASE_VERSION = os.environ.get(_version_env_var, '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
if os.environ.get('CONDA_BUILD', '0') == '1':
    try:
        RELEASE_VERSION = 'v' + open(
            '__conda_version__.txt', 'rt').readline().rstrip()
    except IOError:
        pass

release_py_path = _path_under_setup(pkg_name, '_release.py')

if len(RELEASE_VERSION) > 1:
    if RELEASE_VERSION[0] != 'v':
        raise ValueError("$%s does not start with 'v'" % _version_env_var)
    TAGGED_RELEASE = True
    __version__ = RELEASE_VERSION[1:]
else:  # read __version__ attribute from _release.py:
    TAGGED_RELEASE = False
    exec(open(release_py_path).read())

classifiers = [
    "Development Status :: 4 - Beta",
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
]

tests = [
    '%s.tests' % pkg_name,
]

with io.open(_path_under_setup(pkg_name, '__init__.py'),
             'rt', encoding='utf-8') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
assert 10 < len(short_description) < 255
long_descr = io.open(_path_under_setup('README.rst'), encoding='utf-8').read()
assert len(long_descr) > 100

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=short_description,
    long_description=long_descr,
    classifiers=classifiers,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://github.com/bjodah/' + pkg_name,
    license='GPLv3',
    packages=[pkg_name] + tests,
    include_package_data=True,
    install_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    setup_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    extras_require={'docs': ['Sphinx', 'sphinx_rtd_theme']},
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
