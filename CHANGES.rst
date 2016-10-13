v0.6.0
======
- Changed development status from alpha to beta.
- Refactored to use AnyODE base class (share code with pycvodes & pygslodeiv2)

v0.5.2
======
- Fixes to setu.py

v0.5.1
======
- ``nsteps`` kwarg added (maximum number of steps).
- More robust setup.py

v0.5.0
======
- Changes to info dict: rename 'nrhs' -> 'nfev', 'njac' -> 'njev', added 'cpu_time', 'success'

v0.4.1
======
- Added support for (first) derivative in output
- Min and max step now allowed to be set
- Check against using dx0=0.0

v0.4.0
======
- New function signature: integrate_predefined and integrate_adaptive now
  also return an info dict containing ``nrhs`` and ``njac`` conatining
  number of calls to each function made during last integration.
- Expose ``gslodeiv2.steppers`` tuple.
- check_callbable and check_indexing kwargs now defaults to False

v0.3.3
======
- Fix minor memory leak
- Made y read-only in Python callbacks
- Do not overwrite Python error string when callback raises Exception.

v0.3.2
======
- Ship tests with package (e.g.: python -m pytest --pyargs pygslodeiv2)

v0.3.1
======
- Less strict callback checks on python side.
- Minor C++ API clean up.


v0.3.0
======
- Jacobian callback only need for steppers using it.

v0.2.0
======
- integrate_predefined added. More extensive tesing of steppers.

v0.1
====
- Integration using adaptive step-size supported.
