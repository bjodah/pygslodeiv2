v0.3.4
======
- expose ``gslodeiv2.steppers`` tuple.

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
