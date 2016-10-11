from math import exp


def test_PyDecay():
    import pyximport
    pyximport.install()
    from _gsl_odeiv2_anyode import PyDecay

    pd = PyDecay(1.0)
    tout, yout = pd.adaptive(1.0, 1.0)
    for t, y in zip(tout, yout):
        assert abs(y - exp(-t)) < 2e-9


if __name__ == '__main__':
    test_PyDecay()
