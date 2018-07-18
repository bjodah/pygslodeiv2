import os

env = {
    'BLAS': 'gslcblas',
    'GSL_LIBS': 'gsl' if os.name == 'nt' else 'gsl,m'
}
