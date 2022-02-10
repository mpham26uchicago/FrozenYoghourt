import numpy as np
from sympy import *
from sympy.physics.quantum import TensorProduct
from functools import reduce
from qiskit.quantum_info import *
from scipy.optimize import minimize
from numpy.testing import assert_array_almost_equal as aae


def view(mat, rounding = 10):
    display(Matrix(np.round(mat, rounding)))

class Mode:
    representation = 'numpy'

    @classmethod
    def toggle(cls):
        cls.representation = {'numpy': 'sympy',
                              'sympy': 'numpy'}[cls.representation]
        print(f'\033[1m{cls.representation.capitalize()}\033[0m mode activated')

    @classmethod
    def now(cls):
        print(f'\033[1m{cls.representation.capitalize()}\033[0m mode is \033[1mon\033[0m')
