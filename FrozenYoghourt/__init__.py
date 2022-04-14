import numpy as np
import scipy
from sympy import *
from sympy.physics.quantum import TensorProduct
from functools import reduce
from qiskit.quantum_info import *
from scipy.optimize import minimize
from numpy.testing import assert_array_almost_equal as aae
from typing import *
import matplotlib.pyplot as plt
import pyperclip

def default_import(copy = True, ):
    """Just import everything"""
    
    import_lst = [
    "from FrozenYoghourt import *",
    "from FrozenYoghourt import *",
    "from FrozenYoghourt.mode import *",
    "from FrozenYoghourt.gates import *",
    "from FrozenYoghourt.gates.single import *",
    "from FrozenYoghourt.gates.double import *",
    "from FrozenYoghourt.gates.multi import *",
    "from FrozenYoghourt.maths import *",
    "from FrozenYoghourt.quantum import *",
    "from FrozenYoghourt.quantum.decomposition import *",
    "from FrozenYoghourt.circuit import *",
    "from FrozenYoghourt.visualization import *",
    ]
    
    import_str = '\n'.join(import_lst)
    
    if copy:
        pyperclip.copy(import_str)
        print('\033[1;34mCopied to clipboard!\033[1;32m')
    else:
        print(import_str)
    
def view(mat, rounding = 10):
    """
    Display matrix in Latex format with the appropriate rounding
    
    Parameters
    ---------
    mat: numpy.ndarray, sympy.matrices.dense.MutableDenseMatrix
        Matrix to be viewed
    rounding: int, optional
        Number of decimal places to round off (default is 10)
    """
    
    if type(mat) == np.ndarray:
        display(Matrix(np.round(mat, rounding)))
    else:
        display(Matrix(np.round(np.array(mat.evalf()).astype(complex), rounding)))