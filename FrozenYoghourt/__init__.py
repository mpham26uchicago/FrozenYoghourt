import numpy as np
from sympy import *
from sympy.physics.quantum import TensorProduct
from functools import reduce
from qiskit.quantum_info import *
from scipy.optimize import minimize
from numpy.testing import assert_array_almost_equal as aae
from typing import *
import matplotlib.pyplot as plt

def default_import(items:list = ['a', 'r', 'g', 'm', 'q', 'c', 'v']):
    
    """
    Print the import statement for the inputted modules
    
    Parameter
    ---------
    items: list, optional
        List of modules to import
    """
    
    import_statements = {'a': 'from FrozenYoghourt import *\n', 
                         'r': 'from FrozenYoghourt.mode import *\n', 
                         'g': 'from FrozenYoghourt.gates import *\n', 
                         'm': 'from FrozenYoghourt.maths import *\n', 
                         'q': 'from FrozenYoghourt.quantum import *\n', 
                         'c': 'from FrozenYoghourt.circuit import *\n', 
                         'v': 'from FrozenYoghourt.visualization import *\n'}
    
    
    import_str = ''
    for item in items:
        import_str += import_statements[item]
        
    print(import_str)
    
def view(mat, rounding = 10):
    
    """
    Display matrix in Latex format with the appropriate rounding
    
    Parameter
    ---------
    mat: numpy.ndarray, sympy.matrices.dense.MutableDenseMatrix
        Matrix to be viewed
    rounding: int, optional
        Number of decimal places to round off (default is 10)
    """
    
    if type(mat) == np.ndarray:
        display(Matrix(np.round(mat, rounding)))
    else:
        display(mat.evalf(rounding))
        
   