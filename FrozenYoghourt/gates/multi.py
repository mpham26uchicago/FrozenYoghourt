from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.single import *
from FrozenYoghourt.double import *

def local_ops(num_qubits = 1, no_ops = 1, unimodular = False):
    
    if type(no_ops) == bool:
        unimodular = no_ops; no_ops = 1
    
    if no_ops == 1:
        return tp(*[u2(unimodular, i) for i in range(num_qubits)])
    else:
        return [tp(*[u2(unimodular, i) for i in range(j*num_qubits, (j+1)*num_qubits)]) for j in range(no_ops)]

def pauli(val, mult = 1):
    
    """
    Return a Pauli matrix or tensor product of matrices
    
    Parameter
    ---------
    val: int, str, list
        Value or list of values of matrices in either integer (0, 1, 2, 3)
        or string ('i', 'x', 'y', 'z') form. If a list is given, compute 
        the tensor product of the pauli's in the list from left to right.
    
    mul: int
        Multiplicity of value or list of values
        
    Returns
    -------
    pauli_mat: np.ndarray, sympy.Matrices
        Pauli matrix or tensor product of matrices
        
    """
    
    option = type(val)
    if option == list:
        option = type(val[0])
    
    if option == int:
        num = {0: ID(), 1: X(), 2: Y(), 3:Z()}
        if type(val) == int:
            pauli_mat = tp(num[val], mul = mul)
        else:
            pauli_list = [num[elem] for elem in val]
            pauli_mat = tp(*pauli_list, mul = mul)
      
    elif option == str:
        text = {'i': ID(), 'x': X(), 'y': Y(), 'z': Z()}
        if type(val) == str:
            pauli_mat = tp(text[val], mult = mult)
        else:
            pauli_list = [text[elem] for elem in val]
            pauli_mat = tp(*pauli_list, mult = mult)
            
    return pauli_mat

def Swap(q1, q2, num_qubits):
    
    order_list = []

    for k in range(2**num_qubits):
        bin_lst = list(bin(k)[2:].zfill(num_qubits))
        bin_lst[q1], bin_lst[q2] = bin_lst[q2], bin_lst[q1]
        order_list.append(int(''.join(bin_lst), 2))

    if Mode.representation == 'numerical':   
        swap_matrix = np.eye(2**num_qubits)[:, order_list]
    else:
        swap_matrix = eye(2**num_qubits)[:, order_list]
        
    return swap_matrix