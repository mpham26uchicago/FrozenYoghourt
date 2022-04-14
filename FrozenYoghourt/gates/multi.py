from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.gates.single import *
from FrozenYoghourt.gates.double import *

def local_ops(num_qubits = 1, no_ops = 1, unimodular = False):
    """
    Parameters
    ----------
    num_qubits : int, optional
        The number of qubits.
        If not specified, the number of qubits is set to 1.
    no_ops : int, optional
        The number of local operators to generate.
        If not specified, the number of local operators is set to 1.
    unimodular : bool, optional
        If set to True, the local operators are generated with unimodular angles.
        If set to False, the local operators are generated with random angles.
        If not specified, the local operators are generated with random angles.
    
    Returns
    -------
    local_ops : list or numpy.ndarray or sympy.Matrix
        A list of local operators.
        If no_ops is 1, a single matrix is returned
        Each operator is a tensor product of "num_qubits" unitary matrices.
    """
    
    if type(no_ops) == bool:
        unimodular = no_ops; no_ops = 1
    
    if no_ops == 1:
        
        return tp(*[u2(unimodular, i) for i in range(num_qubits)])
    else:
        return [tp(*[u2(unimodular, i) for i in range(j*num_qubits, (j+1)*num_qubits)]) for j in range(no_ops)]

def pauli(val, mult = 1):
    
    """
    Return a Pauli matrix or tensor product of Pauli matrices
    
    Parameter
    ---------
    val: int, str, list
        Value or list of values of matrices in either integer (0, 1, 2, 3)
        or string ('i', 'x', 'y', 'z') form. If a list is given, compute 
        the tensor product of the pauli's in the list from left to right.
    
    mul: int
        Multiplicity value or list of multiplicity values
        
    Returns
    -------
    pauli_mat: np.ndarray or sympy.Matrix
        Pauli matrix or tensor product of Pauli matrices
        
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
    
    """
    Generate a multi-qubits swap matrix.
    
    Parameters
    ----------
    q1 : int
        The index of the first qubit to be swapped.
    q2 : int
        The index of the second qubit to be swapped.
    num_qubits : int
        The number of qubits in the system.

    Returns
    -------
    swap_matrix : np.ndarray or sympy.Matrix
        The swap matrix.
    """

    
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