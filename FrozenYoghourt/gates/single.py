from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.gates.double import *
from FrozenYoghourt.gates.multi import *

def ID():
    """
    Returns the identity matrix.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The identity matrix.    
    """
    
    if Mode.representation == 'numerical':
        return np.array([[1, 0], [0, 1]])
    else:
        return Matrix([[1, 0], [0, 1]])

def X():
    """
    Return the Pauli-X matrix.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The Pauli-X matrix.    
    """
    
    if Mode.representation == 'numerical':
        return np.array([[0, 1], [1, 0]])
    else:
        return Matrix([[0, 1], [1, 0]])

def Y():
    """
    Return the Pauli-Y matrix.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The Pauli-Y matrix.    
    """
    
    if Mode.representation == 'numerical':
        return np.array([[0, -1j], [1j, 0]])
    else:
        return Matrix([[0, -I], [I, 0]])

def Z():
    """
    Return the Pauli-Z matrix.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The Pauli-Z matrix.    
    """
    
    if Mode.representation == 'numerical':
        return np.array([[1, 0], [0, -1]])
    else:
        return Matrix([[1, 0], [0, -1]])
    
def H():
    """
    Return the Hadamard matrix.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The Hadamard matrix.    
    """
    
    if Mode.representation == 'numerical':
        return 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]])
    else:
        return 1 / sqrt(2) * Matrix([[1, 1], [1, -1]])

def Rx(theta:float):
    """
    Return a 2x2 matrix that represents a rotation around the x-axis by an angle theta.

    Parameters
    ----------
    theta : float
        The angle of rotation.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 2x2 matrix that represents a rotation around the x-axis by an angle theta.
    """
    
    if Mode.representation == 'numerical':
        return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                         [-1j * np.sin(theta / 2), np.cos(theta / 2)]])
    else:
        return Matrix([[cos(theta / 2), -I * sin(theta / 2)],
                       [-I * sin(theta / 2), cos(theta / 2)]])

def Ry(theta:float):
    """
    Return a 2x2 matrix that represents a rotation around the y-axis by an angle theta.

    Parameters
    ----------
    theta : float
        The angle of rotation.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 2x2 matrix that represents a rotation around the y-axis by an angle theta.
    """
    
    if Mode.representation == 'numerical':
        return np.array([[np.cos(theta / 2), -np.sin(theta / 2)],
                         [np.sin(theta / 2), np.cos(theta / 2)]])
    else:
        return Matrix([[cos(theta / 2), -sin(theta / 2)],
                       [sin(theta / 2), cos(theta / 2)]])

def Rz(theta:float):
    """
    Return a 2x2 matrix that represents a rotation around the z-axis by an angle theta.

    Parameters
    ----------
    theta : float
        The angle of rotation.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 2x2 matrix that represents a rotation around the z-axis by an angle theta.
    """
    
    if Mode.representation == 'numerical':
        return np.array([[np.exp(-1j * theta / 2), 0],
                         [0, np.exp(1j * theta / 2)]])
    else:
        return Matrix([[exp(-I * theta / 2), 0],
                       [0, exp(I * theta / 2)]])

def Phase(theta:float):
    """
    Return a 2x2 matrix that represents a phase shift of by an angle theta.

    Parameters
    ----------
    theta : float
        The angle of rotation.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 2x2 matrix that represents a phase shift an angle theta.
    """
    
    if Mode.representation == 'numerical':
        return np.array([[1, 0], [0, np.exp(1j*theta)]])
    else:
        return Matrix([[1, 0], [0, exp(I*theta)]])
    
def u2(angles:list = None, unimodular = False, index = ''):
    """
    Generates a random 2x2 unitary matrix.
    
    Parameters
    ----------
    angles : list, optional
        A list of three angles, theta, phi, and lambda, by which the matrix is generated.
        If not specified, the angles are generated randomly and are returned as symbols.
    unimodular : bool, optional
        If set to True, the matrix is generated with unimodular angles.
        If set to False, the matrix is generated with random angles.
        If not specified, the matrix is generated with random angles.
    index : str, optional
        A string that is added to the end of the names of the angles.
        If not specified, the angles are named theta, phi, and lambda.
    
    Returns
    -------
    U : numpy.ndarray or sympy.Matrix
        A 2x2 unitary matrix.
    """
    
    if type(angles) == bool:
        unimodular = angles
        angles = None
    
    if Mode.representation == 'numerical':
        if unimodular: 
            return to_su(random_unitary(2).data)
        else: 
            return random_unitary(2).data
    else:
        if angles is None:
            theta, phi, lam = symbols(f'theta_{index}, phi_{index}, lambda_{index}', real = True)
        else:
            theta, phi, lam = angles
        if unimodular:
            return Matrix([[exp(I*phi)*cos(theta / 2), -exp(I * -lam) * sin(theta / 2)],
                                          [exp(I * lam) * sin(theta / 2), exp(I * (-phi)) * cos(theta / 2)]])
        else:
            return Matrix([[cos(theta / 2), -exp(I * lam) * sin(theta / 2)],
                                           [exp(I * phi) * sin(theta / 2), exp(I * (phi + lam)) * cos(theta / 2)]])