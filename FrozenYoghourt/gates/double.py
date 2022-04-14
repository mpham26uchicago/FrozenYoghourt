from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.gates.single import *
from FrozenYoghourt.gates.multi import *

def CX():
    """
    Return the matrix of the CX (or CNOT) gate.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The matrix of the CX gate.    
    """
    
    if Mode.representation == 'numerical':
        return np.array([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 0, 1],
                         [0, 0, 1, 0]])
    else:
        return Matrix([[1, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 0, 1],
                       [0, 0, 1, 0]])
    
def CU(U):
    """
    Return the matrix of the control-unitary gate.
    
    Parameters
    ----------
    U : 2x2 matrix
        The unitary matrix to be controlled.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The matrix of the control-unitary gate   
    """
    
    if Mode.representation == 'numerical':
        return np.block([[ID(), np.zeros((2, 2))], 
          [np.zeros((2, 2)), U]])
    
    else:
        return BlockDiagMatrix(ID(), U).as_explicit()

def Gamma():
    """
    Return the matrix of the Gamma matrix for calculating KAK 
    angles from page 7 of https://arxiv.org/pdf/quant-ph/0507171.pdf

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The Gamma matrix.
    """
    
    if Mode.representation == 'numerical':
        return 1/2 * np.array([[1, -1, 1, -1], 
                               [1, 1, -1, -1], 
                               [1, -1, -1, 1], 
                               [1, 1, 1, 1]])
    else:
        return 1/2 * Matrix([[1, -1, 1, -1], 
                             [1, 1, -1, -1], 
                             [1, -1, -1, 1], 
                             [1, 1, 1, 1]])

def Magic():
    """
    Return the Magic matrix useful for KAK decomposition

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The Magic matrix 
    """
    
    if Mode.representation == 'numerical':
        return 1 / np.sqrt(2) * np.array([[1, 1j, 0, 0],
                                          [0, 0, 1j, 1],
                                          [0, 0, 1j, -1],
                                          [1, -1j, 0, 0]])
    else:
        return 1 / sqrt(2) * Matrix([[1, I, 0, 0],
                                     [0, 0, I, 1],
                                     [0, 0, I, -1],
                                     [1, -I, 0, 0]])

def Rxx(theta:float):
    """
    Return a 4x4 matrix that represents a rotation around the X⊗X-axis by an angle theta.

    Parameters
    ----------
    theta : float
        The angle of rotation.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 4x4 matrix that represents a rotation around the X⊗X-axis by an angle theta.
    """
    
    if Mode.representation == 'numerical':
        return np.array([[np.cos(theta / 2), 0, 0, -1j * np.sin(theta / 2)],
                         [0, np.cos(theta / 2), -1j * np.sin(theta / 2), 0],
                         [0, -1j * np.sin(theta / 2), np.cos(theta / 2), 0],
                         [-1j * np.sin(theta / 2), 0, 0, np.cos(theta / 2)]])
    else:
        return Matrix([[cos(theta / 2), 0, 0, -I * sin(theta / 2)],
                       [0, cos(theta / 2), -I * sin(theta / 2), 0],
                       [0, -I * sin(theta / 2), cos(theta / 2), 0],
                       [-I * sin(theta / 2), 0, 0, cos(theta / 2)]])

def Ryy(theta:float):
    """
    Return a 4x4 matrix that represents a rotation around the Y⊗Y-axis by an angle theta.

    Parameters
    ----------
    theta : float
        The angle of rotation.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 4x4 matrix that represents a rotation around the Y⊗Y-axis by an angle theta.
    """
    
    if Mode.representation == 'numerical':
        return np.array([[np.cos(theta / 2), 0, 0, 1j * np.sin(theta / 2)],
                         [0, np.cos(theta / 2), -1j * np.sin(theta / 2), 0],
                         [0, -1j * np.sin(theta / 2), np.cos(theta / 2), 0],
                         [1j * np.sin(theta / 2), 0, 0, np.cos(theta / 2)]])
    else:
        return Matrix([[cos(theta / 2), 0, 0, I * sin(theta / 2)],
                       [0, cos(theta / 2), -I * sin(theta / 2), 0],
                       [0, -I * sin(theta / 2), cos(theta / 2), 0],
                       [I * sin(theta / 2), 0, 0, cos(theta / 2)]])

def Rzz(theta:float):
    """
    Return a 4x4 matrix that represents a rotation around the Z⊗Z-axis by an angle theta.

    Parameters
    ----------
    theta : float
        The angle of rotation.

    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 4x4 matrix that represents a rotation around the Z⊗Z-axis by an angle theta.
    """
    
    if Mode.representation == 'numerical':
        return np.array([[np.exp(-1j * theta / 2), 0, 0, 0],
                         [0, np.exp(1j * theta / 2), 0, 0],
                         [0, 0, np.exp(1j * theta / 2), 0],
                         [0, 0, 0, np.exp(-1j * theta / 2)]])
    else:
        return Matrix([[exp(-I * theta / 2), 0, 0, 0],
                       [0, exp(I * theta / 2), 0, 0],
                       [0, 0, exp(I * theta / 2), 0],
                       [0, 0, 0, exp(-I * theta / 2)]])

def CAN(tx:float, ty:float, tz:float):
    
    """
    Return the two-qubits CAN (Canonical) gate
    
    Parameters
    ----------
    tx: float
        This is the angle of rotation about the X⊗X-axis.
    ty: float
        This is the angle of rotation about the Y⊗Y-axis.
    tz: float
        This is the angle of rotation about the Z⊗Z-axis.
    
    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The 4x4 matrix representing the two-qubits canonical gate.
        
    Notes
    -----
    The two-qubit Canonical gate is locally equivalent to all 4x4 
    unitary matrices. It is given as Rxx(tx)@Ryy(ty)@Rzz(tz).
            
    """
    
    if Mode.representation == 'numerical':
        return np.array([[np.exp(-1j*tz/2)*np.cos((tx-ty)/2), 0, 0, -1j*np.exp(-1j*tz/2)*np.sin((tx-ty)/2)], 
                         [0, np.exp(1j*tz/2)*np.cos((tx+ty)/2), -1j*np.exp(1j*tz/2)*np.sin((tx+ty)/2), 0], 
                         [0, -1j*np.exp(1j*tz/2)*np.sin((tx+ty)/2), np.exp(1j*tz/2)*np.cos((tx+ty)/2), 0], 
                         [-1j*np.exp(-1j*tz/2)*np.sin((tx-ty)/2), 0, 0, np.exp(-1j*tz/2)*np.cos((tx-ty)/2)]])
    else:
        return Matrix([[exp(-I*tz/2)*cos((tx-ty)/2), 0, 0, -I*exp(-I*tz/2)*sin((tx-ty)/2)], 
                       [0, exp(I*tz/2)*cos((tx+ty)/2), -I*exp(I*tz/2)*sin((tx+ty)/2), 0], 
                       [0, -I*exp(I*tz/2)*sin((tx+ty)/2), exp(I*tz/2)*cos((tx+ty)/2), 0], 
                       [-I*exp(-I*tz/2)*sin((tx-ty)/2), 0, 0, exp(-I*tz/2)*cos((tx-ty)/2)]])