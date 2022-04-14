from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.gates import *

def epsilon(psi, num_qubits=2):
    """
    Calculate the epsilon parameter of a given state
    
    Parameters
    ----------
    psi : ndarray
        The state vector of the system.
    num_qubits : int
        The number of qubits in the system.

    Returns
    -------
    float
        The epsilon parameter of the system.

    Notes
    -----
    The epsilon parameter is a measure of the entanglement of the system.
    It is defined as the absolute value of the inner product of the state
    vector with the tensor product of the Pauli Y operator with itself
    "num_qubits" times.
    """
    
    if Mode.representation == 'numerical':
        return np.abs(psi.T @ tp(Y(), no_times=num_qubits) @ psi)[0]
    else:
        return Abs(psi.T @ tp(Y(), no_times=num_qubits) @ psi)[0]

def global_phase(A, B):
    """
    This function checks if two matrices A and B are equivalent up to a global phase.

    Parameters
    ----------
    A : numpy.ndarray
        A 2D numpy array of complex numbers.
    B : numpy.ndarray
        A 2D numpy array of complex numbers.

    Returns
    -------
    numpy.ndarray
        A 1D numpy array of complex numbers, containing the global phase factor.
        If A and B are not equivalent up to a global phase, the function returns False.
    """
    
    D = np.diag(np.conj(A).T @ B)

    if np.isclose(np.min(D), np.max(D)):
        return D[0]
    else:
        print('A and B are not equivalent up to global phase')
        return False

def ymap(U, unimodular = False):
    """
    This function returns the Y-map of a given unitary matrix U. The Y-map is used
    to calculate the double coset of a 4x4 unitary matrix 
    (Page 3 https://web.eecs.umich.edu/~imarkov/pubs/jour/pra04-univ.pdf).
    
    Parameters
    ----------
    U : np.ndarray or sympy.Matrix
        The unitary matrix to be mapped.
    unimodular : bool, optional
        If True, the function will convert U to a special unitary matrix before applying the Y-map.
        The default is False.
    
    Returns
    -------
    np.ndarray or sympy.Matrix
        The Y-map of U.
    """
    
    n = int(np.log2(U.shape[0]))

    if (Mode.representation == 'numerical') and (not np.isclose(np.linalg.det(U), 1)):
        U = to_su(U)
    elif unimodular:
        U = to_su(U)
        
    E = tp(Y(), no_times=n)
    return U @ E @ U.T @ E

def xmap(M, variable=None):
    """
    This function returns the coefficients of the characteristic
    polynomial of a given matrix. The name xmap comes from the 
    "chi invariant" on Page 2 https://arxiv.org/pdf/quant-ph/0308045.pdf.
    
    Parameters
    ----------
    M: np.ndarray or sympy.Matrix
        The unitary matrix.
    variable: None, optional
        Variable to represent the characteristic polynomial.
        
    Returns
    -------
    If mode is "numerical":
    coef: np.ndarray
        Return an array of the characteristic coefficients.
    If mode is "symbolic"
    charpoly: sympy.Poly
        Return the characteristic polynomial. 
    """
    
    dim = M.shape[0]

    coef = np.array([1])
    Mk = np.array(M)

    for k in range(1, dim + 1):
        ak = -Mk.trace() / k
        coef = np.append(coef, ak)
        Mk += np.diag(np.repeat(ak, dim))
        Mk = np.dot(M, Mk)

    if Mode.representation == 'numerical':
        return coef

    else:
        if variable is None:
            return Matrix(coef)
        else:
            var = symbols(variable)
            variable_mat = Matrix([var ** n for n in range(dim + 1)])
            coef = Matrix(coef).T
            char_poly = simplify(coef @ variable_mat)[0]
            return char_poly

def shende_invariant(U, V = None, return_charpoly = False):
    
    """
    Calculate the Shende invariant for a 4x4 unitary matrix U
    given on page 3 of https://web.eecs.umich.edu/~imarkov/pubs/jour/pra04-univ.pdf
    
    Parameters:
    -----------
    U : numpy.ndarray
        A square matrix of size n.
    V : numpy.ndarray
        A square matrix of size n.
    return_charpoly : bool
        If True, the Shende invariant of U and V are returned, 
        even if they are not equal.
        
    Returns:
    --------
    If return_charpoly is False:
        equivalent (bool): If True, V and U are equivalent.
        If False, V and U are not equivalent.

    If return_charpoly is True:
        equivalent (bool): If True, V and U are equivalent.
        If False, V and U are not equivalent.
        charpoly_U (np.matrix): The Shende invariant of U.
        charpoly_V (np.matrix): The Shende invariant of V.
        
    """

    if V is None:
        return xmap(ymap(to_su(U)))
    
    elif Mode.representation == 'numerical':
        charpoly_U = xmap(ymap(to_su(U)))
        charpoly_V = xmap(ymap(to_su(V)))
        
        charpoly_equiv = np.allclose(charpoly_U, charpoly_V)
        
        if return_charpoly:
            return charpoly_equiv, charpoly_U, charpoly_V
        else:
            return charpoly_equiv
        
    else: ## This needs fixing because the determinant function doesn't work
        charpoly_U = xmap(ymap(to_su(U)))
        charpoly_V = xmap(ymap(to_su(V)))
        
        charpoly_equiv = charpoly_U.equals(charpoly_V)
        
        if return_charpoly:
            return charpoly_equiv, charpoly_U, charpoly_V
        else:
            return charpoly_equiv

def huang_invariant(U):
    
    """
    Compute the Huang invariant of a given SU(4) matrix. This invariant
    is used to determine the whether two matrices are locally equivalent.
    It is described on page 7 of https://arxiv.org/pdf/2105.06074.pdf.
    
    Parameters
    ----------
    U : np.ndarray or sympy.Matrix
        An SU(4) matrix.
    
    Returns
    -------
    C, B, A : float
        Coefficients of the characteristic quadratic
    """
        
    V = dagger(Magic())@to_su(U)@Magic()
    
    if Mode.representation == 'numerical':    
        D,O = np.linalg.eig(V@V.T)
        spectrum = np.angle(D)
        spectrum[3] -= np.sum(spectrum)
        
        C = np.prod(np.sin(1/2 * spectrum))
        B = -1/4 * np.sum(np.sin(spectrum))
        A = np.prod(np.cos(1/2 * spectrum))
        
        return C, B, A

    else:
        Real, Imag = Matrix(V).as_real_imag()
        t = symbols('t')
        poly = det(Real + t*Imag)
        character_polynomial, _ = div(poly.args[1], 1/poly.args[0])
        
        return character_polynomial

def is_local(U):
    
    """
    This function checks if a given unitary is local.

    Parameters
    -----------
    U : np.array
        A 4x4 unitary matrix

    Returns
    -------
    local : bool 
        True if the unitary is local, False otherwise
    """
    
    M = Magic()
    V = dagger(M)@to_su(U)@M
    if np.isclose(np.linalg.det(V), 1) and np.all(np.isclose(V@V.T, np.identity(4))):
        return True
    else: 
        return False
    
def kron_decomp(U:np.ndarray):
    """
    Given a 4x4 unitary matrix, this function returns two 2x2 
    unitary matrices, a and b, such that: U = a ⊗ b
    
    Parameters
    ----------
    U : np.ndarray
        4x4 numpy array representing a unitary matrix
    
    Returns
    -------
    a : np.ndarray
        2x2 numpy array representing a unitary matrix
    b : np.ndarray
        2x2 numpy array representing a unitary matrix
    """
    
    assert is_local(U), 'Matrix is not local' # Check if gates is local
    
    m = U.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)

    u, sv, vh = np.linalg.svd(m)

    a = np.sqrt(sv[0]) * u[:, 0].reshape(2, 2)
    b = np.sqrt(sv[0]) * vh[0, :].reshape(2, 2)

    return a, b

def KAK(U:np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:

    """
    This function is an implementation of the KAK decomposition. 
    It takes a 4x4 unitary matrix U and decomposes it into local
    gates wrapped around a Canonical gate
    
    Parameters
    ----------
    U : np.ndarray
        A 4x4 unitary matrix
    
    Returns
    -------
    K2 : np.ndarray
        A 4x4 local gate
    A : np.ndarray
        A 4x4 canonical matrix
    K1 : np.ndarray
        A 4x4 local gate
    phase : float
        A global phase
    """
    
    M = Magic()
    
    phase = global_phase(to_su(U), U)
    V = dagger(M)@to_su(U)@M
    D_squared, P = np.linalg.eig(V.T@V)

    spectrum = np.angle(D_squared) # These steps make sure that D stays unimodular after the square root
    spectrum[3] -= np.sum(spectrum)
    spectrum = (spectrum/2).reshape(4, )

    t_list = Gamma().T@spectrum

    D = np.diag(np.exp(1j*spectrum))
    
    if np.isclose(np.linalg.det(P), -1): # If the determinant of P is -1, turn it to 1
        P[:, 0] = -P[:, 0]
    
    Q_1, Q_2 = P.T, V@P@np.conj(D).T # Extract the outer orthogonal gates

    A = np.exp(1j*t_list[0]/2) * CAN(t_list[1], t_list[2], t_list[3]) # Construct canonical matrix

    K1, K2 = M@Q_1@dagger(M), M@Q_2@dagger(M) # Extract local gates

    return K2, A, K1, phase

def bra(x: Union[str, int], num_qubits = None):
    """
    Returns a bra vector representing the state of a quantum system.
    
    Parameters
    ----------
    x : int or str
        The integer value or binary string representation of the state.
    num_qubits : int, optional
        The number of qubits in the state. If `x` is a string, `num_qubits` is
        not required. If None, defaults to the length of the binary string
        
    Returns
    -------
    full_bra : np.ndarray or sympy.Matrix
        The bra vector representation of the state.
    
    Examples
    --------
    >>> bra('00')
    array([[1., 0., 0., 0.]])
    >>> bra(0, 2)
    array([[1., 0., 0., 0.]])
    >>> bra(0)
    array([[1., 0.]])
    """
    
    if isinstance(x, str):
        num_qubits = len(x)
    elif num_qubits is None:
        num_qubits = len(bin(x))-2
        
    x_val = int(x)
    
    if Mode.representation == 'numerical':
        full_bra = np.zeros((1, 2**num_qubits))
    else:
        full_bra = zeros(1, 2**num_qubits)
        
    full_bra[0, x_val] = 1
    
    return full_bra 

def ket(x: Union[str, int], num_qubits = None):
    
    """
    Returns a ket vector representing the state of a quantum system.
    
    Parameters
    ----------
    x : int or str
        The integer value or binary string representation of the state.
    num_qubits : int, optional
        The number of qubits in the state. If `x` is a string, `num_qubits` is
        not required. If None, defaults to the length of the binary string
        
    Returns
    -------
    full_ket : np.ndarray or sympy.Matrix
        The ket vector representation of the state.
    """
    
    if isinstance(x, str):
        num_qubits = len(x)
    elif num_qubits is None:
        num_qubits = len(bin(x))-2
        
    x_val = int(x)
    
    if Mode.representation == 'numerical':
        full_ket = np.zeros((2**num_qubits, 1))
    else:
        full_ket = zeros(2**num_qubits, 1)
        
    full_ket[x_val, 0] = 1
    
    return full_ket