from FrozenYoghourt import *
from FrozenYoghourt.mode import *

def mm(*lst, mult=1):
    """
    Multiplies a list of matrices together.

    Parameters
    ----------
    lst : list
        A list of matrices to multiply together.
    mult : int
        The number of times to multiply the list of matrices together.

    Returns
    -------
    numpy.ndarray
        The product of the matrices.
    """
    lst *= mult
    return reduce(lambda a, b: a @ b, lst)

def tp(*lst, mult=1):
    """
    Returns the tensor product of the arguments.
    If mult is not 1, the tensor product is repeated "mult" times.
    If the representation is numerical, the tensor product is calculated using numpy.
    If the representation is symbolic, the tensor product is calculated using sympy.
    
    Parameters
    ----------
    lst : list
        A list of matrices or vectors.
    mult : int
        The number of times to repeat the tensor product.
        
    Returns
    -------
    numpy.ndarray or sympy.Matrix
        The tensor product of the arguments.
    """
    lst *= mult
    if Mode.representation == 'numerical':
        return reduce(lambda a, b: np.kron(a, b), lst)
    else:
        return TensorProduct(*lst)

def dagger(*u:Union[Tuple[np.ndarray], Tuple[Matrix]]):
    """
    Returns the conjugate transpose of a matrix or a list of matrices.
    
    Parameters
    ----------
    u : np.ndarray or sympy.Matrix
        A matrix or a list of matrices.
    
    Returns
    -------
    np.ndarray or sympy.Matrix
        The conjugate transpose of the input matrix or list of matrices.
    """
    
    if Mode.representation == 'numerical':
        if len(u) == 1:
            return np.conj(u[0].T)
        else:
            dagger_list = [np.conj(mat.T) for mat in u]
            return dagger_list
    else:
        if len(u) == 1:
            return u[0].T.conjugate()
        else:
            dagger_list = [mat.T.conjugate() for mat in u]
            return dagger_list
    
def svd(A:Union[np.ndarray, Matrix]):
    """
    Computes the singular value decomposition of a matrix.

    Parameters
    ----------
    A : np.ndarray or sympy.Matrix
        The matrix to decompose.

    Returns
    -------
    U : np.ndarray or sympy.Matrix
        The left singular unitary matrix

    S : np.ndarray or sympy.Matrix
        The singular diagonal matrix

    V : Manp.ndarray or sympy.Matrixtrix
        The right singular unitary matrix
    """
    
    if Mode.representation == 'numerical':
        return np.linalg.svd(A)
    else:
        AH = A.H
        m, n = A.shape
        if m >= n:
            V, S = (AH * A).diagonalize()

            S_diag = [S[i, i] for i in range(S.rows)]
            ranked = []
            for i, x in enumerate(S.diagonal()):
                if not x.is_zero:
                    ranked.append(i)

            V = V[:, ranked]
            S = Matrix.diag([sqrt(S[i, i]) for i in range(S.rows) if i in ranked])

            V, _ = V.QRdecomposition()
            U = A * V * S.inv()
        else:
            U, S = (A * AH).diagonalize()

            ranked = []
            for i, x in enumerate(S.diagonal()):
                if not x.is_zero:
                    ranked.append(i)

            U = U[:, ranked]
            S = Matrix.diag([sqrt(S[i, i]) for i in range(S.rows) if i in ranked])

            U, _ = U.QRdecomposition()
            V = AH * U * S.inv()

        return U, S, V

def to_su(*u:Union[Tuple[np.ndarray], Tuple[Matrix]]):
    """
    Converts a matrix or a list of matrices to special unitary form.

    Parameters
    ----------
    u : np.ndarray or sympy.Matrix
        A matrix or a list of matrices to be converted to special unitary form.

    Returns
    -------
    np.ndarray or sympy.Matrix
        The matrix or the list of matrices in special unitary form.
    """
    
    if Mode.representation == 'numerical':
        if len(u) == 1:
            return u[0] * complex(np.linalg.det(u[0])) ** (-1 / np.shape(u[0])[0])
        else:
            to_su_list = [mat * complex(np.linalg.det(mat)) ** (-1 / np.shape(mat)[0]) for mat in u]
            return to_su_list
    else:
        if len(u) == 1:
            return (u[0] * complex(u[0].det()) ** (-1 / u[0].shape[0])).evalf()
        else:
            to_su_list = [(mat * complex(mat.det()) ** (-1 / mat.shape[0])).evalf() for mat in u]
            return to_su_list

def fast_substitution(matrix:Matrix, variables, values, to_numpy=False):
    """
    Substitutes values into a matrix.
    
    Parameters
    ----------
    matrix : sympy.Matrix
        The matrix to substitute values into.
    variables : list or tuple or Symbol
        The variables to substitute.
    values : list or tuple or float
        The values to substitute.
    to_numpy : bool
        Whether to convert the matrix to a numpy array.
        
    Returns
    -------
    sympy.Matrix or numpy.ndarray
        The substituted matrix.
    """
    
    if (type(variables) == list) or (type(variables) == tuple):
        substituted_matrix = matrix.subs(list(zip(variables, values)))
    else:
        substituted_matrix = matrix.subs(variables, values)
        
    if to_numpy:
        return np.array(substituted_matrix).astype(complex)
    else:
        return substituted_matrix.evalf()

def optimizing_identities(lhs_fun, rhs_fun, no_var, var_split, 
                          random=(0, 2 * np.pi), num_qubits=2, error=1e-7):
    """
    This function is used to optimize the parameters of a given identity.

    Parameters
    ----------
    lhs_fun : function
        The function that returns the left hand side of the identity.
    rhs_fun : function
        The function that returns the right hand side of the identity.
    no_var : int
        The total number of variables in the identity.
    var_split : int
        The number of variables in the left hand side of the identity.
    random : tuple, optional
        The range of the random numbers to be generated. The default is (0, 2 * np.pi).
    num_qubits : int, optional
        The number of qubits in the identity. The default is 2.
    error : float, optional
        The error tolerance for the optimization. The default is 1e-7.
    
    Returns
    -------
    job : scipy.optimize.optimize.OptimizeResult
        The result of the optimization.
    """
    
    compute_rhs = True
    if var_split == no_var:
        compute_rhs = False

    def cost(x, compute_rhs):
        LHS = lhs_fun(*x[:var_split])

        if compute_rhs:
            RHS = rhs_fun(*x[var_split:])
        else:
            RHS = rhs_fun

        M = (LHS) @ np.conj(RHS).T
        M = M / np.exp(1j * np.angle(M[0, 0]))

        loss = np.linalg.norm(M) - 2 ** (num_qubits / 2)

        return loss

    if random is not None:
        x = np.random.uniform(*random, no_var)
    else:
        x = np.zeros(no_var)

    job = minimize(fun=cost, x0=x, args=(compute_rhs))

    if np.abs(job.fun) < error:
        print("Optimization Sucessful")
    else:
        print("Optimization Unsucessful")

    return job
