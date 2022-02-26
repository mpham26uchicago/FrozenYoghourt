from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.gates import *

def local_ops(num_qubits = 1, no_ops = 1, unimodular = False):
    
    if type(no_ops) == bool:
        unimodular = no_ops; no_ops = 1
    
    if no_ops == 1:
        return tp(*[u2(unimodular, i) for i in range(num_qubits)])
    else:
        return [tp(*[u2(unimodular, i) for i in range(j*num_qubits, (j+1)*num_qubits)]) for j in range(no_ops)]

def epsilon(psi, num_qubits=2):
    if Mode.representation == 'numpy':
        return np.abs(psi.T @ tp(Y(), no_times=num_qubits) @ psi)[0]
    else:
        return Abs(psi.T @ tp(Y(), no_times=num_qubits) @ psi)[0]

def global_phase(A, B):
    D = np.diag(np.conj(A).T @ B)

    if np.isclose(np.min(D), np.max(D)):
        return D[0]
    else:
        print('A and B are not equivalent up to global phase')
        return False

def ymap(U, unimodular = False):
    n = int(np.log2(U.shape[0]))

    if (Mode.representation == 'numpy') and (not np.isclose(np.linalg.det(U), 1)):
        U = to_su(U)
    elif unimodular:
        U = to_su(U)
        
    E = tp(Y(), no_times=n)
    return U @ E @ U.T @ E

def xmap(M, variable=None):
    dim = M.shape[0]

    coef = np.array([1])
    Mk = np.array(M)

    for k in range(1, dim + 1):
        ak = -Mk.trace() / k
        coef = np.append(coef, ak)
        Mk += np.diag(np.repeat(ak, dim))
        Mk = np.dot(M, Mk)

    if Mode.representation == 'numpy':
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

    if V is None:
        return xmap(ymap(to_su(U)))
    
    elif Mode.representation == 'numpy':
        charpoly_U = xmap(ymap(to_su(U)))
        charpoly_V = xmap(ymap(to_su(V)))
        
        if return_charpoly:
            return np.all(np.isclose(charpoly_U, charpoly_V)), charpoly_U, charpoly_V
        else:
            return np.all(np.isclose(charpoly_U, charpoly_V))
        
    else: ## This needs fixing because the determinant function doesn't work
        charpoly_U = xmap(ymap(to_su(U)))
        charpoly_V = xmap(ymap(to_su(V)))
        
        if return_charpoly:
            return charpoly_U.equals(charpoly_V), charpoly_U, charpoly_V
        else:
            return charpoly_U.equals(charpoly_V)


def huang_invariant(U):
    V = dagger(Magic())@to_su(U)@Magic()
    
    if Mode.representation == 'numpy':    
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
    M = Magic()
    V = dagger(M)@to_su(U)@M
    if np.isclose(np.linalg.det(V), 1) and np.all(np.isclose(V@V.T, np.identity(4))):
        return True
    else: 
        return False
    
def kron_decomp(U:np.ndarray):
    
    assert is_local(U), 'Matrix is not local' # Check if gates is local
    
    m = U.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)

    u, sv, vh = np.linalg.svd(m)

    a = np.sqrt(sv[0]) * u[:, 0].reshape(2, 2)
    b = np.sqrt(sv[0]) * vh[0, :].reshape(2, 2)

    return a, b

def KAK(U:np.ndarray):
    
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