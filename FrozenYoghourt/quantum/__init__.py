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

def kron_decomp(U:np.ndarray):
    m = U.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)

    u, sv, vh = np.linalg.svd(m)

    a = np.sqrt(sv[0]) * u[:, 0].reshape(2, 2)
    b = np.sqrt(sv[0]) * vh[0, :].reshape(2, 2)

    return a, b


def gamma(U, unimodular = False):
    n = int(np.log2(U.shape[0]))

    if (Mode.repsentation == 'numpy') and (not np.isclose(np.linalg.det(U), 1)):
        U = to_su(U)
    elif unimodular:
        U = to_su(U)
        
    E = tp(Y(), no_times=n)
    return U @ E @ U.T @ E

def chi(M, variable=None):
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
        return chi(gamma(to_su(U)))
    
    elif Mode.representation == 'numpy':
        charpoly_U = chi(gamma(to_su(U)))
        charpoly_V = chi(gamma(to_su(V)))
        
        if return_charpoly:
            return np.all(np.isclose(charpoly_U, charpoly_V)), charpoly_U, charpoly_V
        else:
            return np.all(np.isclose(charpoly_U, charpoly_V))
        
    else: ## This needs fixing because the determinant function doesn't work
        charpoly_U = chi(gamma(to_su(U)))
        charpoly_V = chi(gamma(to_su(V)))
        
        if return_charpoly:
            return charpoly_U.equals(charpoly_V), charpoly_U, charpoly_V
        else:
            return charpoly_U.equals(charpoly_V)


def canonical_class_vector(x, y, z, return_phase = True):
    
    # Step 1 and 2
    param = np.array([x, y, z])
    shift = param // np.pi
    phase = np.exp(-1j*(np.pi/2)*(np.sum(shift)%2)) # Keep track of whether to add phase to correct for KAK
    k = np.sort(param - np.pi*shift)[::-1]

    # Step 3
    if k[0] + k[1] > np.pi:
        k = np.sort(np.array([np.pi - k[1], np.pi - k[0], k[2]]))[::-1]

    # Step 4
    if np.isclose(k[2], 0) and (k[0] > np.pi/2):
        k[0] = np.pi - k[0]
        k = np.sort(k)[::-1]
        phase *= 1j
    
    if return_phase:
        return k[0], k[1], k[2], phase
    else:
        return k[0], k[1], k[2]

### Add huang_invariant tomorrow