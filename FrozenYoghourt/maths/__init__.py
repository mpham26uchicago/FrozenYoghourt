from FrozenYoghourt import *
from FrozenYoghourt.mode import *

close = lambda a, b: np.all(np.isclose(a, b))

def mm(*lst, no_times=1):
    lst *= no_times
    return reduce(lambda a, b: a @ b, lst)

def tp(*lst, no_times=1):
    lst *= no_times
    if Mode.representation == 'numpy':
        return reduce(lambda a, b: np.kron(a, b), lst)
    else:
        return TensorProduct(*lst)

def dagger(*u:Union[Tuple[np.ndarray], Tuple[MutableDenseMatrix]]):
    if Mode.representation == 'numpy':
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
    
def svd(A:np.ndarray):
    if Mode.representation == 'numpy':
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

def to_su(*u:Union[Tuple[np.ndarray], Tuple[MutableDenseMatrix]]):
    if Mode.representation == 'numpy':
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
    if (type(variables) == list) or (type(variables) == tuple):
        substituted_matrix = matrix.subs(list(zip(variables, values)))
    else:
        substituted_matrix = matrix.subs(variables, values)
        
    if to_numpy:
        return np.array(substituted_matrix).astype(complex)
    else:
        return substituted_matrix.evalf()

def optimizing_identities(lhs_fun, rhs_fun, no_var, var_split, random=(0, 2 * np.pi), num_qubits=2, error=1e-7):
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
