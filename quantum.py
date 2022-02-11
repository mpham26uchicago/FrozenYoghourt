from mode import *
from maths import *


class Quantum:

    def random_local_gates(num_qubits):
        if Mode.representation == 'numpy':
            return Maths.tp(*[random_unitary(2) for i in range(num_qubits)])
        else:
            return Maths.tp(*[Gates.U(*symbols(f'theta_{i}, phi_{i}, lambda_{i}')) for i in range(num_qubits)])

    def epsilon(psi, num_qubits=2):
        if Mode.representation == 'numpy':
            return np.abs(psi.T @ Maths.tp(Gates.Y(), no_times=num_qubits) @ psi)[0]
        else:
            return Abs(psi.T @ Maths.tp(Gates.Y(), no_times=num_qubits) @ psi)[0]

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
    
   
    def gamma(U):
        n = int(np.log2(U.shape[0]))

        E = Maths.tp(Gates.Y(), no_times=n)
        return U @ E @ U.T @ E

    def chi(M, variable='x'):
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
            var = symbols(variable)
            variable_mat = Matrix([var ** n for n in range(dim + 1)])
            coef = Matrix(coef).T
            char_poly = simplify(coef @ variable_mat)[0]
            return char_poly

    def double_cosets(g, h):
        E = Gates.Magic()

        g, h = Quantum.to_su(g), Quantum.to_su(h)

        u, v = np.conj(E).T @ g @ E, np.conj(E).T @ h @ E

        D1, P1 = np.linalg.eig(u @ u.T)
        D2, P2 = np.linalg.eig(v @ v.T)

        idx = D1.argsort()[::-1]
        D1 = D1[idx]
        P1 = P1[:, idx]

        idx = D2.argsort()[::-1]
        D2 = D2[idx]
        P2 = P2[:, idx]

        if not np.isclose(D1[0], D2[0]):
            D2 = -D2
            idx = D2.argsort()[::-1]
            D2 = D2[idx]
            P2 = P2[:, idx]

        a = P1.T
        b = P2.T

        c = np.conj(v).T @ b.T @ a @ u

        left_coset = E @ b.T @ a @ np.conj(E).T
        right_coset = E @ c.T @ np.conj(E).T

        return left_coset, right_coset
