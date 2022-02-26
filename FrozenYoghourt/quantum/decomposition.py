from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.gates import *

def step1(k):
    
    k1 = np.copy(k)
    Pauli = [ID(), X(), Y(), Z()]
    shift = k1 // np.pi
    
    # k1eep track1 of whether to add phase to correct for k1Ak1
    phase = np.exp(-1j*(np.pi/2)*(np.sum(shift)%2)) 

    k1 -= np.pi*shift

    pauli_shift = (shift %2).astype(int)
    L = mm(*[tp(Pauli[pauli_val*(index+1)], Pauli[pauli_val*(index+1)]) 
                       for index, pauli_val in enumerate(pauli_shift)])

    # Fix global phase if different phase is -1
    phase *= global_phase(CAN(*k), phase*L@CAN(*k1)) 
    
    return k1, L, phase

def step2(k):
    
    L_list = [np.identity(4)]
    R_list = [np.identity(4)]

    if k[0] < k[1]:
        k[0], k[1] = k[1], k[0]
        L1, R1 = tp(Rz(np.pi/2), no_times = 2), tp(Rz(-np.pi/2), no_times = 2)
        L_list.append(L1), R_list.append(R1)

    if k[1] < k[2]:
        k[1], k[2] = k[2], k[1]
        L2, R2 = tp(Rx(np.pi/2), no_times = 2), tp(Rx(-np.pi/2), no_times = 2)
        L_list.append(L2), R_list.append(R2)
        final = True
    else:
        final = False
    
    if final:
        if k[0] < k[1]:
            k[0], k[1] = k[1], k[0]
            L3, R3 = tp(Rz(np.pi/2), no_times = 2), tp(Rz(-np.pi/2), no_times = 2)
            L_list.append(L3), R_list.append(R3)
    else:
        if k[0] < k[2]:
            k[0], k[2] = k[2], k[0]
            L3, R3 = tp(Ry(np.pi/2), no_times = 2), R@tp(Ry(-np.pi/2), no_times = 2)
            L_list.append(L3), R_list.append(R3)
            
    L = mm(*L_list) # The order goes L1 L2 L3 C3 R3 R2 R1
    R = mm(*R_list[::-1])
    
    return k, L, R

def step3(k):
    L3, R3 = np.eye(4), np.eye(4)
    if k[0] + k[1] > np.pi:
        L3a = -tp(Y(), Y())@tp(X(), X())@tp(Z(), ID())@tp(Rz(np.pi/2), Rz(np.pi/2))
        R3a = tp(Rz(-np.pi/2), Rz(-np.pi/2))@tp(Z(), ID()) # Bring to (pi-ky, pi-kx, kz)

        k, L3b, R3b = step2(np.array([np.pi - k[1], np.pi - k[0], k[2]]) ) # Sorting
        
        L3, R3 = L3a@L3b, R3b@R3a

    return k, L3, R3

def step4(k):
    L4 , R4 = np.eye(4), np.eye(4)
    if np.isclose(k[2], 0) and (k[0] > np.pi/2):
        L4, R4 = -1j*tp(X(), X())@tp(Y(), ID()), tp(Y(), ID())
        k[0] = np.pi - k[0]
        
    return k, L4, R4

def KAK(U, return_canonical = False):
    M = Magic()

    phase1 = global_phase(to_su(U), U)
    V = dagger(M)@to_su(U)@M
    D_squared, P = np.linalg.eig(V.T@V)

    spectrum = np.angle(D_squared) # These steps make sure that D stays unimodular after the square root
    spectrum[3] -= np.sum(spectrum)
    spectrum = (spectrum/2).reshape(4, )

    k_list = Gamma().T@spectrum

    # Step 1, bring each kx, ky, kz to [0, pi) (done) (Check 10_000 values)
    k, L1, phase2 = step1(k_list[1:])

    # Step 2, set kx ≤ ky ≤ kz (done) (Check 20_000 values)
    k, L2, R2 = step2(k)

    # Step 3, Fix kx + ky > pi (done) (Check 20_000 values)
    k, L3, R3 = step3(k)

    # Step 4, fix if kz = 0 and kx > pi/4 
    k, L4, R4 = step4(k)
    
    if np.isclose(np.linalg.det(P), -1): # If the determinant of P is -1, turn it to 1
        P[:, 0] = -P[:, 0]

    K1, K2 = M@P.T@dagger(M), M@V@P@np.diag(np.exp(-1j*spectrum))@dagger(M) # Extract local gates
    
    phase = phase1*phase2 # Collect Phase
    Can = CAN(*k) # Collect Canonical Gate
    L, R = K2@L1@L2@L3@L4, R4@R3@R2@R1@K1 # Collect local gates
    
    if not is_local(L): # Add global phase for factoring
        L*=1j
        R*=-1j
    
    if return_canonical:
        return phase, L, Can, R, k
    else:
        return phase, L, Can, R