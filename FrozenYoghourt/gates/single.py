from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.double import *
from FrozenYoghourt.multi import *

def ID():
    if Mode.representation == 'numerical':
        return np.array([[1, 0], [0, 1]])
    else:
        return Matrix([[1, 0], [0, 1]])

def X():
    if Mode.representation == 'numerical':
        return np.array([[0, 1], [1, 0]])
    else:
        return Matrix([[0, 1], [1, 0]])

def Y():
    if Mode.representation == 'numerical':
        return np.array([[0, -1j], [1j, 0]])
    else:
        return Matrix([[0, -I], [I, 0]])

def Z():
    if Mode.representation == 'numerical':
        return np.array([[1, 0], [0, -1]])
    else:
        return Matrix([[1, 0], [0, -1]])
    
def H():
    if Mode.representation == 'numerical':
        return 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]])
    else:
        return 1 / sqrt(2) * Matrix([[1, 1], [1, -1]])

def Rx(theta:float):
    if Mode.representation == 'numerical':
        return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                         [-1j * np.sin(theta / 2), np.cos(theta / 2)]])
    else:
        return Matrix([[cos(theta / 2), -I * sin(theta / 2)],
                       [-I * sin(theta / 2), cos(theta / 2)]])

def Ry(theta:float):
    if Mode.representation == 'numerical':
        return np.array([[np.cos(theta / 2), -np.sin(theta / 2)],
                         [np.sin(theta / 2), np.cos(theta / 2)]])
    else:
        return Matrix([[cos(theta / 2), -sin(theta / 2)],
                       [sin(theta / 2), cos(theta / 2)]])

def Rz(theta:float):
    if Mode.representation == 'numerical':
        return np.array([[np.exp(-1j * theta / 2), 0],
                         [0, np.exp(1j * theta / 2)]])
    else:
        return Matrix([[exp(-I * theta / 2), 0],
                       [0, exp(I * theta / 2)]])

def Phase(theta:float):
    if Mode.representation == 'numerical':
        return np.array([[1, 0], [0, np.exp(1j*theta)]])
    else:
        return Matrix([[1, 0], [0, exp(I*theta)]])
    
def u2(angles:list = None, unimodular = False, index = ''):
    
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