from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.single import *
from FrozenYoghourt.multi import *

def CX():
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
    
    # Control-unitary gate
    
    if Mode.representation == 'numerical':
        return np.block([[ID(), np.zeros((2, 2))], 
          [np.zeros((2, 2)), U]])
    
    else:
        return BlockDiagMatrix(ID(), U).as_explicit()

def Gamma():
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