from FrozenYoghourt import *
from FrozenYoghourt.mode import *
from FrozenYoghourt.gates import *
from FrozenYoghourt.maths import *
from FrozenYoghourt.quantum import *

class NumericalCircuit:
    
    def __init__(self, num_qubits):
        self.num_qubits = num_qubits
        self.unitary = np.eye(2**num_qubits)
        
    # class method
    def x(self, qubit: Union[int, list]):
        if isinstance(qubit, int):
            gate_unitary_list = [X() if index == qubit else ID() for index in range(self.num_qubits)]
        elif isinstance(qubit, list):
            qubit.sort()
            gate_unitary_list = [X() if index in qubit else ID() for index in range(self.num_qubits)]
        else:
            raise TypeError("Qubit(s) to apply operation must be an integer or a list")
            
        # Calculate Unitary Matrix of the operation
        gate_unitary = tp(*gate_unitary_list)
        
        # Multiply the following matrix to the circuit unitary
        self.unitary = gate_unitary@self.unitary
        
    def y(self, qubit: Union[int, list]):
        if isinstance(qubit, int):
            gate_unitary_list = [Y() if index == qubit else ID() for index in range(self.num_qubits)]
        elif isinstance(qubit, list):
            qubit.sort()
            gate_unitary_list = [Y() if index in qubit else ID() for index in range(self.num_qubits)]
        else:
            raise TypeError("Qubit(s) to apply operation must be an integer or a list")
            
        # Calculate Unitary Matrix of the operation
        gate_unitary = tp(*gate_unitary_list)
        
        # Multiply the following matrix to the circuit unitary
        self.unitary = gate_unitary@self.unitary
        
    def z(self, qubit: Union[int, list]):
        if isinstance(qubit, int):
            gate_unitary_list = [Z() if index == qubit else ID() for index in range(self.num_qubits)]
        elif isinstance(qubit, list):
            qubit.sort()
            gate_unitary_list = [Z() if index in qubit else ID() for index in range(self.num_qubits)]
        else:
            raise TypeError("Qubit(s) to apply operation must be an integer or a list")
            
        # Calculate Unitary Matrix of the operation
        gate_unitary = tp(*gate_unitary_list)
        
        # Multiply the following matrix to the circuit unitary
        self.unitary = gate_unitary@self.unitary
        
    def h(self, qubit: Union[int, list]):
        if isinstance(qubit, int):
            gate_unitary_list = [H() if index == qubit else ID() for index in range(self.num_qubits)]
        elif isinstance(qubit, list):
            qubit.sort()
            gate_unitary_list = [H() if index in qubit else ID() for index in range(self.num_qubits)]
        else:
            raise TypeError("Qubit(s) to apply operation must be an integer or a list")
            
        # Calculate Unitary Matrix of the operation
        gate_unitary = tp(*gate_unitary_list)
        
        # Multiply the following matrix to the circuit unitary
        self.unitary = gate_unitary@self.unitary
        
    def rx(self, angle: Union[float, list], qubit: Union[int, list]):
        if isinstance(qubit, int):
            
            # Validate the type of "angle"
            assert isinstance(angle, float) or isinstance(angle, int), 'Angle must be a number'
            
            # Create list of gates on each qubit
            gate_unitary_list = [Rx(angle) if qubit_index == qubit else ID() 
                                 for qubit_index in range(self.num_qubits)]
            
        elif isinstance(qubit, list):
            
            # Validate the type of "angle"
            assert isinstance(angle, list), 'Angles must be a list'
            
            gate_unitary_list = [ID()]*self.num_qubits
            for qubit_index, angle_param in enumerate(angle):
                gate_unitary_list[qubit[qubit_index]] = Rx(angle_param)
                
        else:
            raise TypeError("Qubit(s) to apply operation must be an integer or a list")
            
        # Calculate Unitary Matrix of the operation
        gate_unitary = tp(*gate_unitary_list)
        
        # Multiply the following matrix to the circuit unitary
        self.unitary = gate_unitary@self.unitary
        
        def ry(self, angle: Union[float, list], qubit: Union[int, list]):
            if isinstance(qubit, int):

                # Validate the type of "angle"
                assert isinstance(angle, float) or isinstance(angle, int), 'Angle must be a number'

                # Create list of gates on each qubit
                gate_unitary_list = [Ry(angle) if qubit_index == qubit else ID() 
                                     for qubit_index in range(self.num_qubits)]

            elif isinstance(qubit, list):

                # Validate the type of "angle"
                assert isinstance(angle, list), 'Angles must be a list'

                gate_unitary_list = [ID()]*self.num_qubits
                for qubit_index, angle_param in enumerate(angle):
                    gate_unitary_list[qubit[qubit_index]] = Ry(angle_param)

            else:
                raise TypeError("Qubit(s) to apply operation must be an integer or a list")

            # Calculate Unitary Matrix of the operation
            gate_unitary = tp(*gate_unitary_list)

            # Multiply the following matrix to the circuit unitary
            self.unitary = gate_unitary@self.unitary
            
        def rz(self, angle: Union[float, list], qubit: Union[int, list]):
            if isinstance(qubit, int):

                # Validate the type of "angle"
                assert isinstance(angle, float) or isinstance(angle, int), 'Angle must be a number'

                # Create list of gates on each qubit
                gate_unitary_list = [Rz(angle) if qubit_index == qubit else ID() 
                                     for qubit_index in range(self.num_qubits)]

            elif isinstance(qubit, list):

                # Validate the type of "angle"
                assert isinstance(angle, list), 'Angles must be a list'

                gate_unitary_list = [ID()]*self.num_qubits
                for qubit_index, angle_param in enumerate(angle):
                    gate_unitary_list[qubit[qubit_index]] = Rz(angle_param)

            else:
                raise TypeError("Qubit(s) to apply operation must be an integer or a list")

            # Calculate Unitary Matrix of the operation
            gate_unitary = tp(*gate_unitary_list)

            # Multiply the following matrix to the circuit unitary
            self.unitary = gate_unitary@self.unitary
            
        def p(self, angle: Union[float, list], qubit: Union[int, list]):
            if isinstance(qubit, int):

                # Validate the type of "angle"
                assert isinstance(angle, float) or isinstance(angle, int), 'Angle must be a number'

                # Create list of gates on each qubit
                gate_unitary_list = [Phase(angle) if qubit_index == qubit else ID() 
                                     for qubit_index in range(self.num_qubits)]

            elif isinstance(qubit, list):

                # Validate the type of "angle"
                assert isinstance(angle, list), 'Angles must be a list'

                gate_unitary_list = [ID()]*self.num_qubits
                for qubit_index, angle_param in enumerate(angle):
                    gate_unitary_list[qubit[qubit_index]] = Phase(angle_param)

            else:
                raise TypeError("Qubit(s) to apply operation must be an integer or a list")

            # Calculate Unitary Matrix of the operation
            gate_unitary = tp(*gate_unitary_list)

            # Multiply the following matrix to the circuit unitary
            self.unitary = gate_unitary@self.unitary
            
        def cx(self, control = Union[int, list], target = Union[int, list]):
            if isinstance(control, int) and isinstance(target, int):
                gate_unitary_list = [CU(control, target, X(), self.num_qubits)]

            elif isinstance(control, list) and isinstance(target, list):
                # Validate size of control and target list
                assert len(control) == len(target), "List of controls and list of targets don't have the same size"
                gate_unitary_list = [CU(control[i], target[i], X(), self.num_qubits) for i in range(len(control))]
                gate_unitary_list.reverse()

            else:
                raise TypeError("Control and target bit(s) to apply operation must be an integer or a list")

            # Calculate Unitary Matrix of the operation
            gate_unitary = mm(*gate_unitary_list)

            # Multiply the following matrix to the circuit unitary
            self.unitary = gate_unitary@self.unitary