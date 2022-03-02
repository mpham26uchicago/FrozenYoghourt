Change Log
----------

0.0.1 (10/02/2022)
------------------

1. Wrote README file
2. Add P gates method to gates.py
3. Add CU method to gates.py
4. Add view method to mode.py
5. Add log.txt for keeping track of changes

0.0.3 (11/02/2022)
------------------

1. Move to_su to maths.py
2. Move kron_decomp to quantum.py
3. Change default variable in the chi method to "x"

0.0.7 (12/02/2022)
------------------

1. Import gates to quantum

0.0.8 (12/02/2022)
------------------

1. Fix Quantum.double_cosets by importing the correct packages
2. Add default_import method to allow for faster import prompt
3. Change random_local_gates to random_local_ops and allow for creating more operation at the same time.
4. Allow for doing to_su on list of matrices.

0.0.10.1 (15/02/2022)
---------------------

1. Add view method to visualize numerical matrices
2. Add CAN method
3. Add Gamma gates

0.0.11 (17/02/2022)
-------------------

1. Delete Class from files so now g.CAN will just be CAN. Although the Mode class is kept.
2. Fix default_import to match the change in 1.

0.0.12 (17/02/2022)
-------------------

1. Allow for custom custom mode toggle
2. Add color to toggle / now

0.0.12.2 (17/02/2022)
---------------------

1. Add a dagger function that can operate on multiple matrices input
2. Fix the to_su function so that it can operate on multiple matrices input
3. Delete double_cosets to be replaced with the KAK
4. Add printing parameter to toggle to allow for not printing results

0.0.13 (18/02/2022)
-------------------

1. Change P gate to Phase
2. Change chi to include coefficients return of symbolic matrices.
3. Change random_local_ops to local_ops
4. Replace U function with u2 function which now allows for special unitary gates
5. Fix local_ops so that it nows includes a special unitary option
6. Add evaluate numerical returns for to_su
7. Add canonical_class_vector method

0.0.14 (19/02/2022)
-------------------

1. Write documentation for mode
2. Allow single variable input for fast_substitution
3. Add angles parameters to u2 function
4. Create two new modules circuit and visualization
5. Fix the default_import to reflect the above changes.
6. Add scatter to visualization
7. Fix gamma, implicitly convert matrix to unimodular if numpy and optionally if sympy
8. Add close to math.py to compare matrices and get boolean values

0.0.15 (22/02/2022)
----------------------

1. Add huang_invariant to quantum
2. Add KAK to quantum!
3. Change Id method to ID to avoid collision
4. Add a pauli method to gates to compute tensor product of pauli matrices
5. Change no_times argument in tp and mm to mult (for multiplicity)
6. Change view in __init__ to allow for displaying sympy matrices with rounding
7. Change gamma method to ymap and chi method to xmap
8. Fix is_local to return boolean value
9. Create a new decomposition.py file
10. Delete canonical_class_vector for step1-4 in decomposition

0.0.16 (06/03/2022)
---------------------

1. Change 'numpy' to 'numerical' and 'sympy' and 'symbolic'
2. Fix CU from Id to ID
3. Add a NumericalCircuit method to circuit. Includes 'x, y, z, h, rx, ry, rz, cx' method.