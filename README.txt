FrozenYoghourt is a collections of useful modules for working with circuit theory. Our library offers a tool for both numerical calculation with numpy that is optimized for speed and symbolic manipulation with sympy that is helpful for studying closed form circuit.

We divide our library into 4 main modules.

1. mode: used for switching between numerical and symbolic representations.
2. gates: contains many useful one and two qubits gates. The symbolic representation are especially conducive to analyzing parameterized families.
3. maths: this module contains many mathematical methods that are generally useful for quantum computing
4. quantum: this module contains specific methods for working with circuit decomposition
5. circuit: we are developping this module to optimize matrix multiplication and tensor product in the context of quantum circuits. This should also allows for easy analysis of quantum state and isometry.