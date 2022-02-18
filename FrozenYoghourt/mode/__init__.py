from FrozenYoghourt import *
    

class Mode:
    
    representation = 'numpy'
    
    @classmethod
    def toggle(cls, specific_mode = None, printing = True):

        if specific_mode is None:
            cls.representation = {'numpy': 'sympy',
                                  'sympy': 'numpy'}[cls.representation]
            if printing:
                print(f'\033[1;34m{cls.representation.upper()} \033[1;32mactivated')
        else:
            cls.representation = {'n': 'numpy', 
                                  'numpy': 'numpy', 
                                  's': 'sympy', 
                                  'sympy':'sympy'}[specific_mode]
            if printing:
                print(f'\033[1;34m{cls.representation.upper()} \033[1;32mactivated')
            
    @classmethod
    def now(cls):
        print(f'\033[1;34m{cls.representation.upper()} \033[1;30mis\033[1;30m \033[1;32mon')