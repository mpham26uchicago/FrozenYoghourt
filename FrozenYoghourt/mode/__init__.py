from FrozenYoghourt import *
    

class Mode:
    """
    A class to help with switching between Numerical mode and Symbolic mode
    
    Attributes
    ----------
    representation: str
        a string to represent the current mode ('numerical' or 'symbolic')
        
    Methods
    -------
    toggle(specific_mode = None, printing = True)
        Switch to an alternate mode, or to a specific inputted mode
    now()
        Print out the current mode
    
    """
    
    representation = 'numerical'
    
    @classmethod
    def toggle(cls, specific_mode = None, printing = True):
        """
        Switch to an alternate mode, or to a specific inputted mode
        
        If the argument specific_mode is not passed in, the default None is used
        which switch to the alternate mode
        
        Parameters
        ----------
        specific_mode: str, optional
            Mode to switch to (default is None)
        printing: bool, optional
            Option to print a message when switching mode (default is True)
        
        """

        if specific_mode is None:
            cls.representation = {'numerical': 'symbolic',
                                  'symbolic': 'numerical'}[cls.representation]
            if printing:
                print(f'\033[1;34m{cls.representation.upper()} \033[1;32mactivated')
        else:
            cls.representation = {'n': 'numerical', 
                                  'numpy':'numerical', 
                                  'numerical':'numerical',
                                  's': 'symbolic', 
                                  'sympy':'symbolic', 'symbolic':'symbolic'}[specific_mode]
            if printing:
                print(f'\033[1;34m{cls.representation.upper()} \033[1;32mactivated')
            
    @classmethod
    def now(cls):
        """
        Print out the current mode
        """
        
        print(f'\033[1;34m{cls.representation.upper()} \033[1;30mis\033[1;30m \033[1;32mon')