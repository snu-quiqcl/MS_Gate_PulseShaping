# -*- coding: utf-8 -*-
"""
Created on Sun May 26 14:19:41 2024

@author: alexi
"""

import numpy as np
import json
import PulseCalculation as PC
import normalModeQuadratic

class MOP:
    C_matrix : np.array
    def __init__ (self):
        pass
    
    @classmethod
    def setGlobalVariables(cls, file_path : str) -> None:
        with open(file_path, "r") as file:
            data = json.load(file)
            
def optimizeParameters() -> None:
    """
    Optimize the parameters (Omega and delta) to minimize the cost function.
    """
    pass

def doCalculation():
    normalModeQuadratic.doCalculation()
    
    # Optimize the parameters
    optimizeParameters()
    
if __name__ == "__main__":
    normalModeQuadratic.doCalculation()
    MOP.setGlobalVariables(file_path = "")
    doCalculation()