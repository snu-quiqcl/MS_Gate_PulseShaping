# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 19:53:09 2024

@author: alexi
"""
import numpy as np
import json

class Hamiltonian:
    n_ion : int = None
    n_motional : int = None
    n_state : int = None
    
    omega_state : list[float] = None
    omega_motional : list[float] = None
    eta : list[float] = None
    b : list[list[float]] = None
    
    Omega_r : float = None
    Omega_b : float = None
    delta_r : float = None
    delta_b : float = None
    phi_r : float = None
    phi_b : float = None
    
    H : np.array = None
    time_offset : float = None
    
    def __init__(self):
        pass
    
    @classmethod
    def setGlobalVariables(cls, file_path : str):
        with open(file_path, "r") as file:
            data = json.load(file)
            
            cls.n_ion = int(data["n_ion"])
            cls.n_motional = int(data["n_motional"])
            cls.n_state = int(data["n_state"])
            cls.omega_motional = [float(x) for x in data["omega_motional"]]
            cls.omega_state = [float(x) for x in data["omega_state"]]
            cls.eta = [float(x) for x in data["eta"]]
            
            if cls.n_motional != len(cls.omega_motional):
                raise Exception("n_motional and number of omega_motional "
                                "is different")
            if cls.n_state != len(cls.omega_state):
                raise Exception("n_state and number of omega_state "
                                "is different")
            
            cls.H = np.zeros((cls.n_motional * cls.n_state, cls.n_motional * cls.n_state))
    
    @classmethod
    def updateHamiltonian(
            cls, 
            Omega_r : float, 
            Omega_b : float, 
            delta_r : float,
            delta_b : float,
            phi_r : float,
            phi_b : float,
            time_offset : float
        ) -> np.array:
        cls.Omega_r = Omega_r
        cls.Omega_b = Omega_b
        cls.delta_r = delta_r
        cls.delta_b = delta_b
        cls.phi_r = phi_r
        cls.phi_b = phi_b
        cls.time_offset = time_offset