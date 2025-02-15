# -*- coding: utf-8 -*-
"""
Created on Mon May 20 17:14:22 2024

@author: alexi
"""

import numpy as np
import json
import matplotlib.pyplot as plt

from scipy.integrate import cumulative_trapezoid

ALIGN_LENGTH=66
IMAGE_DPI = 200

# This is Physical Parameters
class TPP:
    # fixed physcial parameters
    mode_num : int = None
    max_j : int = None
    omega : list[float]  = None# omega_{k} = omega[k] mode frequency
    eta : list[float]  = None# eta_{k} = eta[k] lamb dicke parameter
    b : list[list[float]]  = None# b_{j}^{k} = b[k][j]
    tau : float = None
    tau_array : np.array = None
    N : int = None # Division number of pulse
    target_ion_index1 : int = None
    target_ion_index2 : int = None
    max_Omega : float = None
    max_delta : float = None
    min_Omega : float = None
    min_delta : float = None
    n_0 : int = None
    manual_mode_frequency : bool = None
    
    # variables
    Omega : np.array  = None
    delta : np.array  = None
    
    # Numerical Calculation parameters
    integral_div : int = None # Division for integration

    def __init__(cls):
        pass
    
    @classmethod
    def setGlobalVariables(cls, file_path : str) -> None:
        """
        Discretize the time interval from 0 to t into N points.
        
        t : float
            The upper limit of the time interval.
        N : int
            The number of discretized points in the pulse.
        eta_k : float
            Lamb Dicke parameter of ion
        b_k_j : float
            Mode coefficient of ion
        theta_k : np.array
            Accumulated phase of ion mode
        Omega : np.array
            The Rabi frequency.
        tau : np.array
            The discretized time points. Note that this function returns only
            integral value, not the list
        
        Returns
        -------
        None
        
        """
        hbar : float = 6.626e-34/(2*np.pi) # Plank
        M : float = 170.936323 * (1e-3)/6.02e23 # Mass of 171Yb in [kg]
        
        cls.max_j = 2
        cls.omega : list[float] = [] # mode frequency in [rad/s]
        cls.b : dict[int : dict[str : list[float]]]= dict() # mode vector in [A.U]
        cls.eta : dict[int : float] = dict() # lamb dicke parameter
        
        log : str = ""
        
        with open("NormalModeCalculationReport.json", "r") as file:
            data = json.load(file)
            cls.mode_num = int(data["N"]) # number of mode. Note that this is different from ADAM's N
            
            for i in range(cls.mode_num):
                normal_mode = data[f"{i}_th_mode"]
                cls.omega.append(float(normal_mode["w"]) * 1e6 * np.pi * 2) # convert [MHz] to [rad/s]
                cls.b[i] = {
                    "x" : [float(comp) for comp in normal_mode["b"]],
                    "y" : [0.0 for comp in normal_mode["b"]]
                }
            
            
        with open(file_path, "r") as file:
            data = json.load(file)
            
            if str(data["manual_mode_frequency"]) == "True": 
                cls.manual_mode_frequency = True
            else:
                cls.manual_mode_frequency = False
            
            if cls.manual_mode_frequency:
                omega_list = []
                b_list = dict()
                i = 0
                for freq, mode in data["manual_mode"].items():
                    omega_list.append(float(freq) * 1e6 * np.pi * 2)
                    b_list[i] = {
                        "x" : [float(x) for x in mode["x"]],
                        "y" : [float(y) for y in mode["y"]]
                    }
                    i = i + 1
                cls.omega = omega_list
                cls.b = b_list
                cls.mode_num = len(cls.omega)
            
            cls.n_0 = int(data["n_0"])
            cls.Deltak = {
                "x" : float(data["Deltak"]["x"]),
                "y" : float(data["Deltak"]["y"])
            }
            cls.target_ion_index1 = int(data["target_ion_index1"])
            cls.target_ion_index2 = int(data["target_ion_index2"])
            cls.integral_div = int(data["integral_div"])
            cls.tau = float(data["tau"])
            cls.N = int(data["N"])
            cls.max_Omega = float(data["max_Omega"])
            cls.max_delta = float(data["max_delta"])
            cls.min_Omega = float(data["min_Omega"])
            cls.min_delta = float(data["min_delta"])
            cls.delta : np.array = np.zeros(cls.N)
            cls.Omega : np.array = np.zeros(cls.N)
                
            if cls.integral_div == None or cls.integral_div == 0:
                raise Exception("integral division number is not defined")
            
            interval: int = cls.integral_div // cls.N
            cls.integral_div = interval * cls.N
            log += f">> integral_div is reset to {cls.integral_div} ..."
                
            for i in range(cls.mode_num):
                cls.eta[i] = {
                    "x" : cls.Deltak["x"] * np.sqrt(hbar/(2*M*cls.omega[i])),
                    "y" : cls.Deltak["y"] * np.sqrt(hbar/(2*M*cls.omega[i]))
                }
                
            Omega_list : list[float] = [float(x) for x in data["Omega"]]
            delta_list : list[float] = [float(x) for x in data["delta"]]
            
            if( len(Omega_list) != cls.N ):
                raise Exception("Lenth of Omega is different from N")
                
            if( len(delta_list) != cls.N ):
                raise Exception("Lenth of delta is different from N")

            for i in range(cls.N):
                cls.Omega[i] = Omega_list[i]
                cls.delta[i] = delta_list[i]
                
            cls.tau_array = np.linspace(0, cls.tau, cls.integral_div)
        
        log += f">> maximum Omega \n>> {cls.max_Omega}"
        log += f">> maximum delta \n>> {cls.max_delta}"                    
        log += f">> integral interval \n>> {cls.tau_array[1] -cls.tau_array[0]}\n"
        log += f">> number of mode \n>> {cls.mode_num}\n"
        log += f">> number of ions \n>> {cls.max_j}\n"
        log += f">> mode frequency [rad/s] \n>> {[round(x,3) for x in cls.omega]}\n"
        omega_print : list[float] = [x / (2*np.pi*1e6) for x in cls.omega]
        log += f">> mode frequency [MHz] \n>> {[round(x,3) for x in omega_print]}\n"
        for i in range(cls.mode_num):
            log += (f"+------------------------{i} th eigen vector------------------------+\n")
            log += (f"| b_x : {[round(x,3) for x in cls.b[i]['x']]}".ljust(ALIGN_LENGTH) + '|\n')
            log += (f"| b_y : {[round(x,3) for x in cls.b[i]['y']]}".ljust(ALIGN_LENGTH) + '|\n')
            
        log += f">> Lamb dicke parameter \n>> {[round(x,3) for x in cls.eta]}\n"
        log += f">> Gate Time [s] \n>> {cls.tau}\n"
        log += f">> initial Omega \n>> {[round(x,3) for x in cls.Omega]}\n"
        log += f">> initial delta \n>> {[round(x,3) for x in cls.delta]}\n"    
        print(log)  
    
    
###############################################################################
## Calculate function
###############################################################################
def theta(k: int, delta : np.array) -> np.array:
    """
    Compute theta_k for given omega_k, delta and tau.
    
    k : int
        Mode number of ion
    
    Returns
    -------
    np.array
        The array of theta_k values at each time point in tau.
    """
    theta_values = cumulative_trapezoid(TPP.omega[k] - delta, TPP.tau_array, initial=0)
    return theta_values

def alpha(j : int, k : int, Omega : np.array, delta : np.array) -> complex:
    """
    Compute alpha_j^k(tau) for given eta_k, b_k_j, theta_k, Omega and tau.
    
    k : int
        Mode number of ion
    j : int
        index of ion
    
    Returns
    -------
    complex
        The final value of alpha_j^k at the last time point in tau.
    """
    integrand = Omega * np.exp(1j * theta(k , delta))
    alpha_value = (
        (TPP.eta[k]["x"] * TPP.b[k]["x"][j] / 2) * np.trapz(integrand, TPP.tau_array) + 
        (TPP.eta[k]["y"] * TPP.b[k]["y"][j] / 2) * np.trapz(integrand, TPP.tau_array)
    )
    return alpha_value

def alpha_avg(j : int, k : int, Omega : np.array, delta : np.array) -> complex:
    """
    Compute alpha_j^k(tau) for given eta_k, b_k_j, theta_k, Omega and tau.
    
    k : int
        Mode number of ion
    j : int
        index of ion
    
    Returns
    -------
    complex
        The final value of alpha_j^k at the last time point in tau.
    """
    integrand = Omega * np.exp(1j * theta(k, delta))
    alpha_values = (
        (TPP.eta[k]["x"] * TPP.b[k]["x"][j] / 2) * 
        cumulative_trapezoid(integrand, TPP.tau_array, initial=0) + 
        (TPP.eta[k]["y"] * TPP.b[k]["y"][j] / 2) * 
        cumulative_trapezoid(integrand, TPP.tau_array, initial=0)
    )
    alpha_integral = np.trapz(alpha_values, TPP.tau_array)
    return alpha_integral

def alphaArray(j : int, k : int, Omega : np.array, 
               delta : np.array) -> np.array:
    """
    Compute alpha_j^k(tau) for given eta_k, b_k_j, theta_k, Omega and tau.
    
    k : int
        Mode number of ion
    j : int
        index of ion
    
    Returns
    -------
    complex
        The final value of alpha_j^k at the last time point in tau.
    """
    integrand = Omega * np.exp(1j * theta(k, delta))
    alpha_values = (
        (TPP.eta[k]["x"] * TPP.b[k]["x"][j] / 2) * 
        cumulative_trapezoid(integrand, TPP.tau_array, initial=0) + 
        (TPP.eta[k]["y"] * TPP.b[k]["y"][j] / 2) * 
        cumulative_trapezoid(integrand, TPP.tau_array, initial=0)
    )
    return alpha_values

def totalTheta(Omega : np.array, delta : np.array) -> float:
    """
    Compute Θ(τ) based on the given formula.
    
    Omega : np.array
        The array of Omega values at each time point in tau.
    delta : np.array
        The array of delta values at each time point in tau.
    
    Returns
    -------
    float
        The computed Θ(τ) value.
    """
    if (len(Omega) != len(delta)):
        raise Exception("Length of Omega and delta is different")
    
    result : float = 0
    for k in range(TPP.mode_num):
        theta_k : np.array = theta(k, delta)
        
        integrand_cos = Omega * np.cos(theta(k, delta))
        integrand_sin = Omega * np.sin(theta(k, delta))
        integral_cos = cumulative_trapezoid(integrand_cos, TPP.tau_array, initial=0)
        integral_sin = cumulative_trapezoid(integrand_sin, TPP.tau_array, initial=0)
        integral_cos = integral_cos * Omega * np.sin(theta_k)
        integral_sin = integral_sin * Omega * np.cos(theta_k)
        
        
        eta_k : float = TPP.eta[k]["x"]
        b_k_j1 : float = TPP.b[k]["x"][TPP.target_ion_index1]
        b_k_j2 : float = TPP.b[k]["x"][TPP.target_ion_index2]

        result += (0.5 * eta_k * eta_k * b_k_j1 * b_k_j2 * 
                   np.trapz(integral_cos - integral_sin, TPP.tau_array))
        
        eta_k : float = TPP.eta[k]["y"]
        b_k_j1 : float = TPP.b[k]["y"][TPP.target_ion_index1]
        b_k_j2 : float = TPP.b[k]["y"][TPP.target_ion_index2]
        
        result += (0.5 * eta_k * eta_k * b_k_j1 * b_k_j2 * 
                   np.trapz(integral_cos - integral_sin, TPP.tau_array))
    
    return result

def totalThetaArray(Omega : np.array, delta : np.array) -> np.array:
     """
     Compute Θ(τ) based on the given formula.
     
     Omega : np.array
         The array of Omega values at each time point in tau.
     delta : np.array
         The array of delta values at each time point in tau.
     
     Returns
     -------
     float
         The computed Θ(τ) value.
     """
     if (len(Omega) != len(delta)):
         raise Exception("Length of Omega and delta is different")
     
     result : float = 0
     for k in range(TPP.mode_num):
         theta_k : np.array = theta(k, delta)
         
         integrand_cos = Omega * np.cos(theta(k, delta))
         integrand_sin = Omega * np.sin(theta(k, delta))
         
         integral_cos = cumulative_trapezoid(integrand_cos, TPP.tau_array, initial=0)
         integral_sin = cumulative_trapezoid(integrand_sin, TPP.tau_array, initial=0)
         
         integral_cos = integral_cos * Omega * np.sin(theta_k)
         integral_sin = integral_sin * Omega * np.cos(theta_k)
         
         eta_k : float = TPP.eta[k]["x"]
         b_k_j1 : float = TPP.b[k]["x"][TPP.target_ion_index1]
         b_k_j2 : float = TPP.b[k]["x"][TPP.target_ion_index2]

         result = (
             result + (
                 0.5 * eta_k**2 * b_k_j1 * b_k_j2 * 
                 cumulative_trapezoid(
                     integral_cos - integral_sin, 
                     TPP.tau_array, initial=0
                 )
             )
         )
         
         eta_k : float = TPP.eta[k]["y"]
         b_k_j1 : float = TPP.b[k]["y"][TPP.target_ion_index1]
         b_k_j2 : float = TPP.b[k]["y"][TPP.target_ion_index2]

         result = (
             result + (
                 0.5 * eta_k**2 * b_k_j1 * b_k_j2 * 
                 cumulative_trapezoid(
                     integral_cos - integral_sin, 
                     TPP.tau_array, initial=0
                 )
             )
         )
     
     return result           
        
def extendVariable(original_variable : np.array) -> np.array:
    extended_variable = np.zeros(TPP.integral_div)
    pulse_interval = TPP.integral_div // TPP.N
    
    for i in range(TPP.N):
        start_index = i * pulse_interval
        end_index = start_index + pulse_interval
        if end_index > TPP.integral_div:
            end_index = TPP.integral_div
        extended_variable[start_index:end_index] = original_variable[i]
    return extended_variable
    
def plotComplex(name : str , complex_array : np.array, 
                gate_index : int = 0) -> None:
    real_parts = np.real(complex_array)
    imag_parts = np.imag(complex_array)
    color_len : int = 31
    color_interval = len(complex_array) // color_len
    
    cmap = plt.get_cmap('rainbow')
    colors = [cmap(i/(color_len-1)) for i in range(color_len)]
    
    plt.figure(figsize=(8, 6),dpi=IMAGE_DPI)
    
    for i in range(color_len):
        start_index : int = i * color_interval
        end_index : int = start_index + color_interval
        if end_index > len(complex_array):
            end_index = len(complex_array)
        plt.plot(real_parts[start_index:end_index], 
                 imag_parts[start_index:end_index], 
                 color=colors[i], linewidth=1)
            
    plt.scatter(real_parts[-1], imag_parts[-1], 
                color='black', marker='x', s=100, label='Last Point')
    # plt.scatter(real_parts[gate_index], imag_parts[gate_index], 
    #             color='red', marker='x', s=100, label='Gate Point')
    plt.xlabel('Real Part')
    plt.ylabel('Imaginary Part')
    plt.title(name)
    plt.grid(True)
    plt.legend()
    plt.savefig(f"figures/{name}.png", dpi=IMAGE_DPI)
    plt.show()
    
def plotTimeEvolution(
        name : str , variable : np.array, gate_index : int = 0) -> None:
    plt.figure(dpi=IMAGE_DPI)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(TPP.tau_array*1000000, variable)
    plt.xlabel('time [us]')
    plt.ylabel('Variable')
    if name == "Theta":
        if variable[gate_index] > 0:
            ax.axhline(y=np.pi / 4, color='r', linestyle='--', label='Θ = π/4')
        if variable[gate_index] < 0:
            ax.axhline(y=-np.pi / 4, color='r', linestyle='--', label='Θ = -π/4')
            
    if name != "Theta":
        plt.ylim(0,1)
    
    # ax.axvline(x=TPP.tau_array[gate_index]*1000000, color = 'r',  linestyle='--', 
    #            label=f'gate time = {TPP.tau_array[gate_index]*1000000}')
    
    ax.legend()
    
    plt.title(name)
    plt.savefig(f"figures/{name}.png", dpi=IMAGE_DPI)
    plt.show()

def plotSimpleTimeEvolution(name : str , variable : np.array) -> None:
    plt.figure(dpi=IMAGE_DPI)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(TPP.tau_array*1000000, variable)
    plt.xlabel('time [us]')
    plt.ylabel('Variable')
    
    plt.title(name)
    plt.savefig(f"figures/{name}.png", dpi=IMAGE_DPI)
    plt.show()

def plotVariable(name : str , variable : np.array) -> None:
    plt.figure(dpi=IMAGE_DPI)
    positions = np.arange(len(variable))
    plt.bar(positions, variable/(1e6*2*np.pi))
    plt.xlabel('index')
    plt.ylabel('Variable [MHz]')
    plt.title(name)
    plt.savefig(f"figures/{name}.png", dpi=IMAGE_DPI)
    plt.show()


