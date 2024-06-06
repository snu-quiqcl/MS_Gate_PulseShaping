# -*- coding: utf-8 -*-
"""
Created on Sun May 26 14:28:27 2024

@author: alexi
"""

import numpy as np
import PulseCalculation as PC
import normalModeQuadratic
import json
import matplotlib.pyplot as plt
from typing import Tuple, Literal

ALIGN_LENGTH = 66

class ADAM(PC.TPP):
    derivative : float
    learning_rate : float
    beta1 : float
    beta2 : float
    epsilon : float
    batch : int
    rescale_factor1 : float 
    cost_values : list[float] = []
    log : str = ""
    opt_mode : Literal["AM", "FM", "AWG"]
    initial_value_shuffle : bool = None
    parameter_constraint : bool = None
    
    def __init__(self):
        pass
    
    @classmethod
    def setGlobalVariables(cls, file_path : str) -> None:
        with open(file_path, "r") as file:
            data = json.load(file)
            ion_number : int = int(data["ion_number"])
            w_axial : float = float(data["w_axial"])
            w_transverse : float = float(data["w_transverse"])
            
            normalModeQuadratic.doCalculation(
                ion_number = ion_number,
                w_axial = w_axial,
                w_transverse = w_transverse
            )
            
            cls.derivative = float(data["derivative"])
            cls.learning_rate = float(data["learning_rate"])
            cls.beta1 = float(data["beta1"])
            cls.beta2 = float(data["beta2"])
            cls.epsilon = float(data["epsilon"])
            cls.batch = int(data["batch"])
            cls.rescale_factor1 = float(data["rescale_factor1"])
            cls.opt_mode = str(data["opt_mode"])
            
            if str(data["initial_value_shuffle"]) == "True": 
                cls.initial_value_shuffle : bool = True
            else:
                cls.initial_value_shuffle : bool = False
            if str(data["parameter_constraint"]) == "True":
                cls.parameter_constraint : bool = True
            else:
                cls.parameter_constraint : bool = False
            if not cls.opt_mode in ["AM", "FM", "AWG"]:
                raise Exception(f">> Not unknown optimization "
                                "mode {cls.opt_mode}")
            
        PC.TPP.setGlobalVariables(file_path)
        
        if cls.initial_value_shuffle == True:
            if cls.opt_mode in ["AM", "AWG"]:
                PC.TPP.Omega = np.random.uniform(
                    min(PC.TPP.Omega) * 0.9, max(PC.TPP.Omega) * 1.1, PC.TPP.N)
            if cls.opt_mode in ["FM", "AWG"]:
                PC.TPP.delta = np.random.uniform(
                    min(PC.TPP.delta) * 0.9, max(PC.TPP.delta) * 1.1, PC.TPP.N)
    
    @classmethod
    def getGradient(
            cls, 
            Omega_original : np.array, 
            delta_original : np.array
        ) -> Tuple[np.array, np.array]:
        """
        Compute the numerical gradient of the cost function with respect to
        Omega_original and delta_original using finite differences.
    
        Parameters:
        Omega_original (np.array): The original Omega parameters.
        delta_original (np.array): The original delta parameters.
    
        Returns:
        Tuple[np.ndarray, np.ndarray]: The gradient with respect to 
        Omega_original and delta_original.
        """
        Omega_original = np.array(Omega_original)
        delta_original = np.array(delta_original)
    
        grad_Omega = np.zeros_like(Omega_original)
        grad_delta = np.zeros_like(delta_original)
    
        # Calculate the baseline cost
        baseline_cost = cls.costFunction(Omega_original, delta_original)
    
        # Compute the gradient with respect to Omega_original
        if cls.opt_mode in ["AM", "AWG"]:
            for i in range(len(Omega_original)):
                Omega_perturbed = Omega_original.copy()
                derivative = Omega_perturbed[i] * cls.derivative
                Omega_perturbed[i] += derivative
                perturbed_cost = cls.costFunction(Omega_perturbed, delta_original)
                grad_Omega[i] = (perturbed_cost - baseline_cost) / derivative
    
        # Compute the gradient with respect to delta_original
        if cls.opt_mode in ["FM", "AWG"]:
            for i in range(len(delta_original)):
                delta_perturbed = delta_original.copy()
                derivative = delta_perturbed[i] * cls.derivative
                delta_perturbed[i] += derivative
                perturbed_cost = cls.costFunction(Omega_original, delta_perturbed)
                grad_delta[i] = (perturbed_cost - baseline_cost) / derivative
            
        return grad_Omega, grad_delta
    
    @classmethod
    def addLog(cls, log:str):
        cls.log += log
        
    @classmethod
    def saveLog(cls):
        with (open("CostOptimizationReport.log", "w") as file):
            file.write(cls.log)
        
    @classmethod
    def optimizeParameters(cls) -> None:
        """
        Optimize the parameters (Omega and delta) to minimize the cost function
        using ADAM algorithm.
        """
        cost : float = cls.costFunction(Omega_original = PC.TPP.Omega,
                                     delta_original = PC.TPP.delta)
        Omega_grad : np.array = np.zeros(len(PC.TPP.Omega))
        delta_grad : np.array = np.zeros(len(PC.TPP.delta))
        m_Omega_t : np.array = np.zeros(len(PC.TPP.Omega))
        m_delta_t : np.array = np.zeros(len(PC.TPP.delta))
        v_Omega_t : np.array = np.zeros(len(PC.TPP.Omega))
        v_delta_t : np.array = np.zeros(len(PC.TPP.delta))
        log:str = "\n"
        
        for i in range(cls.batch):
            log += "\n"
            log +=(
                "##################################################################\n"
                f"##                        {i} th gradient                       ##\n"
                "##################################################################\n")
            
            log += (f"+------------------------{i} th omega value------------------------+\n")
            log += ('[')
            j : int = 0
            for x in PC.TPP.Omega:
                log +=(f"    {str(x)}")
                j += 1
                if j != len(PC.TPP.Omega):
                    log += ",\n"
            log += "]\n"
            log += (f"+------------------------{i} th delta value------------------------+\n")
            log += ('[')
            j : int = 0
            for x in PC.TPP.delta:
                log +=(f"    {str(x)}")
                j += 1
                if j != len(PC.TPP.Omega):
                    log += ",\n"
            log += "]\n"
            
            (Omega_grad, delta_grad) = cls.getGradient(
                Omega_original = PC.TPP.Omega, delta_original = PC.TPP.delta)
            
            m_Omega_t = cls.beta1 * m_Omega_t + (1-cls.beta1) * Omega_grad
            m_delta_t = cls.beta1 * m_delta_t + (1-cls.beta1) * delta_grad
            
            v_Omega_t = cls.beta2 * v_Omega_t + (1-cls.beta2) * Omega_grad * Omega_grad
            v_delta_t = cls.beta2 * v_delta_t + (1-cls.beta2) * delta_grad * delta_grad
            
            m_Omega_t_bc = m_Omega_t/(1-(cls.beta1 ** (i+1)))
            m_delta_t_bc = m_delta_t/(1-(cls.beta1 ** (i+1)))
            
            v_Omega_t_bc = v_Omega_t/(1-(cls.beta2 ** (i+1)))
            v_delta_t_bc = v_delta_t/(1-(cls.beta2 ** (i+1)))
            
            PC.TPP.Omega = (PC.TPP.Omega - cls.learning_rate * m_Omega_t_bc/(np.sqrt(v_Omega_t_bc) + cls.epsilon))
            PC.TPP.delta = (PC.TPP.delta - cls.learning_rate * m_delta_t_bc/(np.sqrt(v_delta_t_bc) + cls.epsilon))
            cost = cls.costFunction(Omega_original = PC.TPP.Omega,
                                delta_original = PC.TPP.delta)
            cls.cost_values.append(cost)
            
            log += (f">> v_Omega_t : {v_Omega_t}\n")
            log += (f">> m_Omega_t : {m_Omega_t}\n")
            log += (f">> v_delta_t : {v_delta_t}\n")
            log += (f">> m_delta_t : {m_delta_t}\n")
            
            log += (f">> Omega grad : {cls.learning_rate * m_Omega_t_bc/(np.sqrt(v_Omega_t_bc) + cls.epsilon)}\n")
            log += (f">> delta grad : {cls.learning_rate * m_delta_t_bc/(np.sqrt(v_delta_t_bc) + cls.epsilon)}\n")
            
            log += (f">> cost : {cost}\n")
            cls.addLog(log)
            print(log)
        
        log += "\n"
        log +=(
            "##################################################################\n"
            "##                      Final Omgea Report                      ##\n"
            "##################################################################\n")
        log += ('[\n')
        j : int = 0
        for x in PC.TPP.Omega:
            log +=(f"    \"{str(x)}\"")
            j += 1
            if j != len(PC.TPP.Omega):
                log += ",\n"
        log += "\n]\n"
        log += "\n"
        log +=(
            "##################################################################\n"
            "##                      Final delta Report                      ##\n"
            "##################################################################\n")
        log += ('[\n')
        j : int = 0
        for x in PC.TPP.delta:
            log +=(f"    \"{str(x)}\"")
            j += 1
            if j != len(PC.TPP.Omega):
                log += ",\n"
        log += "\n]\n"
        log += "\n"
        
        cls.addLog(log)
        print(log)
        cls.saveLog()
            
    @classmethod
    def costFunction(
            cls,
            Omega_original : np.array, 
            delta_original : np.array
        ) -> float:
        """
        Cost function to be minimized.
        This function calculates |alpha(j, k)|^2 for all j and k and sums them up.
        """
        Omega = np.zeros(PC.TPP.integral_div)
        delta = np.zeros(PC.TPP.integral_div)

        pulse_interval = PC.TPP.integral_div // PC.TPP.N
        
        for i in range(PC.TPP.N):
            start_index = i * pulse_interval
            end_index = start_index + pulse_interval
            if end_index > PC.TPP.integral_div:
                end_index = PC.TPP.integral_div
            Omega[start_index:end_index] = Omega_original[i]
            delta[start_index:end_index] = delta_original[i]
            
        cost = 0
        for k in range(PC.TPP.max_k):
            alpha_value = PC.alpha(PC.TPP.target_ion_index1, k, Omega, delta)
            cost += np.abs(alpha_value)**2
            
        for k in range(PC.TPP.max_k):
            alpha_value = PC.alpha(PC.TPP.target_ion_index2, k, Omega, delta)
            cost += np.abs(alpha_value)**2
        
        cost = cost * cls.rescale_factor1
        
        theta_value = PC.totalTheta(Omega, delta)
                
        cost += np.abs(theta_value+np.pi/4)**2  * cls.rescale_factor1
        
        if cls.parameter_constraint == True:
            for x in Omega_original:
                cost += max(-(x - PC.TPP.min_Omega),0) * cls.rescale_factor1 / max(Omega_original)
            for x in Omega_original:
                cost += max((x - PC.TPP.max_Omega),0) * cls.rescale_factor1 / max(Omega_original)
        
        # I found positive theta value through this code...
        # cost -= sum(np.log(Omega)) * cls.rescale_factor1
        
        return cost

def doCalculation() -> None:
    ADAM.optimizeParameters()
    plotResult()

    
def plotCostFunction() -> None:
    plt.figure(dpi=300)
    plt.plot(ADAM.cost_values)
    plt.xlabel('Iteration')
    plt.ylabel('Cost function value')
    plt.title('Cost function value during optimization')
    plt.savefig(f"figures/Cost.png", dpi=300)
    plt.show()
    
def plotResult() -> None:
    plotCostFunction()
    
    for i in range(PC.TPP.max_k):
        PC.plotComplex(
            f"{PC.TPP.target_ion_index1} ion, {i} mode", 
            PC.alphaArray(
                PC.TPP.target_ion_index1,
                i,
                PC.extendVariable(PC.TPP.Omega), 
                PC.extendVariable(PC.TPP.delta)
            )
        )
        PC.plotComplex(
            f"{PC.TPP.target_ion_index2} ion, {i} mode", 
            PC.alphaArray(
                PC.TPP.target_ion_index2,
                i,
                PC.extendVariable(PC.TPP.Omega), 
                PC.extendVariable(PC.TPP.delta)
            )
        )
    PC.plotTimeEvolution("Theta",PC.totalThetaArray(
        PC.extendVariable(PC.TPP.Omega), 
        PC.extendVariable(PC.TPP.delta))
    )
    
    PC.plotVariable("Omega",PC.TPP.Omega)
    PC.plotVariable("delta",PC.TPP.delta)
    
if __name__ == "__main__":
    ADAM.setGlobalVariables(file_path = "OptimizationConfiguration.json")
    doCalculation()
