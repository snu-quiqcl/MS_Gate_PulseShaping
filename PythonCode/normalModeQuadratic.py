import numpy as np
import json
import copy as cp
import inspect

from typing import Literal
from scipy.optimize import fsolve

ALIGN_LENGTH = 66

class ImmutableMeta(type):
    def __setattr__(cls, key, value):
        # Check the call stack to see if this method is being called from within the class
        stack = inspect.stack()
        # Stack index 1 is the caller of __setattr__, index 2 is the caller of the caller, and so on
        if len(stack) > 1 and stack[1].function in cls.__dict__:
            super().__setattr__(key, value)
        else:
            if hasattr(cls, key):
                raise Exception("You cannot change constant value")
            else:
                raise Exception(f"You do not have constant {key}")
            super().__setattr__(key, value)

class NMQ(metaclass=ImmutableMeta):
    m_ion : float = None
    l_scale : float = None
    w_transverse : float = 2*np.pi*2*1E6
    w_axial : float = 2*np.pi*0.5*1E6
    z_array : list[float] = None
    u_array : list[list[float]]= [[],
               [],
               [-0.62996, 0.62996],
               [-1.0772, 0, 1.0772],
               [-1.4368, -0.45348, 0.45348, 1.4368],
               [-1.7429, -0.8221, 0, 0.8221, 1.7429],
               [-2.0123, -1.1361, -0.36992, 0.36992, 1.1361, 2.0123],
               [-2.2545, -1.4129, -0.68694, 0, 0.68694, 1.4129, 2.2545],
               [-2.4758, -1.6621, -0.96701, -0.31802, 0.31802, 0.96701, 1.6621, 3.4758],
               [-2.6803, -1.8897, -1.2195, -0.59958, 0, 0.59958, 1.2195, 1.8897, 2.6803],
               [-2.8708, -2.10003, -1.4504, -0.85378, -0.2821, 0.2821, 0.85378, 1.4504, 2.10003, 2.8708]]

    c : float = 2.99792458*1E10     # cm/s
    e : float = (1.602176634*1e-19)*c/10   # SI to cgs conversion (esu)
    
    log : str = ''
    data : dict = {}
    
    def __init__(self):
        pass
    
    @classmethod
    def setPhysicalParameter(cls, w_transverse : float, w_axial : float) -> None:
        cls.w_transverse = w_transverse
        cls.w_axial = w_axial
        cls.m_ion = cls.returnMass()
        cls.l_scale = cls.returnLength()
        cls.z_array = cls.returnZArray()
        cls.log = ''
        
        log = (">> length scale = {} um\n".format(round(cls.l_scale,3)))
        cls.writeLog(log)
        print(log)

    @classmethod
    def returnMass(cls) -> float:
        N_ion = 171
        m_p = 1.67262192*1e-27                      # kg
        cls.m_ion = N_ion*m_p    
        cls.m_ion = cls.m_ion*1E3                           # g                   
        return cls.m_ion
    
    @classmethod
    def returnLength(cls) -> float:
        cls.l_scale = (cls.e**2/(cls.m_ion*(cls.w_axial**2)))**(1/3)   # cm = 10 mm
        cls.l_scale = cls.l_scale*1E4                      # um
        return cls.l_scale   
    
    @classmethod
    def returnZArray(cls) -> np.array:
        cls.z_array = cp.deepcopy(cls.u_array)
        for N in range(2, len(cls.z_array)):
            cls.z_array[N] = [i*cls.l_scale for i in cls.z_array[N]]
        return cls.z_array
    
    @classmethod
    def writeLog(cls, log : str) -> None:
        cls.log += log
        return
    
    @classmethod
    def saveLog(cls) -> None:
        with (open("NormalModeCalculationReport.log", "w") as file):
            file.write(cls.log)
            
    @classmethod
    def writeData(cls, data : dict) -> None:
        cls.data.update(data)
        
    @classmethod
    def saveData(cls) -> None:
        file_path = "NormalModeCalculationReport.json"
        with open(file_path, "w") as file:
            json.dump(cls.data, file, indent=4)



def returnUArray(N : int) -> tuple[int, np.array]:
    init_guess = []
    init_ratio = 10     # important
    
    for i in range(N):
        init_guess.append(1 + i/init_ratio)
        
    root = fsolve(eqFunc, init_guess, args = N)
    
    for r in range(len(root)):
        root[r] = round(root[r], 4)
        
    return N, np.array(root)

def eqFunc(x : list[float], N : int) -> list[float]:
    eq_set = list()
    for i in range(N):
        sum = 0
        for j in range(N):
            if j < i:
                sum -= 1/(abs(x[i] - x[j]))**2
            if j > i:
                sum += 1/(abs(x[i] - x[j]))**2
        y = x[i] + sum
        eq_set.append(y)
    return eq_set

def constructK(N : int, mode : Literal['transverse', 'axial']) -> np.array:
    K = np.zeros((N,N), dtype=float)
    me_factor = ((NMQ.e**2)/NMQ.m_ion)*1E12
    
    if mode == 'transverse':
        c_factor = 1
        w_factor = NMQ.w_transverse
        
    elif mode == 'axial':
        c_factor = -2
        w_factor = NMQ.w_axial
        
    for n in range(N):
        for m in range(N):
            if m == n:
                sum_ = 0
                for i in range(N):
                    if i != n:
                        sum_ += 1/(abs(NMQ.z_array[N][n] - NMQ.z_array[N][i]))**3
                K[n][m] = w_factor**2 - c_factor*me_factor*sum_
            else:
                single = 1/(abs(NMQ.z_array[N][n] - NMQ.z_array[N][m]))**3
                K[n][m] = c_factor*me_factor*single

    return K
    
# Main Calculation
def diagK(N : int, mode : Literal['transverse', 'axial'], 
           K : np.array) -> None:
    
    log = ''
    data : dict = {}
    w, v = np.linalg.eig((K))
    ev_list : list = list()
    neg_count : int = 0
    
    log += (
        "##################################################################\n"
        "##                     Eigen Vector Report                      ##\n"
        "##################################################################\n")
    
    for i in range(N):
        if w[i] > 0:
            ev_list.append((round((np.sqrt(w[i])/(2*np.pi))*1e-6, 3), v[:, i]))
        else:
            neg_count += 1
    
    ev_list.sort(key=getEigenvalue)
    
    for i in range(len(ev_list)):
        log += (f"+------------------------{i} th eigen value------------------------+\n")
        log += (f'| w_k : {ev_list[i][0]} [MHz]'.ljust(ALIGN_LENGTH-1) + '|\n')
        
        ev_list_str : list(str) = []
        ev_list_full_str : list(str) = []
        for j in range(len(ev_list[i][1])):
            ev_list_str.append(round(ev_list[i][1][j],3))
            ev_list_full_str.append(str(ev_list[i][1][j]))
        log += (f'| eigenvector: {ev_list_str}'.ljust(ALIGN_LENGTH-1) + '|\n')
        data[f"{i}_th_mode"] = {
            "w" : f"{ev_list[i][0]}",
            "b" : ev_list_full_str
        }
        NMQ.writeData(data)
    
    log += "+"
    for i in range(ALIGN_LENGTH-2):
        log += "-"
    log += "+\n"
            
    
    log += (">> Negative eigenvalue counts: {}\n".format(neg_count))
    print(log)
    NMQ.writeLog(log)

def getEigenvalue(pair):
    return pair[0]

def printInfo(N : int, mode : Literal['transverse', 'axial']) -> None:
    spacing_list = []
    log : str = ""
    
    for i in range(N-1):
        spacing = abs(NMQ.z_array[N][i] - NMQ.z_array[N][i + 1])
        spacing_list.append(round(spacing, 3))
        
    log +=("##################################################################\n")
    log +=("##                     Ion Placement Report                     ##\n")
    log +=("##################################################################\n")
    log +=(">> ")
    
    for i in range(N-1):
        log +=(f"i{i+1}<-")
        log +=(str(spacing_list[i]))
        log +=("->")
    log +=(f"i{N}\n\n")

    ta_ratio = NMQ.w_transverse/NMQ.w_axial
    zz_threshold = 0.77*N/(np.log(N))**(1/2)
    
    log +=("##################################################################\n")
    log +=("##                       Zig-Zag Report                         ##\n")
    log +=("##################################################################\n")
    
    log +=('>> transverse-axial ratio = {}\n'
          '>> zig-zag threshold = {}\n'.format(
        round(ta_ratio, 3), round(zz_threshold, 3)))
    
    if ta_ratio > zz_threshold:
        log +=('>> zig-zag configuration does not occur\n')
    else:
        log +=('>> zig-zag configuration occurs\n')

    if mode == 'transverse':
        w = NMQ.w_transverse/(2*np.pi*1E6)
    elif mode == 'axial':
        w = NMQ.w_axial/(2*np.pi*1E6)
    log +=('>> {} normal modes with single {} frequency {} MHz\n'.format(
        mode, mode, round(w, 3)))
    
    print(log)
    NMQ.writeLog(log)
    

def doCalculation(
        w_transverse : float = 2*np.pi*2*1E6, 
        w_axial : float = 2*np.pi*0.5*1E6,
        ion_number : int = 5
    ) -> None:
    log : str = ''
    data : dict = {}
    N : int = ion_number
    mode : str = 'transverse'
    NMQ.setPhysicalParameter(w_transverse = w_transverse, 
                             w_axial = w_axial)
    data["M"] = str(NMQ.returnMass())
    data["N"] = str(N)
    data["mode"] = mode
    NMQ.writeData(data)

    N, root = returnUArray(N)
    
    log +=("##################################################################\n")
    log +=("##                   Array Calculation Report                   ##\n")
    log +=("##################################################################\n")
    
    u_array_str = np.array2string(root, precision=3, separator=', ')
    z_array_str = np.array2string(NMQ.l_scale * root, precision=3, separator=', ')
    
    log +=(('>> ion number N = {}\n'
           '>> axial freq = {} [MHz]\n'
           '>> u_array = {} [A.U]\n'
           '>> z_array = {} [um]\n').format(
               N, NMQ.w_axial/(2*np.pi*1E6), u_array_str, z_array_str))
               
    print(log)
    NMQ.writeLog(log)
    
    printInfo(N, mode)
    K = constructK(N, mode)
    diagK(N, mode, K)
    
    print(">> Normal mode calculation ended without error...")
    NMQ.saveLog()
    NMQ.saveData()
    print(">> Calculated result is saved...")
    
if __name__ == "__main__":
    doCalculation(w_transverse = 2*np.pi*2*1E6, 
                  w_axial = 2*np.pi*0.5*1E6,
                  ion_number = 2)