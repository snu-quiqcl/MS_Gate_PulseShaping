# MS_Gate_PulseShaping

## Overview 
This is Ion Trap MS Gate pulse optimization and its simulation program. This program consists of two main components: ```PulseCalculation.py``` and ```CostOptimization.py```.

+ ```PulseCalculation.py```: This script calculates alpha and Theta based on the input pulse. The program stores parameters such as the input Rabi frequency, detuning, and the number of modes.
+ ```CostOptimization.py```: Based on the calculations from ```PulseCalculation.py```, this script optimizes the pulse. The optimization is performed using OpenAI's ADAM algorithm, following the same methodology as described in IonQ's paper.
When running ```CostOptimization.py```, the pulse is optimized according to the specified epochs and learning rate. It then calculates the trajectory of $\alpha_{j}^{k}$ based on the input pulse and computes $\Theta$ over time, displaying the results in a graph.
+ ```normalModeQuadratic.py``` : The program calculates the degenerated mode frequencies and corresponding mode vectors by considering the coulomb interaction based on the number of input modes, and whether they are axial or transverse, as well as the initial mode frequency.

## PulseCalculation
The program calculates the trajectories of $\alpha$, $\Theta$ based on this below unitary operator.

$$ U_{MS}(\tau) = \exp\left(\sum_{j_{1},j_{2}}\sum_{k}(\alpha_{j}^{k}a_{k}^{\dagger}-\alpha_{j}^{k*}a_{k})\sigma_{x}^{j}\right)\exp\left( i\Theta(\tau) \sigma_{x}^{j_{1}}\sigma_{x}^{j_{2}}\right) $$

$$ \Theta(\tau) =\frac{1}{2}\sum_{k}\eta_{k}^{2}b_{j_{1}}^{k}b_{j_{2}}^{k}\int_{0}^{\tau}dt\Omega(t)\sin\theta_{k}(t)\int_{0}^{t}dt'\Omega(t')\cos\theta_{k}(t') -\frac{1}{2}\sum_{k}\eta_{k}^{2}b_{j_{1}}^{k}b_{j_{2}}^{k}\int_{0}^{\tau}dt\Omega(t)\cos\theta_{k}(t)\int_{0}^{t}dt'\Omega(t')\sin\theta_{k}(t') $$

$$ \alpha_{j}^{k}(\tau) = \frac{\eta_{k} b^{k}_ {j}}{2} \int_{0}^{\tau}\Omega(t)\exp(i\theta_{k}(t))dt $$

$$ \theta_{k}(\tau) = \omega_{k}-\int_{0}^{\tau}\delta(t)dt $$

The parameters of above equations are as follows:

+ $\alpha_{j}^{k}$ : Displacement of $j$ th ion,  $k$ th mode.
+ $\eta_{k}$ : Lamb dicke parameter of $k$ th mode. This is calculated by $\sqrt{\frac{\hbar}{2M\omega_{k}}}$
+ $\Omega(t)$ : Rabi frequency of raman transition. Note that this is not rabi frequency of single laser. It is well known that it can be calculated from $\Omega_{r}\Omega_{b}/\Delta$, but it is substituted directly in this program.
+ $\delta(t)$ : Detuning of raman transition
+ $\tau$ : Gate time

## CostOpimization
The ultimate goal of the MS Gate is to converge α to 0 and Θ to π/4 by the end of the gate operation. The program defines these values as the cost function. This definition is chosen because the value obtained by subtracting this cost function from 1 approximates the fidelity.

$$ \varepsilon = \varepsilon_{\alpha} + \varepsilon_{\Theta} $$

$$ \varepsilon_{\alpha} = \sum_{k}(|\alpha_{j_{1}}^{k}|^{2}+|\alpha_{j_{2}}^{k}|^{2}) $$

$$ \varepsilon_{\Theta}=\left(\Theta-\frac{\pi}{4}\right)^{2} $$

```CostOptimization.py``` optimizes the process using the gradient descent method. The input values for this process are specified in the ```OptimizationConfiguration.json``` file located in the ```PythonCode```  folder. Example json files is shown as follows:

```json
{
    "opt_mode" : "AM",
    "initial_value_shuffle" : "False",
    "w_transverse" : "13395751.074906878", 
    "w_axial"  : "3141592.653589793",
    "manual_mode_frequency" : "True",
    "omega" : [
        "12880529.879718151",
        "13395751.074906878"
    ],
    "n_0" : "20",
    "ion_number" : "2",
    "parameter_constraint" : "True",
    "Deltak" : "28339146.473469555",
    "tau" : "200e-6",
    "N" :"1",
    "target_ion_index1" : "0",
    "target_ion_index2" : "1",
    "rescale_factor1" : "1",
    "integral_div" : "1000000",
    "epoch" : "0",
    "derivative" : "1e-9",
    "learning_rate" : "1e4",
    "beta1" : "0.9",
    "beta2" : "0.999",
    "epsilon" : "1e-8",
    "max_Omega" : "628318.5307179586",
    "min_Omega" : "-6283.185307179586",
    "max_delta" : "12566370.614359172",
    "min_delta" : "-6283.185307179586",
    "Omega" : [
        "455370.5040424259"
    ],
    "delta" : [
        "13322148.047112534"
    ]
}
```
+ ```opt_mode```: Specifies the pulse modulation mode. ```AM```, ```FM```, and ```AWG``` mode is possible. Literally, ```AM``` optimizes amplitude within gate time and ```FM``` optimizes freqeuncy within gate time. ```AWG``` optimizes both of them.
+ ```initial_value_shuffle```: Boolean flag indicating whether to make randomness in initial parameter values of $\Omega$ and $\delta$.
+ ```w_transverse```: Transverse motional mode frequency in Hertz.
+ ```w_axial```: Axial motional mode frequency in Hertz.
+ ```manual_mode_frequency```: Boolean flag for using manual mode frequencies. If this flag is ```False```, ```normalModeQuadratic.py``` calculates mode frequency and mode vector based on axial and transverse mode frequencies. On the other hand, if this flag is ```True```, mode frequencies specified in ```omega``` will be used in optimization or population calculation.
+ ```omega```: Array of motional mode frequencies. It is meaningful only when ```manual_mode_frequency = True```
+ ```n_0```: Average thermal mode number. This value does not influence the optimization process. However, they demonstrate the effect of the thermal state in the population calculations performed after the optimization.
+ ```ion_number```: Number of ions in the system.
+ ```parameter_constraint```: Boolean flag indicating whether $\Omega$ or $\delta$ is constrained in the ```max_Omega```, ```min_Omega```, ```max_delta```, and  ```min_delta```. This flag does not guarantee that $\Omega$ and $\delta$ are completely constrained between maximum and minimum parameters. However, by incorporating a ReLU function into the cost function, it introduces a tendency for $\Omega$ and $\delta$ to be optimized within the desired bounds.
+ ```Deltak```: $\Delta k$ value in unit of $[rad/m]$. This values is used when Lamb-Dicke parameter $\eta = \Delta k \sqrt{\hbar/2M\omega}$ is calculated.
+ ```tau```: Pulse duration in seconds.
+ ```N```: Number of divided pulses. Note that this is different from ```integral_div```.
+ ```target_ion_index1```: Index of the first target ion. Note that ion index is in the range of ```0 ~ N-1```
+ ```target_ion_index2```: Index of the second target ion.
+ ```rescale_factor1```: Rescale factor for cost calculation. This parameter does not affect the optimization process but provides convenience in viewing the values.
+ ```integral_div```: Number of divisions for integral calculation.
+ ```epoch```: Number of epochs for optimization(i.e. number of iteration).
+ ```derivative```: Derivative step size for gradient calculation.
+ ```learning_rate```: Learning rate for the ADAM optimizer. This is recommended to set 0.001 in paper. However for fast optimization, it is recommended to set ```1E3 ~ 1E6``` at initial.
+ ```beta1```: Beta1 parameter for the ADAM optimizer, controlling the exponential decay rate for the first moment estimates. This is recommended to set to 0.001, and do not change if there is no specific reason.
+ ```beta2```: Beta2 parameter for the ADAM optimizer, controlling the exponential decay rate for the second moment estimates.This is recommended to set to 0.999, and do not change if there is no specific reason.
+ ```epsilon```: Epsilon parameter for the ADAM optimizer, preventing division by zero.
+ ```max_Omega```: Maximum value for $\Omega$ during optimization.
+ ```min_Omega```: Minimum value for $\Omega$ during optimization.
+ ```max_delta```: Maximum value for $\delta$ during optimization.
+ ```min_delta```: Minimum value for $\delta$ during optimization.
+ ```Omega```: Array of initial values for $\Omega$.
+ ```delta```: Array of initial values for $\delta$.



