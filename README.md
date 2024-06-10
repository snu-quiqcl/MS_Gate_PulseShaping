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

```CostOptimization.py``` optimizes the process using the gradient descent method. The input values for this process are specified in the ```OptimizationConfiguration.json``` file located in the ```PythonCode```  folder. The parameters in this file are as follows:




