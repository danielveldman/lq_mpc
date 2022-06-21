# Model Predictive Control
A simple implementation of Model Predictive Control (MPC) of unconstrained linear dynamics such as (discretized) wave and heat equations with quadratic cost functionals. 
This was used to do the numerical simulations in https://arxiv.org/abs/2206.01097. 

Main features:
* Comparison the MPC control to the inifnite horizon optimal control and the MPC limit (as defined in https://arxiv.org/abs/2206.01097).
* Imperfections in the plant model can be included. 
* Numerical validation of convergence rates for $(T-\tau, \tau) \rightarrow (\infty, 0)$.

![MPCw_wave](figures/MPCw_T=41250_tau=1250.jpeg)
Figure 1: MPC control obtained by running MPC_wave_w.m
![MPCwX_wave](figures/MPCwX_T=41250_tau=1250.jpeg)
Figure 2: Norm of the state trajectory resulting from the application of the control in Figure 1. 

![test](figures/fig6a.jpg)
Figure 3: Convergence analysis for MPC when varying $T - \tau$

<!-- ![MPCconvw_wave](figures/Fig6b.PNG)
![MPCconvA_wave](figures/Fig6c.PNG) -->
