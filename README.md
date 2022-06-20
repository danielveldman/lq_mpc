# Model Predictive Control
A simple implementation of Model Predictive Control (MPC) of unconstrained linear dynamics such as (discretized) wave and heat equations with quadratic cost functionals. 
This repository contains code for the application of MPC to linear systems in which the plant model used by the MPC controller can differ from the plant that is being controlled. It also contains codes that show the convergence of MPC in\ $T - \tau $ and\ $ \tau $. This was used to do the numerical simulations in https://arxiv.org/abs/2206.01097. 



<!-- [MPC_wave](figures/MPC_T=41250_tau=1250.jpeg)
%[MPCX_wave](figures/MPCX_T=41250_tau=1250.jpeg)

%![MPCA_wave](figures/MPCA_T=41250_tau=1250.jpeg)
%![MPCAX_wave](figures/MPCAX_T=41250_tau=1250.jpeg) -->

![MPCw_wave](figures/MPCw_T=41250_tau=1250.jpeg)
Figure 1: MPC control obtained by running 
![MPCwX_wave](figures/MPCwX_T=41250_tau=1250.jpeg)
Figure 2: Norm of the state trajectory resulting from the application of the control in Figure 1. 

<p align="center">
![MPCconv_wave](figures/Fig6a.PNG)
Figure 3: Convergence analysis for MPC when varying $T - \tau$
</p>

<!-- ![MPCconvw_wave](figures/Fig6b.PNG)
![MPCconvA_wave](figures/Fig6c.PNG) -->
