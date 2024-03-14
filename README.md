# On-Kernel-Design-for-Regularized-Non-Causal-System-Identification
 
- Here is the formula of the NCSI kernel for the Automatica paper ``On Kernel Design for Regularized  Non-Causal System Identification'', see 
[the link to the preprint version](https://arxiv.org/abs/2307.13999 "On Kernel Design for Regularized  Non-Causal System Identification") This is only partial code on kernel construction functions, and the GP codes are not provided but similar to impulseest() in MATLAB.  The integrated codes will be provided in a developing toolbox.

- main function ncsi_kernel.m
- tc_kernel_SI_bs_test3.m (equivalent to tc_kernel_SI_bs.m) is the file consistent with the paper's notations, but it is inefficient. (tc_kernel_SI_bs.m is more efficient)

- generate_linear_system_randomly.m is the file to generate Monte Carlo experiments of the non-minimum phase systems (D2-D3 in the paper)
