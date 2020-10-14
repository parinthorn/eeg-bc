# eeg-bc
This repository is for EEG-bc project version control.

The project is to develop MATLAB codes for estimating Granger Causality in source space from EEG time series. Our group contains three members:

Parinthorn Manomaisaowapak and Anawat Nartkulpat and Jitkomut Songsiri
Department of Electrical Engineering, Faculty of Engineering
Chulalongkorn University, Bangkok, Thailand 
e-mail: parinthorn@gmail.com, anawat.nart0@gmail.com, jitkomut.s@chula.ac.th
 
The detail of mathematical formulation is described in our manuscript: 

MS ID#: BIORXIV/2020/329276

MS TITLE: Granger Causality Inference in EEG Source Connectivity Analysis: A State-Space Approach

The folder contains
data_generation: generate state-space models whose parameter correspond to have some sparse Granger causality pattern.

source_selection: estimation of state-space model with sparse prior on the rows of C (source output matrix).

gc_computation: calculate estimated Granger causality matrix.

pvo_subspace: original necessary files for subspace identification by Peter Van Overschee.

input_data: examples of EEG time series, head model, all required inputs to run experiment.

experiment: codes for each experiment explained in the paper.

Dependencies: Our codes rely on several sources in the following. 
1) Based codes for stochastic subspace identification by Peter Van OVerschee
https://homes.esat.kuleuven.be/~smc/sysid/software/
2) Model generation by S. Haufe et.al. This includes the calculation of head model and EEG signal generation. Our implementation extends the codes from Haufe. 
http://bbci.de/supplementary/EEGconnectivity/BBCB.html
3) CVX: Our program require CVX installed in MATLAB (to solve the noise covariance estimation problem, which is a convex program)
Available at http://cvxr.com/cvx/download/
4) Some buit-in MATLAB commands from control or DSP toolboxes, e.g, solving RICCATI equation, generating pinknoise.

