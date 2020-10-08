%% GC_ESTIMATION_SCRIPT
% 
% This script runs a demo of all processes in GC estimation
% 
% Input: EEG time series
% Output: Estimated GC matrix
% 

clear
clc

% path to load data example
addpath('./input_data/eeg_simulated_data/');
% if you install CVX (required to solve noise estimation problem) on a Linux
% machine, you may startup CVX first. Change your MATLAB path if needed
% run /usr/local/MATLAB/R2019a/cvx/cvx_startup.m

load('eegdata_act_10_deep_0_1') % load one example of EEG time series
load('model_act_10_deep_0_1') % load the corresponding ground-truth model

%% Estimation

% this is L (lead-field) used in estimation when its number of columns
% (m_tilde) may or may not be equal to the actual number of true sources

% IND_CHAN is the index of selected EEG channels

% select to use 108 electrodes
L = model.('L0'); ind_chan = 1:size(L,1);
% select to use 61 electrodes
% L = model.('L10'); ind_chan = model.('L10_ind');
% select to use 32 electrodes
% L = model.('L20_32'); ind_chan = model.('L20_ind_32');
% select to use 19 electrodes
% L = model.('L20'); ind_chan = model.('L20_ind');


r = size(L,1); 

data = eegdata.EEG_data(ind_chan,:,1); % EEG time series  of size num_channels x time points
% estimate state-space parameters
[sys,C_out] = subid_eeg_Lpq(data,L,10,-1);  % i = 10 and let the program choose 'n'
                    
C_ind = C_out.ind_chosen_Lpq.bic; % index of penalty parameter chosen by BIC
A = sys.A;
C = C_out.C_Lpq_CLS(:,:,C_ind);
E = C_out.E_Lpq_CLS(:,:,C_ind);
W = C_out.W_Lpq_CLS(:,:,C_ind); 

% noise covariance estimation
[N,V,obj] = noisecovest(E,L,'homo');

% calculate GC matrix
[F,~] = calgcss(A,C,W,N,[]);

%% Extract the ground-truth parameters (not available in practice)
Ftrue = model.F0;
Ctrue = model.source_model0.C;

%% Plot estimated C
figure(1)
subplot(121)
imagesc(log10((abs(Ctrue))));colormap(flipud(hot));
title('ground-truth $$C$$','Interpreter','latex','fontsize',18)
subplot(122)
imagesc(log10(abs(C)));colormap(flipud(hot));
title('$$\hat{C}_{\mathrm{bic}}$$','Interpreter','latex','fontsize',18)

%% Plot graph of GC estimation
figure(2);
subplot(121);
imagesc(model.F0)
title('true GC matrix')
axis('square')
colormap(flipud(hot))

subplot(122);
imagesc(F)
title('estimated GC matrix')
axis('square')
colormap(flipud(hot))

                    