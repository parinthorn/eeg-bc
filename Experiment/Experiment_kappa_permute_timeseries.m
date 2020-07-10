%% Experiment kappa selection from time-series chopping
clear
close all
clc
load('DATA\ss_model_1_eeg_40_10')
load('MODEL\ss_model_1')


%% Extract V, W
nstate = 5;

T = 1;
r=y.r;
M = SYSTEM(T).sys;
L = y.L;
data = y.data(:,end-2000+1:end,T);

% 
[V_1,W_1] =  ...
    extractVW(data,L,[],2*ceil(nstate/r),nstate);

%% Rough estimate of C using L22, L2,0.5
GridSize = 50;

AA = L'*V_1*W_1'; % For alpha_max computation
alpha_max = max(norms(AA,2,2));
Lipschitz = norm(L,2)^2*norm(W_1,2)^2;
alpha = [0,logspace(-5,0,GridSize-1)*alpha_max];
for ii=1:GridSize
if ii==1
    [C_L21(:,:,ii),~] = nmAPG_ss(L,W_1,V_1,alpha(ii),1,0.9/Lipschitz,1);
else
    [C_L21(:,:,ii),~] = nmAPG_ss(L,W_1,V_1,alpha(ii),1,0.9/Lipschitz,1,C_L21(:,:,ii-1));
end
end
for ii=1:GridSize
    fprintf('Estimation progress: %d / %d \n',ii,GridSize)
    [C_Lpq(:,:,ii),~] = nmAPG_ss(L,W_1,V_1,alpha(ii),0.5,0.9/Lipschitz,1,C_L21(:,:,ii));
end

%% Kappa selection

kappa = kappa_selection_timeseries(data,0.5,500,L,alpha,nstate,C_Lpq);
%%
plot(kappa)
















