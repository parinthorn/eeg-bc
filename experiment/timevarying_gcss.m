
% This is MATLAB script file to test "state-space GC estimation" 
% when the ground-truth system is time-varying, contradicting to
% the assumption of the model which is time-invariant.
% 
% To do so, we generate time-varying VAR models, time-invariant diagonal
% filter (MA part), and constant C
% 
% z(t) = A z(t) + w(t) % state
% x(t) = C z(t) + n(t) % output of state-space
% 
% The results of estimated GC (Fhat matrix) at zero location is robust to
% time-varying components. Fhat at locations with time-varying component
% has bias toward low values. Fhat does not seem to alter with the rate of
% time-varying coefficients.
%
% Copyright to Jitkomut Songsiri
% Reference: Granger Causality Inference in EEG Source Connectivity Analysis: A State-Space Approach
% 
% Parinthorn Manomaisaowapak, Anawat Nartkulpat, Jitkomut Songsiri
% 
% https://www.biorxiv.org/content/10.1101/2020.10.07.329276v3

clear all; close all; clc;

addpath('../data_generation/');
addpath('../gc_computation');
addpath('../pvo_subspace/subfun');
addpath('../source_selection');
addpath('../others');

%   Generate VAR model parameters
n = 2; Ts = 1; sigma = 0.1 ; % noise of w(t)
T = 1000; % number of time points
K = 2; % number of transition time during (0,T) that AR model switches 
NUMRUNS = 500; % number of runs to repeat
Fhat = zeros(n,n,NUMRUNS);

figfilename = 'timevaryingGC_decreasing_08to01K20';
% figfilename = 'timevaryingGC_decreasing_08to01K20';
% figfilename = 'timevaryingGC_increasing_01to08K20';

%   Generate an invertible stable diagonal filter and convert to ss model
filter_tf = gen_diagfilter(0,1,n);
sys_filter = ss(filter_tf,'minimal');

a12 = linspace(0,0.8,K); % a12 varies in the interval
A = zeros(n,n,K); sys = cell(K,1); F0 = zeros(n,n,K); % true FC matrix
for k=1:K
    A(:,:,k) = [0.9 a12(k);0 0.8]; % a set of VAR coefficients
    
    %   Convert VAR to state-space
    sys_var = ss(A(:,:,k),eye(n),eye(n),0,Ts);
    
    %   Obtain a series connected state-space model of sys_var and sys_filter
    sys{k} = series_ss2ss(sys_var,sys_filter);
    
    %   Calculated the GC matrix of sys
    [nstate,~] = size(sys{k}.A);
    F0(:,:,k) = calgcss(sys{k}.A,sys{k}.C,sys{k}.B*sigma*sys{k}.B',sys{k}.D*sigma*sys{k}.D',zeros(nstate,n)); 
end

disp('GC matrix is');
F0

%% Generate Time series and estimate GC on non-stationary data

for ii = 1:NUMRUNS
z0 = rand(nstate,1); timepointk = repmat(ceil(T/K),1,K-1); timepointk(K) = T-ceil(T/K)*(K-1);
x = [];
for k=1:K
    noise_vec = sqrt(sigma)*randn(timepointk(k),n);
    sys2sim = ss(sys{k}.A,sys{k}.B,eye(nstate),zeros(nstate,n),Ts); % change C to I so that we can sim the state
    
    ztmp = lsim(sys2sim,noise_vec,(0:(timepointk(k)-1))',z0)';

    z0 = ztmp(:,end-1);
    x = [x sys{k}.C*ztmp];
end


% Estimation of GC
[Ahat,~,Chat,~,Khat,Rhat] = subid(x,[],2*n,n,[],'CVA',1);

F = calgcss(Ahat,Chat,Khat*Rhat*Khat',Rhat,Khat*Rhat);
Fhat(:,:,ii) = F;

disp('Estimated GC is');
F
end

%% Plot graph

figure(1); fig = tiledlayout(1,2); fig.TileSpacing = 'compact'; fig.Padding = 'compact';

nexttile;
timeindex = [(0:ceil(T/K):(K-1)*ceil(T/K))' ; T] ; % transition time
% plot(timeindex,squeeze(F0(1,2,:)),'o-',timeindex,F(1,2)*ones(K,1),'s-',timeindex,F(2,1)*ones(K,1),'x-'); 

F02plot = [squeeze(F0(1,2,:)) ; F0(1,2,end)]; % add the last value again for plotting with 'stairs'
Ftmp = squeeze(Fhat(1,2,:)); % collect all samples
Fhat2plot = repmat(Ftmp,1,K+1);
fig1 = plot(timeindex,Fhat2plot); hold on;
stairs(timeindex,F02plot,'linewidth',2); 

xlabel('Time'); ylabel('GC from $$x_2$$ to $$x_1$$','interpreter','latex');
legend('$$F_{12}$$','$$\hat{F}_{12}$$','interpreter','latex','location','northeast');
set(gca,'FontSize',28);
fig1 = get(gca); 

nexttile;
h = histogram(Ftmp,50,'normalization','probability','orientation','horizontal');
title('Histogram');
xlabel('Probability'); ylabel('$$\hat{F}_{12}$$','interpreter','latex');
set(gca,'YLim',fig1.YLim,'FontSize',28);

set(gcf,'WindowState','fullscreen')


% eval(['save ',figfilename, ' T K F0 Fhat']);
% print(figfilename,'-painters','-depsc','-r300');
