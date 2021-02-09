% This is MATLAB script file to test "EEG GC framework" if it has bias in GC
% estimation when the ground-truth system is time-varying, contradicting to
% the assumption of the model which is time-invariant.
% 
% To do so, we generate time-varying VAR models, time-invariant diagonal
% filter (MA part), and constant C,L
% 
% z(t) = A z(t) + w(t) % state
% x(t) = C z(t) + n(t) % source
% y(t) = L x(t) + e(t) % EEG
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

% The path may need an adjustment upon how the main functions are stored in the
% local folder

addpath('../pvo_subspace/subfun');
addpath('./estimation_functions');
tostartupcvx

%%   Generate VAR model parameters
n = 3; p = 2 ; % AR parameters 
m = 5; % inactive sources = m-n
r = 10; % number of electrodes
Ts = 1; sigma = 0.1 ; % variance of noise w(t)
T = 5000; % number of time points
NUMRUNS = 500; % number of repetition to get histogram of F
K = 10; % Number of switching VAR models

figfilename = 'timevarying_eegGC_decreasing_K2_v1';
% figfilename = 'timevarying_eegGC_increasing_K1_v1';

%   Generate an invertible stable diagonal filter and convert to ss model
filter_tf = gen_diagfilter(0,1,n);
sys_filter = ss(filter_tf,'minimal');

% indFnz = [1 2; 3 1];
indFnz = [1 3;2 3;2 1]; % example 1
% indFnz = [2 1;2 3;3 1]; % example 2

if K==1
    Anz = 0.2; Aii = 0.5;
else
%     Anz = linspace(0,0.2,K); Aii = 0.5; % A is increasing
    Anz = linspace(0.2,0,K); Aii = 0.5; % A is decreasing
end

L = randn(r,m); % lead-field matrix don't vary
F0 = zeros(n,n,K); eegsys = cell(K,1); eegsys2sim = cell(K,1);
for ii=1:K
	A = zeros(n,n,p);
    
    A(:,:,1) = diag(Aii*ones(n,1)); A(:,:,2) = diag(Aii*0.5*ones(n,1));
        for k = 1:size(indFnz,1)-1
            A(indFnz(k,1),indFnz(k,2),1) = Anz(ii);
            A(indFnz(k,1),indFnz(k,2),2) = Anz(ii)*0.5;
        end
        A(indFnz(end,1),indFnz(end,2),1) = max(Anz); % this index, A is not varying
        A(indFnz(end,1),indFnz(end,2),2) = max(Anz)*0.5;

%     [A,indnz] = gen_sparseVAR_simple(n,p,0.8);
    stateA = [A(:,:,1) A(:,:,2); eye(n) zeros(n)];
    if any(abs(eig(stateA)) >=1 )
            error('AR is unstable');
    end

    %   Convert VAR to state-space
    sys_var = ss(stateA,[eye(n); zeros(n)],[eye(n) zeros(n)],0,Ts);
    
    %   Obtain a series connected state-space model of sys_var and sys_filter
    sys = series_ss2ss(sys_var,sys_filter);
    
    %   Calculated the GC matrix of sys
    [nstate,~] = size(sys.A);
    F0(:,:,ii) = calgcss(sys.A,sys.C,sys.B*sigma*sys.B',sys.D*sigma*sys.D',zeros(nstate,n)); 
    
    % Generate system with EEG observation equation
    % padd zero to C (inactive source). Initially, there are 'n' sources
    nstate = size(sys.A,1);

    C = zeros(m,nstate); % Assume that C is still time-invariant
    C(1:n,:) = sys.C; % put the last (m-n) rows to zeros--> indicating inactive source
    
    eegsys{ii} = ss(sys.A,sys.B,L*C,0,Ts); % eegsys{ii} only differ by sys.A (sizes do not change)
    eegsys2sim{ii} = ss(sys.A,sys.B,eye(nstate),0,Ts); % when used with LSIM to get z(t) as output
end

disp('GC matrix is');
F0


%% Generate Time series and estimate GC on non-stationary data

Fhat = zeros(m,m,NUMRUNS); % contain NUMRUNS of estimated GC

parfor ii = 1:NUMRUNS
    
    z0 = rand(nstate,1); noisey_vec = sqrt(sigma)*randn(T,r);
    
    if K==1
        y = lsim(eegsys{1},sqrt(sigma)*randn(T,n),(0:(T-1))',z0) + noisey_vec;    % y is T x r 

    else
        timepointk = repmat(ceil(T/K),1,K-1); timepointk(K) = T-ceil(T/K)*(K-1); % number of time points in each chunk
    
        y = [];
        for k=1:K
            noisew_vec = sqrt(sigma)*randn(timepointk(k),n);
            noisey_vec = sqrt(sigma)*randn(timepointk(k),r);
        
            ztmp = lsim(eegsys2sim{k},noisew_vec,(0:(timepointk(k)-1))',z0)'; % after transpose, size is nstate x T

            z0 = ztmp(:,end-1); % pass the last value of this chunk as the inital value for another segment
            y = [y ; (eegsys{k}.C*ztmp)'+noisey_vec];
        end
    end
[estsys,C_out] = subid_eeg_L21(y,L,10,-1);  % i = 10 and let the program choose 'nstate'
C_ind = C_out.ind_chosen_L21.bic;
Ahat = estsys.A; % This is nstate x nstate
Chat = C_out.C_L21_CLS(:,:,C_ind); % This is m x n 
Ehat = C_out.E_L21_CLS(:,:,C_ind); % covariance of measurement noise, e(t)
What = C_out.W_L21_CLS(:,:,C_ind); % covariance of state noise, w(t)

Nhat = noisecovest(Ehat,L,'homo'); % covariance of source noise, n(t)

F = calgcss(Ahat,Chat,What,Nhat);
Fhat(:,:,ii) = F;

% disp('Estimated GC is');
% F
% F0
end

%% Plot graphs
% n^2-n = 6 ; use (2,3) 
load jsscolor

figure(1); fig = tiledlayout(2,3); fig.TileSpacing = 'compact'; fig.Padding = 'compact';
timeindex = [(0:ceil(T/K):(K-1)*ceil(T/K))' ; T] ; % transition time

tileindex = 0;
for ii=1:n
    for jj=1:n
        if ii ~= jj
            tileindex = tileindex + 1;
            pooledF = squeeze(Fhat(ii,jj,:)) ; % add the last value again
            
            nexttile;
            histogram(pooledF,50,'normalization','probability','FaceColor',jsscolor.bl); 
            
            if (tileindex == 1 || tileindex == 4 )
                ylabel('Probability');
            end
            hold on;
            
            
            yyaxis right 
            pooledF0 = [F0(ii,jj,1) ; squeeze(F0(ii,jj,:)) ]; % varying F0_{ij}
            stairs(pooledF0,timeindex,'linewidth',2);
            xlabel(['$$F_{',num2str(ii),num2str(jj),'}$$'],'interpreter','latex');
             if (tileindex == 3 || tileindex == 6 )
                ylabel('Time');
            end           
            set(gca,'XLim',[-0.05 max(F0(:))*1.1],'YTick',timeindex,'fontsize',20);
        end
    end
end

set(gcf,'WindowState','fullscreen')

% eval(['save ',figfilename, ' T K F0 Fhat']);
% print(figfilename,'-painters','-depsc','-r300');