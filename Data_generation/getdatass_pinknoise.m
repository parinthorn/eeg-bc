%==========================================================================
%
%	This function generates data from the input state-space system with
%	length = timelength. The noise w(t) is assumed to be a zero mean Gaussian
%	with autocovariance sigma and n(t) is a pink noise.
%
%          z(t+1)  =   A*z(t) + w(t)
%          x(t)    =   C*z(t) + n(t)
%
%	INPUTS
%           system      =   struct obtained from gengcss.m or struct
%                           contain at least the following fields:
%                   sys_ss    : matlab ss object
%                   sigma_w   : noise covariance of AR process
%                                         source
%           timelength	=   length of time series data
%           snr         =   signal per noise ratio between simulated data
%                           and pink noise
%           perturball  =   1 for perturbing all channel with pink noise
%                           0 for perturbing only inactive channel
%	OUTPUTS
%           data        =   time series data generated from system. 
%                           data(:,:,j) represents time series data at
%                           trial j.
%
%========================================================================== 
function [data] = getdatass_pinknoise(system,timelength,snr,perturball)
[nS,~] = size(system.sigma_w); [nA,~] = size(system.sys_ss.A);
m = system.PARAMETER.m;
m_inactive = m-system.PARAMETER.m_active;
ind_active = system.PARAMETER.ind_active;
noise_w = mvnrnd(zeros(1,nS),system.sigma_w,timelength);
%     if system.sigma_eta == 0
%         system.sigma_eta = zeros(nN);
%     end
%     noise_eta = mvnrnd(zeros(1,nN),system.sigma_eta,timelength);
t = 0:timelength-1;
z0 = rand(1,nA);
data_sim = lsim(system.sys_ss,noise_w,t,z0);
if perturball
    pn = dsp.ColoredNoise('pink',timelength,m);
    % disp(size(data_sim)); disp(size(pn()));
    data(:,:) = snr*data_sim'/norm(data_sim,'Fro')+(1-snr)*pn()'/norm(pn(),'Fro');
else
    pn = dsp.ColoredNoise('pink',timelength,m_inactive);
    data(:,:) = snr*data_sim'/norm(data_sim,'Fro');
    ind = setdiff([1:m]',ind_active);
    data(ind,:) = data(ind,:) + (1-snr)*pn()'/norm(pn(),'Fro');
end
end

