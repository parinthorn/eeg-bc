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
%           factorC     =   a factor to multiply to C in order to correct
%                           the result from multiplying SNR dividing by
%                           norm
%
%	Written by ANAWAT NARTKULPAT, 2020
%
% Licenses
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%==========================================================================
function [signal_noise,pinknoise_cov,factorC] = getdatass_pinknoise(system,timelength,SNR,is_pink,perturball)

if nargin < 5
    perturball = 0;
end

[nS,~] = size(system.PARAMETER.sigma_w); [nA,~] = size(system.source_model0.A);
m = system.PARAMETER.m;
ind_active = system.PARAMETER.ind_active;
noise_w = mvnrnd(zeros(1,nS),system.PARAMETER.sigma_w,timelength+100);
t = [0:timelength+99]*system.source_model0.Ts;
z_init = rand(1,nA);
signal_sim = lsim(system.source_model0,noise_w,t,z_init);
signal_sim = signal_sim(101:end,:)';

pn = dsp.ColoredNoise('pink',timelength,m); pinknoise = pn()'; % dimension = #source x timepoint
if ~is_pink
    pinknoise = randn(m,timelength);
end
if ~perturball
    pinknoise(ind_active,:) = 0;
end

factorC = SNR/norm(signal_sim,'Fro');
signal_noise = SNR*signal_sim/norm(signal_sim,'Fro');
pinknoise_cov = zeros(m);

if SNR < 1 && norm(pinknoise,'Fro') > 0
    pn_normalized = (1-SNR)*pinknoise/norm(pinknoise,'Fro');
    signal_noise = signal_noise + pn_normalized;
    % calculate variance for each channel of pink noise assuming diagnal
    % some elements in the diagonal are zeros
    pinknoise_cov = diag(diag(cov(pn_normalized')));
end

end