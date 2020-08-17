function generate_eegdata_from_model(nbatch,sys,foldername,plotting,postfix,sa,is_pink)
%
%   Jitkomut Songsiri, Parinthorn Manomaisaowapak, Anawat Nartkulpat, 2020
%
%   This function generate an EEG model and EEG data from the
%   generated model in 'nbatch' batchs.
%
%   INPUT:
%           nbatch      =   number of data batch for the generated model
%           foldername  =   name of the folder for saving data
%           plotting    =   1 for plotting everything
%                           0 for no plotting (not save any figure)
%           postfix     =   add string postfix to outputfile name
%
%   OUTPUT:
%       No direct output from this function. But the model parameters and
%       EEG data will be saved in the 'data' folder with the directory and
%       files described as below;
%
%   data --> 'foldername' 
%                |
%                |-> eegdata_postfix.mat   :   contain eeg and source data matrices
%                |-> model_postfix.mat     :   contain model parameters
%                |
%                --> figure
%                
%
%==========================================================================
%
%   ORIGINAL README BY HAUFE
%
% Stefan Haufe, 2014, 2015
% stefan.haufe@tu-berlin.de
%
% Arne Ewald, Sept. 2015, mail@aewald.net
%
% If you use this code for a publication, please ask Stefan Haufe for the
% correct reference to cite.

% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.

close all; 
% addpath tools
% addpath tools/matlab_bgl
% addpath tools/export_fig
% addpath TestingCode
% addpath pvo_subspace/subfun

% load head model and some miscallaneous data
load('data/miscdata')

% create directory to store data in
mkdir(['data/' foldername])

if nargin <4
    postfix = '';
end
if nargin<7
    is_pink = 1;
end

%%
% initialize random number generator
rng('default');
rng('shuffle');
rng_seed = rng;

% truth.dataset = foldername;
% sys = model.sys;
truth = sys.truth;

% parameter of state-space model of sources REVISE


PARAMETER = sys.PARAMETER;
n_source = PARAMETER.m; %(m in the document)
% n_source_active = PARAMETER.m_active;
% number of electrodes
EEG_M = length(sa.EEG_clab_electrodes);

% number of cluster of sources (must be less than n_source)
% n_source_cluster = PARAMETER.n_source_cluster;


% spatial standard deviation of the sources (along cortical manifold) in mm
% truth.sigma_range = [10 10];

% SNR range to sample from
% truth.snr_range = [0.1 0.9];

% sample snr
%   truth.snr = truth.snr_range(1) + diff(truth.snr_range)*rand(1); % Can be changed later (for randomly choosing snr value)
% truth.snr = 0.95;

% SNR for sensor
% truth.snr_sensor = 0.95;

%number of biological noise source
% n_noise_sources = 500;

% sampling frequency
fs = 100;
% PARAMETER.fs = fs;

% length of recording in sec
% truth.len = 3*60; 

% resulting number of samples
N = fs*truth.len;

% band of interest, in which interaction takes place
% truth.bandpass = [8 13];

% highpass filter coefficients %% == Why 0.1 hertz? can it be changed?
[b_high a_high] = butter(3, 1/fs*2, 'high');


% loop over datasets to be generated

% create folder for dataset idata
%   mkdir(['data/' foldername '/EEG/dataset_' num2str(idata)])
%   mkdir(['data/' foldername '/truth/dataset_' num2str(idata)]

%% Model generation

% generate source time series restricted to band of interest with the
% parameter of the model

% sys = gengcss_eeg(PARAMETER,truth.bandpass,fs,0);  % system of source dynamics
imagesc(sys.F0)
axis('square')
pause(0.1)

%% Data generation

% pre-allocate source data and EEG data matrices
sources = zeros(n_source,N,nbatch);
EEG_data = zeros(EEG_M,N,nbatch);
% sys.source_model = cell(nbatch,1);
% sys.PARAMETER.pinknoise_cov = zeros(PARAMETER.m,PARAMETER.m,nbatch);
% sys.F = zeros(PARAMETER.m,PARAMETER.m,nbatch);
% factorC = zeros(nbatch,1);


%=================== start batch generation loop ==========================

for ibatch = 1:nbatch
%% Sources data generation
fprintf('Generating batch number: %d \n',ibatch)

% generate sources data from the source system and corrupt inactive sources with pink noise
[sources(:,:,ibatch),~,~] = getdatass_pinknoise(sys,N,truth.snr,is_pink);

% recalculate F using the sample covariance of the generated pink noise
% sys.F(:,:,ibatch) = calgcss(sys.source_model0.A,sys.source_model0.C,sys.PARAMETER.sigma_w,sys.PARAMETER.pinknoise_cov(:,:,ibatch));

% sys.source_model{ibatch} = sys.source_model0; 
% sys.source_model{ibatch}.C = sys.source_model{ibatch}.C*factorC(ibatch);

%% time series generation

% multiple the corrupted sources data with the lead field matrix
EEG_brain_signal_noise = truth.EEG_field_pat*sources(:,:,ibatch);  % the pink noise is already corrupted in the source

% add sensor noise (from eeg probes) as white noise
EEG_sensor_noise = randn(EEG_M, N);
EEG_sensor_noise = EEG_sensor_noise ./ norm(EEG_sensor_noise, 'fro');

% overall noise is dominated by biological noise
EEG_data(:,:,ibatch) = truth.snr_sensor*EEG_brain_signal_noise + (1-truth.snr_sensor)*EEG_sensor_noise;

% apply high-pass filter
EEG_data(:,:,ibatch) = filtfilt(b_high, a_high, EEG_data(:,:,ibatch)')';

%% generate pseudo-baseline EEG/MEG without sources 
% everything as above except that no signal is added at all

% JSS baseline equation
pn = dsp.ColoredNoise(1,N,n_source); purenoise = pn()'; % pursenoise is n_source x time points
purenoise(sys.PARAMETER.ind_active,:)  = 0;   % assume that there is no pink noise at active source channels
if norm(purenoise,'fro') > 0
EEG_brain_noise = truth.EEG_field_pat*purenoise;
EEG_brain_noise = EEG_brain_noise/norm(EEG_brain_noise,'fro');   
else
EEG_brain_noise = truth.EEG_field_pat*purenoise;
end

% white sensor noise
EEG_sensor_noise = randn(EEG_M, N);
EEG_sensor_noise = EEG_sensor_noise ./ norm(EEG_sensor_noise, 'fro');

EEG_baseline_data = truth.snr_sensor*EEG_brain_noise + (1-truth.snr_sensor)*EEG_sensor_noise;  
EEG_baseline_data = filtfilt(b_high, a_high, EEG_baseline_data')';

%% Plotting
idata = 1;
if plotting && ibatch == 1

load('tools/cm17');

% freq_inds = (truth.bandpass(1)*2+1):(truth.bandpass(2)*2+1); 
% [psi, stdpsi] = data2psi2(truth.sources_int', fs*2, fs*4, freq_inds);
% figure; imagesc(psi./stdpsi); colorbar
% max(max(abs(psi./stdpsi)))

k = randperm(n_source); k = k(1); % source index
    
% k-source amplitude distribution
ma = max(truth.source_amp(:, k));
allplots_cortex(sa, truth.source_amp(sa.cortex2K.in_to_cortex75K_geod, k), ...
    [0 ma], cm17a, 'A.U.', 1, ['figures/' foldername '/truth/dataset' num2str(idata) '/source1']);
          
% all source amplitude distributions
ma = max(sum(truth.source_amp, 2));
allplots_cortex(sa, sum(truth.source_amp(sa.cortex2K.in_to_cortex75K_geod, :), 2), ...
    [0 ma], cm17a, 'A.U.', 1, ['figures/' foldername '/truth/dataset' num2str(idata) '/sources']);

% dipole patterns
ma = max(abs(truth.EEG_field_pat(:, k)));
allplots_head(sa, sa.EEG_elec2head*truth.EEG_field_pat(:, k), [-ma ma], cm17, 'A.U.', ['figures/' foldername '/truth/dataset' num2str(idata) '/EEG_pat1'], sa.EEG_locs_3D(:, 1:3));
    
ma = max(abs(sum(truth.field_pat, 2)));
allplots_head(sa, sum(truth.field_pat, 2), [-ma ma], cm17, 'A.U.', ['figures/' foldername '/truth/dataset' num2str(idata) '/pats']);

export_fig(['figures/' foldername '/truth/dataset' num2str(idata) '/EEG_psd'], '-r150', '-a2');

% plot svd spectrum
P1 = svd(EEG_data);      
Pb = svd(EEG_baseline_data);
figure; subplot(3, 1, [1 2])
plot(P1, 'linewidth', 2)
hold on    
plot(Pb, 'r', 'linewidth', 2)
set(gca, 'fontsize', 18)
ylabel('Singular value')
xlabel('Data component')
title('EEG')
grid on
legend(truth.dataset, 'Baseline')
axis tight
export_fig(['figures/' foldername '/truth/dataset' num2str(idata) '/EEG_svd'], '-r150', '-a2');


% plot power spectrum and EEG data compare to baseline

    figure;
    
    no = sqrt(sum(EEG_data.^2, 2));
    [~, in_no] = max(no); % choosing channel with maximum norm eeg (===REVISE===)
    ss = std(EEG_data(in_no, :));
    [P1, f1] = pwelch(EEG_data(in_no, :)/ss, hanning(fs), [], fs, fs);     
    [Pb, fb] = pwelch(EEG_baseline_data(in_no, :)/ss, hanning(fs), [], fs, fs);
    
    subplot(3, 3, [1, 4])
    semilogy(f1, P1, 'linewidth', 2)
    hold on     
    semilogy(fb, Pb, 'r', 'linewidth', 2)
    set(gca, 'fontsize', 18)
    ylabel('Power [dB]')
    xlabel('Frequency [Hz]')
    title(['EEG ',num2str(in_no)])
    grid on
    legend(truth.dataset, 'Baseline')
    axis tight
    xlim([0 45])

    % plot normalize time series
    subplot(3, 3, 7)
    ss = std(EEG_data(in_no, 1:1000));
    plot((1:1000)/fs, 10 + EEG_data(in_no, 1:1000)/ss)
    hold on     
    plot((1:1000)/fs, EEG_baseline_data(in_no, 1:1000)/ss, 'r')
    set(gca, 'fontsize', 10, 'ytick', [])
    xlabel('Time [s]')
    ylabel('Voltage [AU]')
    grid on
    legend(truth.dataset, 'Baseline')
    axis tight
    
    no_temp = no;
    disp(in_no)
    
    % plot from the 2nd and 3rd largest index in 'no'
    
    for i = 1:2
    
    % plot power spectrum
    no_temp(in_no) = -inf;
    [~, in_no] = max(no_temp); % choosing channel with maximum norm eeg (===REVISE===)
    disp(in_no)
    ss = std(EEG_data(in_no, :));
    [P1, f1] = pwelch(EEG_data(in_no, :)/ss, hanning(fs), [], fs, fs);     
    [Pb, fb] = pwelch(EEG_baseline_data(in_no, :)/ss, hanning(fs), [], fs, fs);
    
    subplot(3, 3, [i+1, i+4])
    semilogy(f1, P1, 'linewidth', 2)
    hold on     
    semilogy(fb, Pb, 'r', 'linewidth', 2)
    set(gca, 'fontsize', 18)
    ylabel('Power [dB]')
    xlabel('Frequency [Hz]')
    title(['EEG ',num2str(in_no)])
    grid on
    legend(truth.dataset, 'Baseline')
    axis tight
    xlim([0 45])

    % plot normalize time series
    subplot(3, 3, i+7)
    ss = std(EEG_data(in_no, 1:1000));
    plot((1:1000)/fs, 10 + EEG_data(in_no, 1:1000)/ss)
    hold on     
    plot((1:1000)/fs, EEG_baseline_data(in_no, 1:1000)/ss, 'r')
    set(gca, 'fontsize', 10, 'ytick', [])
    xlabel('Time [s]')
    ylabel('Voltage [AU]')
    grid on
    legend(truth.dataset, 'Baseline')
    axis tight
        
    end
    
    % plot source spectrum from source data (chossing source randomly)
    
    figure;
    k = randi(length(sys.PARAMETER.ind_active));
    [Psource, fsource] = pwelch(truth.sources(sys.PARAMETER.ind_active(k), :)/ss, hanning(fs), [], fs, fs);
%     [Psource, fsource] = pwelch(testdata/ss, hanning(fs), [], fs, fs);
    subplot(2, 1, [1, 2])
    semilogy(fsource, Psource, 'r', 'linewidth', 2)
    set(gca, 'fontsize', 18)
    ylabel('Power [dB]')
    xlabel('Frequency [Hz]')
    title(['Source ',num2str(sys.PARAMETER.ind_active(k))])
    grid on
    legend(truth.dataset, 'Baseline')
    axis tight
    xlim([0 45])

    export_fig(['figures/' foldername '/truth/dataset' num2str(idata) '/EEG_timeseries'], '-r150', '-a2');

end

end
%================== End batch generation loop =============================

%% Saving data

eegdata = struct('EEG_data',EEG_data,'source_data',sources,'sampling_frequency',fs,'snr_source',truth.snr,...
                 'snr_sensor',truth.snr_sensor);
save(['data/',foldername,'/eegdata_', postfix], 'eegdata');

end
%================== End function ==========================================