function generate_eegdata(nbatch,foldername,plotting)
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
%                           0 for no plotting (and no saving any plotting)
%
%   OUTPUT:
%       No direct output from this function. But the model parameters and
%       EEG data will be saved in the 'data' folder with the directory and
%       files described as below;
%
%   data --> 'foldername' 
%                |
%                |-> eegdata.mat   :   contain eeg and source data matrices
%                |-> model.mat     :   contain model parameters
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
addpath tools
addpath tools/matlab_bgl
addpath tools/export_fig
addpath TestingCode
addpath pvo_subspace/subfun

% load head model and some miscallaneous data
load('data/sa')
load('data/miscdata')

% create directory to store data in
mkdir(['data/' foldername])

%%
% initialize random number generator
rng('default');
rng('shuffle');
rng_seed = rng;

truth.dataset = foldername;

% parameter of state-space model of sources REVISE
n_source = 10; %(m in the document)
n_source_active = 5;
PARAMETER.m = n_source;                   % #SOURCES
PARAMETER.m_active = n_source_active;            % #Active sources
PARAMETER.m_inactive = PARAMETER.m-PARAMETER.m_active;  % #Inactive sources (Redundant)
PARAMETER.lag = 2;                  % #lag
PARAMETER.nbutter = 2;              % order of Butterworth bandpass filter
PARAMETER.density = 0.8;            % Density of GC matrix
PARAMETER.group_density = 0.8;
PARAMETER.sigma_ar = 1*eye(PARAMETER.m_active);       % VAR process noise
PARAMETER.rng_seed = rng_seed;      %  save seed of the random number generator

% number of electrodes
EEG_M = length(sa.EEG_clab_electrodes);
PARAMETER.r = EEG_M;

% number of cluster of sources (must be less than n_source)
n_source_cluster = 3;
PARAMETER.n_source_cluster = n_source_cluster;

% spatial standard deviation of the sources (along cortical manifold) in mm
truth.sigma_range = [10 40];

% SNR range to sample from
truth.snr_range = [0.1 0.9];

% sample snr
%   truth.snr = truth.snr_range(1) + diff(truth.snr_range)*rand(1); % Can be changed later (for randomly choosing snr value)
truth.snr = 0.95;

% SNR for sensor
truth.snr_sensor = 0.95;

%number of biological noise source
% n_noise_sources = 500;

% sampling frequency
fs = 100;
PARAMETER.fs = fs;

% length of recording in sec
truth.len = 3*60; 

% resulting number of samples
N = fs*truth.len;

% band of interest, in which interaction takes place
truth.bandpass = [8 13];

% highpass filter coefficients %% == Why 0.1 hertz? can it be changed?
[b_high a_high] = butter(3, 1/fs*2, 'high');

% bandpass filter coefficients
if ~isempty(truth.bandpass)
[b_band a_band] = butter(3, truth.bandpass/fs*2);
end

% loop over datasets to be generated

% create folder for dataset idata
%   mkdir(['data/' foldername '/EEG/dataset_' num2str(idata)])
%   mkdir(['data/' foldername '/truth/dataset_' num2str(idata)]

%% Model generation

% generate source time series restricted to band of interest with the
% parameter of the model

sys = gengcss_eeg(PARAMETER,truth.bandpass,fs);  % system of source dynamics

%% spatial structure definition

ind_cluster_source = sys.PARAMETER.ind_cluster_source;

% sample source rois
ind_roi_cluster = randperm(8); ind_roi_cluster(n_source_cluster+1:end) = []; % assign the ROI number to each cluster of source (size = n_source_cluster x 1)

% assign ROI to each source. Sources in the same cluster are assigned to the same ROI
truth.in_roi = zeros(1,n_source);
for i = 1:n_source_cluster-1
truth.in_roi(ind_cluster_source(i):ind_cluster_source(i+1)-1) = ind_roi_cluster(i);
end
truth.in_roi(ind_cluster_source(n_source_cluster):n_source) = ind_roi_cluster(n_source_cluster);

truth.rois = rois(truth.in_roi);

% sample source extents
truth.sigmas = truth.sigma_range(1) + (truth.sigma_range(2)-truth.sigma_range(1))*rand(n_source, 1);

% sample source centers within rois
for i = 1:n_source
  truth.in_centers(i) = inds_roi_inner_2K{truth.in_roi(i)}(ceil(length(inds_roi_inner_2K{truth.in_roi(i)})*rand(1)));
end
truth.centers = sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K(truth.in_centers), :);

% calculate source amplitude distribution
cortex2K = [];
cortex2K.vc = sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K, :);
cortex2K.tri = sa.cortex2K.tri;

% generate Gaussian distributions on cortical manifold % loop 1:n_sources
for i = 1:n_source
[~, truth.source_amp(:, i)] = graphrbf(cortex2K, truth.sigmas(i), truth.in_centers(i));
end

% set activity outside of source octant to zero %%%==========================================
for i = 1:n_source
    % index that is not in INDS_ROI_OUTER_2K
    ind_tmp = setdiff(1:size(truth.source_amp,1), inds_roi_outer_2K{ truth.in_roi(i) });   
    truth.source_amp( ind_tmp, i) = 0;
end

for k=1:n_source
  TMP = norm(truth.source_amp(:,k)); % since the source of some channel can be entirely zero
  if TMP > 1e-10
    truth.source_amp(:, k) = truth.source_amp(:, k) ./ TMP;
  end
end

% obtain the lead field matrix
truth.EEG_field_pat = sa.cortex75K.EEG_V_fem_normal(:, sa.cortex2K.in_from_cortex75K)*truth.source_amp; % (n_electrode x 2K) x (2K x n_source)
sys.L0 = truth.EEG_field_pat;
sys.L = zeros([size(sys.L0),nbatch]);

%% Data generation

% pre-allocate source data and EEG data matrices
truth.sources = zeros(n_source,N,nbatch);
EEG_data = zeros(PARAMETER.r,N,nbatch);
sys.source_model = cell(nbatch,1);
sys.PARAMETER.pinknoise_cov = zeros(PARAMETER.m,PARAMETER.m,nbatch);
sys.F = zeros(PARAMETER.m,PARAMETER.m,nbatch);
factorC = zeros(nbatch,1);


%=================== start batch generation loop ==========================

for ibatch = 1:nbatch

%% Sources data generation

% generate sources data from the source system and corrupt inactive sources with pink noise
[truth.sources(:,:,ibatch),sys.PARAMETER.pinknoise_cov(:,:,ibatch),factorC(ibatch)] = getdatass_pinknoise(sys,N,truth.snr);

% recalculate F using the sample covariance of the generated pink noise
sys.F(:,:,ibatch) = calgcss(sys.source_model0.A,sys.source_model0.C,sys.PARAMETER.sigma_w,sys.PARAMETER.pinknoise_cov(:,:,ibatch));

sys.source_model{ibatch} = sys.source_model0; 
sys.source_model{ibatch}.C = sys.source_model{ibatch}.C*factorC(ibatch);

%% time series generation

% multiple the corrupted sources data with the lead field matrix
EEG_brain_signal_noise = truth.EEG_field_pat*truth.sources(:,:,ibatch);  % the pink noise is already corrupted in the source
sys.L(:,:,ibatch) = sys.L0*truth.snr/norm(EEG_brain_signal_noise,'fro'); % normalized

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

eegdata = struct('EEG_data',EEG_data,'source_data',truth.sources,'sampling_frequency',fs,'snr_source',truth.snr,...
                 'snr_sensor',truth.snr_sensor);
             
model = struct('PARAMETER',sys.PARAMETER,'source_model0',sys.source_model0,'source_model',sys.source_model,...
                'F0',sys.F0,'F',sys.F);

save(['data/',foldername,'/eegdata'], 'eegdata');
save(['data/' foldername '/model'], 'model');

end
%================== End function ==========================================