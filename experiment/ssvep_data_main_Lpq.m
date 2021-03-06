%% SSVEP EEG Connectivity experiment
% 
% This experiment is to estimate brain connectivity from SSVEP (Steady state visually evoked potential).
% SSVEP are signals that are natural responses to visual stimulation at specific frequencies.
% When the retina is excited by a visual stimulus ranging from 3.5 Hz to 75 Hz, 
% the brain generates electrical activity at the same (or multiples of) frequency of the visual stimulus.
% 
% We have dataset from Deirel and Pedro (and this was obtained from Turkish professor, Tamer). 
% The data are EEG 30 channels and the lead-field matrix are provided by Deirel
% 
% This experiment uses Source selection with Lpq penalty (non-convex case)
% User should change line 28 or 32 (select the number of sources in
% consideration)
% Written by Jitkomut Songsiri, Parinthorn Manomaisaowapak
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
%

clc;  clear all; close all;

inpath = './input_data/ssvep_eeg';
outpath = './saved_experiment_results/';

%% Model parameters

load('data_ssvep_eeg.mat')
load('ssvep_leadfield.mat') % Kindms stores lead-field matrix, computed from SPM software

% r = size(ssvep_eeg{1},1);
ntrial = length(ssvep_eeg); % data have 3 trials
m = size(L,2);
roi_label = {'OL-L','TL-L','FL-L','OL-R','TL-R','FL-R'}; % read from the Chula report

% choices of selecting 18 sources or 72 sources
% m = 18
sample_sources = [5 6 7 17 18 19 29 30 31 41 42 43  53 54 55 65 66 67];source_cluster_index = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14,15],[16,17,18]};tick_vect = 2:3:17; % 18 entries (3 sources per ROI)

% m = 72
% sample_sources = 1:1:72;source_cluster_index = {[1:12],[13:24],[25:36],[37:48],[49:60],[61:72]};tick_vect = 6:12:78;

%
mhat = length(sample_sources);
Lhat = L(:,sample_sources); % L is 30 x 18 (18 sources) or 30 x 72
select_chan = 1:1:29;
Lhat = Lhat(select_chan,:); % remove redundant row
r = length(select_chan);

%% plot the time series
channel_ind = 20; trial_ind = 2; %
figure(1);
for kk=1:ntrial
    subplot(ntrial,1,kk);plot(ssvep_eeg{trial_ind}(channel_ind,:));
    xlabel('Time');ylabel('EEG'); title(['Trial # ',num2str(kk)]);
end

%% Subspace identification on SSVEP EEG data
% figure(99);
nhat = -1; % or choose from our method
ii = max([ceil(2*nhat/r) 10]); % prevent maximum ii less than 5
% ii = 10;
A = cell(ntrial,1);
Chat = cell(ntrial,1);
W = cell(ntrial,1);
N = cell(ntrial,1);
for kk=1:ntrial
    [sys_est(kk).sys,sys_est(kk).C_out] = subid_eeg_Lpq(ssvep_eeg{kk}(select_chan,:),Lhat,ii,nhat);%    estimate 'ntrial' models
end

save([outpath,'ssvep_eeg_estimated_model_nonconvex.mat'],'sys_est');

%% Compute GC indices & plotting

load([outpath,'ssvep_eeg_estimated_model_nonconvex.mat']);

for kk=1:ntrial
bic_ind =  sys_est(kk).C_out.ind_chosen_Lpq.bic; %bic
aicc_ind =  sys_est(kk).C_out.ind_chosen_Lpq.aicc; %aicc

    A{kk} = sys_est(kk).sys.A;
    aicc.Chat{kk} = sys_est(kk).C_out.C_Lpq_CLS(:,:,aicc_ind);
    E = sys_est(kk).C_out.E_Lpq_CLS(:,:,aicc_ind);
    aicc.W{kk} = sys_est(kk).C_out.W_Lpq_CLS(:,:,aicc_ind);
    [aicc.N{kk},~,~] = noisecovest(E,Lhat,'homo');

    bic.Chat{kk} = sys_est(kk).C_out.C_Lpq_CLS(:,:,bic_ind);
    E = sys_est(kk).C_out.E_Lpq_CLS(:,:,bic_ind);
    bic.W{kk} = sys_est(kk).C_out.W_Lpq_CLS(:,:,bic_ind);
    [bic.N{kk},~,~] = noisecovest(E,Lhat,'homo');


    % A,C,W,N = subspace identification + source selection + noise covariance estimation
    % save estimated model
end

%==========AICC===================
aicc.Fhat = cell(ntrial,1);
aicc.Fhat_roi = cell(ntrial,1);
for kk=1:ntrial
    [aicc.Fhat{kk}] = calgcss(A{kk},aicc.Chat{kk},aicc.W{kk},aicc.N{kk},[]);
    [aicc.Fhat_roi{kk}] = calgcss_block(A{kk},aicc.Chat{kk},aicc.W{kk},aicc.N{kk},[],source_cluster_index);
% save estimated Fhat
end
aicc.Fhat_avg = (aicc.Fhat{1}+aicc.Fhat{2}+aicc.Fhat{3})/3;
aicc.Fhat_roi_avg = (aicc.Fhat_roi{1}+aicc.Fhat_roi{2}+aicc.Fhat_roi{3})/3; % average over 3 trials

%==========BIC===================
bic.Fhat = cell(ntrial,1);
bic.Fhat_roi = cell(ntrial,1);
for kk=1:ntrial
    [bic.Fhat{kk}] = calgcss(A{kk},bic.Chat{kk},bic.W{kk},bic.N{kk},[]);
    [bic.Fhat_roi{kk}] = calgcss_block(A{kk},bic.Chat{kk},bic.W{kk},bic.N{kk},[],source_cluster_index);
% save estimated Fhat
end
bic.Fhat_avg = (bic.Fhat{1}+bic.Fhat{2}+bic.Fhat{3})/3;
bic.Fhat_roi_avg = (bic.Fhat_roi{1}+bic.Fhat_roi{2}+bic.Fhat_roi{3})/3; % average over 3 trials


% RR.sys_est = sys_est;
for kk=1:ntrial
RR.C_Lpq{kk,1} = sys_est(kk).C_out.C_Lpq;
RR.C_Lpq_CLS{kk,1} = sys_est(kk).C_out.C_Lpq_CLS;
RR.aicc_ind{kk,1} = sys_est(kk).C_out.ind_chosen_Lpq.aicc;
RR.bic_ind{kk,1} = sys_est(kk).C_out.ind_chosen_Lpq.bic;
end
RR.gc.Fhat.aicc = aicc.Fhat;
RR.gc.Fhat_roi.aicc = aicc.Fhat_roi;
RR.gc.Fhat_avg.aicc = aicc.Fhat_avg;
RR.gc.Fhat_roi_avg.aicc = aicc.Fhat_roi_avg;
RR.gc.Fhat.bic = bic.Fhat;
RR.gc.Fhat_roi.bic = bic.Fhat_roi;
RR.gc.Fhat_avg.bic = bic.Fhat_avg;
RR.gc.Fhat_roi_avg.bic = bic.Fhat_roi_avg;

RR.sample_sources = sample_sources;
RR.roi_label = roi_label;
RR.select_chan = select_chan;

% uncomment if we want to save 
save([outpath,'ssvep_rawresult_Lpq.mat'],'RR');
%========================================== END COMPUTATION ======================================================
