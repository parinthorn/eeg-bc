%% This script runs the experiment of performing GC estimation compared to Brainstorm
% The EEG data is generated under the setting: 50% deep source, 20% active source, 61 No. of electrodes
% This sctrip converts the reconstructed source from 75K resolution to ROI
% 
% input = reconstructed source in 75K AND matrix transform for averaging
% Method 
% step 1) use Kernel to reconstruct 75K source 
% step 2) average reconstructed 75K source to ROI resolution
% 
% output = averaged reconstructed source in ROI
% 
% Written by Parinthorn Manomaisaowapak
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
clear
clc

inpath_eegdata = './input_data/source_recon_brainstorm/'; % EEG data in brainstorm format
inpath_kernel = './saved_experiment_results/simeeg_sourcerecon/brainstorm_directory/Kernel_roi'; % matrix converts from 75K to ROI resolution
dir_name = './saved_experiment_results/simeeg_sourcerecon/brainstorm_directory/data/'; % to load kernel matrix for source reconstruction
outpath = './saved_experiment_results/simeeg_sourcerecon/brainstorm_avg_recon_source/'; 

%% 
load(inpath_kernel);
mkdir(outpath)
for ii=1:100
    tic;
    ker_dir = [dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\',dir([dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\results_MN_EEG_KERNEL*']).name];
    load(ker_dir,'ImagingKernel')
    eeg_to_roi_MN = bsxfun(@rdivide,(Kernel_roi)*ImagingKernel,sum(Kernel_roi,2));%(diag(1./sum(Kernel_roi,2)))*(Kernel_roi)*ImagingKernel;
    ker_dir = [dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\',dir([dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\results_PNAI_EEG_KERNEL*']).name];
    load(ker_dir,'ImagingKernel')
    eeg_to_roi_LCMV = bsxfun(@rdivide,(Kernel_roi)*ImagingKernel,sum(Kernel_roi,2));%(diag(1./sum(Kernel_roi,2)))*(Kernel_roi)*ImagingKernel;
    ker_dir = [dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\',dir([dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\results_sLORETA_EEG_KERNEL*']).name];
    load(ker_dir,'ImagingKernel')
    eeg_to_roi_sLORETA = bsxfun(@rdivide,(Kernel_roi)*ImagingKernel,sum(Kernel_roi,2));%(diag(1./sum(Kernel_roi,2)))*(Kernel_roi)*ImagingKernel;

    load([inpath_eegdata,'data',int2str(ii)])
    x.MN = eeg_to_roi_MN*EEG_data;
    x.LCMV = eeg_to_roi_LCMV*EEG_data;
    x.sLORETA = eeg_to_roi_sLORETA*EEG_data;
    save([outpath,'recon_source',int2str(ii)],'x')
    toc;
    disp(ii)
    clear x
end
