%% This script generates heat map of reconstructed source from Brainstorm
% The EEG data is generated under the setting: 50% deep source, 20% active source, 61 No. of electrodes
% 
% input = eegdata, kernel matrix
% Method 
% step 1) use Kernel to reconstruct 75K source 
% step 2) take absolute averaging
% 
% output = averaged absolute of reconstructed source (over time)
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
dir_name = './saved_experiment_results/simeeg_sourcerecon/brainstorm_directory/data/'; % to load kernel matrix for source reconstruction
outpath = './saved_experiment_results/simeeg_sourcerecon/brainstorm_reconsource_amplitude/'; 

mkdir(outpath)

%% 
for ii=1:100
    tic;
    ker_dir = [dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\',dir([dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\results_MN_EEG_KERNEL*']).name];
    load(ker_dir,'ImagingKernel')
    KERNEL_MN = ImagingKernel;


    ker_dir = [dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\',dir([dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\results_PNAI_EEG_KERNEL*']).name];
    load(ker_dir,'ImagingKernel')
    KERNEL_LCMV = ImagingKernel;


    ker_dir = [dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\',dir([dir_name,'subject_',int2str(ii),'\@rawdata',int2str(ii),'\results_sLORETA_EEG_KERNEL*']).name];
    load(ker_dir,'ImagingKernel')
    KERNEL_sLORETA = ImagingKernel;

    load([inpath_eegdata,'data',int2str(ii)])
    block_prod = {[1:1:5000],[5001:10000],[10001:15000]};
    tmp_MN = zeros(size(ImagingKernel,1),1);
    tmp_LCMV = zeros(size(ImagingKernel,1),1);
    tmp_sLORETA = zeros(size(ImagingKernel,1),1);
    parfor tt=1:3
        tmp_MN = tmp_MN+sum(abs(KERNEL_MN*EEG_data(:,block_prod{tt})),2);
        tmp_LCMV = tmp_LCMV+sum(abs(KERNEL_LCMV*EEG_data(:,block_prod{tt})),2);
        tmp_sLORETA = tmp_sLORETA+sum(abs(KERNEL_sLORETA*EEG_data(:,block_prod{tt})),2);
    end
%      = KERNEL_MN*EEG_data;
    x.MN = tmp_MN/15000;
    x.LCMV = tmp_LCMV/15000;
    x.sLORETA = tmp_sLORETA/15000;

    save([outpath,'map_',int2str(ii)],'x')
    toc;
    disp(ii)
    clear x
end
