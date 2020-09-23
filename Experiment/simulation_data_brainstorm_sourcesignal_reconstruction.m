%% Simulation data experiment: 50% deep source
%%                             20% active source
%%                             61 No. of electrodes
%% This script generate 75K dimensional source signal by multiplying
%% the kernel from brainstorm to the eeg data
%% Written by: PARINTHORN MANOMAISAOWAPAK
clear
clc
dir_name = 'G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\experiment_brainstorm_dir\act_10_deep_2_nelec_61_constraint75K\data\';

save_loc = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\map_reconstruct_75K\';
mkdir(save_loc)
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

    load(['G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\data_for_brainstorm\data',int2str(ii)])
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

    save([save_loc,'map_',int2str(ii)],'x')
    toc;
    disp(ii)
    clear x
end
