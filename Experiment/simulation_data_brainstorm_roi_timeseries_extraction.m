%% average rois
clear
clc
dir_name = 'G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\experiment_brainstorm_dir\act_10_deep_2_nelec_61_constraint75k\data\';

load('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\Kernel_roi');
save_dir = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\data_roi_reconstruct_75K\';
mkdir(save_dir)
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

    load(['G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\data_for_brainstorm\data',int2str(ii)])
    x.MN = eeg_to_roi_MN*EEG_data;
    x.LCMV = eeg_to_roi_LCMV*EEG_data;
    x.sLORETA = eeg_to_roi_sLORETA*EEG_data;
    save([save_dir,'data',int2str(ii)],'x')
    toc;
    disp(ii)
    clear x
end
