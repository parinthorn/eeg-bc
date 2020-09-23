clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\';
load('F_result_roi.mat')
load('F_result_roi_permind.mat')
%%
E.permind = permind;
for ii=1:100

    load([foldername,'model_','_act_',int2str(10),'_deep_',int2str(2),'_',int2str(ii)])
    groundtruth_source_cluster_index = cell(model.PARAMETER.n_source_cluster,1);
    for kk=1:model.PARAMETER.n_source_cluster-1
        groundtruth_source_cluster_index{kk} = [model.PARAMETER.ind_cluster_source(kk):1:model.PARAMETER.ind_cluster_source(kk+1)-1];
    end
    groundtruth_source_cluster_index{end} = model.PARAMETER.ind_cluster_source(end):1:50;
    tmpf = @(x) x+50;
    estimate_source_cluster_index = [groundtruth_source_cluster_index; ...
        cellfun(tmpf,groundtruth_source_cluster_index,'UniformOutput',false)];
    permute_source_cluster_index = (estimate_source_cluster_index(permind{ii}));
    permute_source_cluster_index= [permute_source_cluster_index{:}];


    F_true_roi_tmp = [M.F_true_roi{1}{ii} zeros(4,4);zeros(4,8)];
    F_true_roi_tmp = F_true_roi_tmp.*(F_true_roi_tmp>1e-6);
    E.F_true_roi{1}{ii} = F_true_roi_tmp(permind{ii},permind{ii});
    E.F_true_roi_ind{1}{ii} = find(M.F_true_roi{1}{ii});
    F_roi_tmp = M.F_roi{1}{ii};
    F_roi_tmp = F_roi_tmp.*(F_roi_tmp>1e-6);
    E.F_roi{1}{ii} = F_roi_tmp(permind{ii},permind{ii});
    E.F_roi_nz_ind{1}{ii} = find(E.F_roi{1}{ii});

    F_true_tmp = [M.F_true{1}{ii} zeros(50,50);zeros(50,100)];
    E.F_true{1}{ii} = F_true_tmp(permute_source_cluster_index,permute_source_cluster_index);
    E.F{1}{ii} = M.F{1}{ii}(permute_source_cluster_index,permute_source_cluster_index);
    E.F_nz_ind{1}{ii} = {find(E.F{1}{ii})};
    E.permind_all{ii,1} = permute_source_cluster_index;
%     figure(1)
%     subplot(1,2,1);imagesc(E.F_true_roi{1}{ii}); subplot(1,2,2);imagesc(E.F_true{1}{ii})
%     figure(2)
%     subplot(1,2,1);imagesc(E.F_roi{1}{ii}); subplot(1,2,2);imagesc(E.F{1}{ii})
%     pause()
end
%%
% save('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\F_result_L_augment','M')
