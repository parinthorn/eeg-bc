%% Simulation data experiment: 50% deep source
%%                             20% active source
%%                             61 No. of electrodes
%% This script compute Granger causality matrix F and reorder the rois
%% Written by: PARINTHORN MANOMAISAOWAPAK

clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\';

%%
Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
n_deep = [2];
density = [10];
density_it = 0;
for n_act = density
    density_it = density_it+1;
    deep_it = 0;
    for nd = n_deep
        deep_it = deep_it+1;
        n_elec = 0;
        for LL=Lname
            n_elec = n_elec+1;

            for mm = 1:n_model
                postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
                load([foldername,'model_',postfix])
%                 load([foldername,'eegdata_',postfix])
%                 E = reduce_L(eegdata.L);
%                 L = E.(Lname{n_elec});
                L = model.(LL{n_elec});
                disp(postfix)
                load([foldername,'result_',postfix,'_',LL{n_elec}],'sys_est')
                r = size(L,1);
                At = model.source_model0.A;
                Ct = model.source_model0.C;
                Wt = model.PARAMETER.sigma_w;
                Nt = model.PARAMETER.pinknoise_cov;

                groundtruth_source_cluster_index = cell(model.PARAMETER.n_source_cluster,1);
                for ii=1:model.PARAMETER.n_source_cluster-1
                    groundtruth_source_cluster_index{ii} = [model.PARAMETER.ind_cluster_source(ii):1:model.PARAMETER.ind_cluster_source(ii+1)-1];
                end
                groundtruth_source_cluster_index{end} = model.PARAMETER.ind_cluster_source(end):1:50;
                F_roi_true = calgcss_block(At,Ct,Wt,Nt,[],groundtruth_source_cluster_index);
                tmpf = @(x) x+size(Ct,1);
%                 estimate_source_cluster_index = groundtruth_source_cluster_index;
                estimate_source_cluster_index = [groundtruth_source_cluster_index; ...
                    cellfun(tmpf,groundtruth_source_cluster_index,'UniformOutput',false)];
                for tt=1:nbatch
                    % estimate
                    C_ind =  sys_est(tt).sys.select_index(1); % bic
                    A = sys_est(tt).sys.A;
                    C = sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind); %bic:1,
                    E = sys_est(tt).C_out.E_Lpq_CLS(:,:,1);
                    W = sys_est(tt).C_out.W_Lpq_CLS(:,:,1);
                    [N,V,obj] = noisecovest(E,L,'homo');
                    [Ftmp,~,~] = calgcss(A,C,W,N,[]);
                    F{mm,tt} = sparse(Ftmp);
                    F_nz_ind{mm,tt} = find(Ftmp);



                    Ftmp_roi = calgcss_block(A,C,W,N,[],estimate_source_cluster_index);
                    F_roi{mm,tt} = sparse(Ftmp_roi);
                    F_roi_nz_ind{mm,tt} = find(Ftmp_roi);


                end
                F_true{mm,1} = model.F0;
                F_true_ind{mm,1} = find(model.F0);
                F_true_roi{mm,1} = F_roi_true;
            end
            M.F_true{deep_it,density_it} = F_true;
            M.F_true_roi{deep_it,density_it} = F_true_roi;
            M.F_true_ind{deep_it,density_it} = F_true_ind;
            M.F_nz_ind{deep_it,density_it,n_elec} = F_nz_ind;
            M.F{deep_it,density_it,n_elec} = F;
            M.F_roi{deep_it,density_it,n_elec} = F_roi;
            M.F_roi_nz_ind{deep_it,density_it,n_elec} = F_roi_nz_ind;
        end
    end
end

%%
% save('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\F_result_L_augment','M')

clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\';
load('F_result_roi.mat')
load('F_result_roi_permind.mat')
%%
for ii=1:100
    F_true_tmp = [M.F_true_roi{1}{ii} zeros(4,4);zeros(4,8)];
    E.F_true_roi{1}{ii} = F_true_tmp(permind{ii},permind{ii});
    E.F_true_roi_ind{1}{ii} = find(M.F_true_roi{1}{ii});
    F_roi_tmp = M.F_roi{1}{ii};
    E.F_roi{1}{ii} = F_roi_tmp(permind{ii},permind{ii});
    E.F_roi_nz_ind{1}{ii} = find(E.F_roi{1}{ii});
end
