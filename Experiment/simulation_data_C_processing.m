%% Simulation data experiment: 50% deep source
%%                             20% active source
%%                             61 No. of electrodes
%% This script rearrange matrix C to be in the accuracy evaluation format
%% Written by: PARINTHORN MANOMAISAOWAPAK
clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\';

%%

Lname = {'L10'};
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
%                 load([foldername,'result_',postfix],'sys_est')
                load([foldername,'result_',postfix,'_m50'],'sys_est')
                for tt=1:nbatch
                    C_ind_bic{mm,tt} =  sys_est(tt).sys.select_index(1); % bic
                    C_bic{mm,tt} = sparse(sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind_bic{mm,tt})); %bic:1,
                    C_nz_ind_bic{mm,tt} = find(norms(C_bic{mm,tt},2,2));

                    C{mm,tt} =  (sys_est(tt).C_out.C_Lpq_CLS);
                    parfor kk=1:75
                        tmpC{kk,1} = find(norms(sys_est(tt).C_out.C_Lpq_CLS(:,:,kk),2,2));
                    end
                    C_nz_ind{mm,tt} = tmpC;
                end
                C_true{mm,1} = model.source_model0.C;
                C_true_ind{mm,1} = find(norms(model.source_model0.C,2,2));
            end
            M.C_true{deep_it,density_it} = C_true;
            M.C_true_ind{deep_it,density_it} = C_true_ind;
            M.C_nz_ind_bic{deep_it,density_it,n_elec} = C_nz_ind_bic;
            M.C_bic{deep_it,density_it,n_elec} = C_bic;
            M.bic_index{deep_it,density_it,n_elec} = C_ind_bic;

            M.C_nz_ind{deep_it,density_it,n_elec} = C_nz_ind;
            M.C{deep_it,density_it,n_elec} = C;
        end
    end
end
% save('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\result_C_augment_L_100realization_with_bic_index','M')
save('G:\Shared drives\MASTER_DRIVE\Journal\Experiment_main\result_C_augment_L_m50_100realization_with_bic_index','M')
%%


clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_fix\';
Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
% n_deep = [0 2 3];
n_deep = [2];
density = [10 20];

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
                L = model.(LL{1});
                disp([postfix])
                load([foldername,'debug_result_',postfix,'_',LL{1}],'sys_est')
                r = size(L,1);
                for tt=1:nbatch
                    C_ind_bic{mm,tt} =  sys_est(tt).sys.select_index(1); % bic
                    C_bic{mm,tt} = sparse(sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind_bic{mm,tt})); %bic:1,
                    C_nz_ind_bic{mm,tt} = find(norms(C_bic{mm,tt},2,2));

                    C{mm,tt} =  (sys_est(tt).C_out.C_Lpq_CLS);
                    parfor kk=1:75
                        tmpC{kk,1} = find(norms(sys_est(tt).C_out.C_Lpq_CLS(:,:,kk),2,2));
                    end
                    C_nz_ind{mm,tt} = tmpC;
                end
                C_true{mm,1} = model.source_model0.C;
                C_true_ind{mm,1} = find(norms(model.source_model0.C,2,2));
            end
            M.C_true{deep_it,density_it} = C_true;
            M.C_true_ind{deep_it,density_it} = C_true_ind;
            M.C_nz_ind_bic{deep_it,density_it,n_elec} = C_nz_ind_bic;
            M.C_bic{deep_it,density_it,n_elec} = C_bic;
            M.bic_index{deep_it,density_it,n_elec} = C_ind_bic;

            M.C_nz_ind{deep_it,density_it,n_elec} = C_nz_ind;
            M.C{deep_it,density_it,n_elec} = C;
        end
    end
end
%

save('debug_result_C_100realization_with_bic_index','M')
