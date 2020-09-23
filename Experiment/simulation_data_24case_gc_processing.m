%%% Simulation data experiment: 3 deep source density: 0%,50%,75%
%%%                             2 active source density: 20%, 40%
%%%                             4 No. of electrodes: 108, 61, 32, 19
%%% This script iterates through the data in foldername with name eegdata_postfix
%%% with groundtruth model_postfix
%%% Written by Parinthorn Manomaisaowapak

clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_fix\';

%%
Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
n_deep = [0 2 3];
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
                disp(postfix)
                load([foldername,'debug_result_',postfix,'_',LL{1}],'sys_est')
                r = size(L,1);
                for tt=1:nbatch
                    C_ind =  sys_est(tt).sys.select_index(1); % bic
                    A = sys_est(tt).sys.A;
                    C = sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind); %bic:1,
                    E = sys_est(tt).C_out.E_Lpq_CLS(:,:,1);
                    W = sys_est(tt).C_out.W_Lpq_CLS(:,:,1);
                    [N,V,obj] = noisecovest(E,L,'homo');

                    [Ftmp,~,~] = calgcss(A,C,W,N,[]);
                    F{mm,tt} = sparse(Ftmp);
                    F_nz_ind{mm,tt} = find(Ftmp);
                end
                F_true{mm,1} = model.F0;
                F_true_ind{mm,1} = find(model.F0);
            end
            M.F_true{deep_it,density_it} = F_true;
            M.F_true_ind{deep_it,density_it} = F_true_ind;
            M.F_nz_ind{deep_it,density_it,n_elec} = F_nz_ind;
            M.F{deep_it,density_it,n_elec} = F;
        end
    end
end

%%
% save('Fresult','M')
