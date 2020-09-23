%%% Simulation data experiment: 50% deep source
%%%                             20% active source
%%%                             61 No. of electrodes
%%% This script iterates through the data in 'foldername' with name eegdata_'postfix'
%%% with groundtruth model_'postfix'
%%% Written by Parinthorn Manomaisaowapak

clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\';
nbatch = 1;
n_model = 100;
n_deep = [2];
density = [10];
%% Divide works
for n_act = density
    for nd = n_deep
        for mm = 1:n_model
            postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
            disp(postfix)
            load([foldername,'eegdata_',postfix])
            load([foldername,'model_',postfix])
            data = eegdata.EEG_data;
            L = model.L10;
            ind_chan = model.L10_ind;
            ind_Ctrue = find(norms(model.source_model0.C,2,2));
            mtilde = size(L,2);
            Timepoints = 15000;
            tic;
            for tt=1:nbatch % trial iteration
                data_in = data(ind_chan,end-Timepoints+1:end,tt);
                [sys_est(tt).sys,sys_est(tt).C_out] = ...
                    subid_eeg_Lpq(data_in ...
                    ,L,ind_Ctrue,10,-1);
                disp([foldername,'result_',postfix,'_',int2str(tt)])
            end
            toc;
            save([foldername,'result_',postfix,'_m50'],'sys_est')
            clear sys_est
        end
    end
end
