%%% Simulation data experiment: 3 deep source density: 0%,50%,75%
%%%                             2 active source density: 20%, 40%
%%%                             4 No. of electrodes: 108, 61, 32, 19
%%% This script iterates through the data in foldername with name eegdata_postfix
%%% with groundtruth model_postfix
%%% Written by Parinthorn Manomaisaowapak
clear
clc
foldername = 'G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_fix\';
nbatch = 1;
n_model = 100;
n_deep = [0 2 4];
density = [10 20];

%% Divide works



Lname = {'L0','L10','L20_32','L20'};
for n_act = density
    for nd = n_deep
        for mm = 1:n_model
            postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
            disp(postfix)
            load([foldername,'eegdata_',postfix])
            load([foldername,'model_',postfix])
            data = eegdata.EEG_data;
            for LL=Lname
                L = model.(LL{1});
                r = size(L,1);
                if strcmp(LL{1},'L0')
                    ind_chan = 1:1:r;
                elseif strcmp(LL{1},'L20_32')
                    ind_chan = model.L20_ind_32;
                else
                    ind_chan = model.([LL{1},'_ind']);
                end
                ind_Ctrue = find(norms(model.source_model0.C,2,2));
                Timepoints = 15000;
                tic;
                for tt=1:nbatch % trial iteration
                    data_in = data(ind_chan,end-Timepoints+1:end,tt);
                    [sys_est(tt).sys,sys_est(tt).C_out] = ...
                        subid_eeg_Lpq_tmp(data_in ...
                        ,L,ind_Ctrue,10,-1);
                    disp([foldername,'result_',postfix,'_',LL{1},int2str(tt)])
                end
                toc;
                save([foldername,'debug_result_',postfix,'_',LL{1}],'sys_est')
                clear sys_est
            end
        end
    end
end
