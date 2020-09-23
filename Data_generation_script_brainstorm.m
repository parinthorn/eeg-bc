
%% data generation for brainstorm & our method
clc
if exist('sa','var')
    clearvars -except sa
else
    clear
    load('G:\Shared drives\MASTER_DRIVE\Journal\DATA_GENERATION\BBCB_code\data\sa')
end
foldername = 'experiment_paper_L_augment';
nbatch = 1;
n_model = 100;
n_deep = [2];
density = [10];


% stop at eegdata__act_20_deep_0_68
for n_act = density
    for nd = n_deep
        for mm=1:n_model
            postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
            disp(postfix)
            generate_eegdata(nbatch,foldername,0,postfix,sa,1,nd,n_act)
        end
    end
end
