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
                [~,permind{mm,1}]=sort([unique(model.truth.in_roi),unique(model.truth.augment_data.in_roi)]);

            end
        end
    end
end
