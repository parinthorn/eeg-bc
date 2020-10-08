%This script collect all source selection result in a cell format same as
%the saved GC matrices in each experiment
% Written by Parinthorn Manomaisaowapak
%
% Licenses
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% Experiment: SIMEEG_VARYFACTOR
clear
clc
savepath = './saved_results_toplot/';
modelpath = './input_data/eeg_simulated_data/';
resultpath = './saved_experiment_results/simeeg_varyfactor/';
Lname = {'L0','L10','L20_32','L20'};
nbatch = 1;
n_model = 100;
n_deep = [0 2 4];
n_active_source = [10 20];
density_it = 0;

for n_act = n_active_source
    density_it = density_it+1;
    deep_it = 0;
    for nd = n_deep
        deep_it = deep_it+1;
        n_elec = 0;
        for LL=1:length(Lname)
            n_elec = n_elec+1;
            C_ind_bic = cell(n_model,nbatch);
            C_bic = cell(n_model,nbatch);
            C_nz_ind_bic = cell(n_model,nbatch);
            C = cell(n_model,nbatch);
            C_nz_ind = cell(n_model,nbatch);
            C_true = cell(n_model,1);
            C_true_ind = cell(n_model,1);
            for mm = 1:n_model
                postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
                load([modelpath,'model_',postfix])
                L = model.(Lname{LL});
                disp(postfix)
                load([resultpath,'simeeg_varyfactor_sysest_',postfix,'_',Lname{LL}],'sys_est')
                r = size(L,1);
                for tt=1:nbatch
                    C_ind_bic{mm,tt} =  sys_est(tt).C_out.ind_chosen_Lpq.bic;% use bic score
                    C_bic{mm,tt} = sparse(sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind_bic{mm,tt})); %bic:1,
                    C_nz_ind_bic{mm,tt} = sys_est(tt).C_out.nz_ind_C_Lpq{C_ind_bic{mm,tt}};
                    C{mm,tt} =  sparse(sys_est(tt).C_out.C_Lpq_CLS);
                    C_nz_ind{mm,tt} = sys_est(tt).C_out.nz_ind_C_Lpq;
                end
                C_true{mm} = model.source_model0.C;
                C_true_ind{mm} = model.PARAMETER.ind_active;%find(norms(model.source_model0.C,2,2));
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
% save([savepath,'result_C_100realization_with_bic_index'],'M')

% get accuracy
density_it = 0;
for n_act = n_active_source
    density_it = density_it+1;
    deep_it = 0;
    for nd = n_deep
        deep_it = deep_it+1;
        n_elec = 0;
        for LL=Lname
            n_elec = n_elec+1;
            tmp = compare_C(M.C_nz_ind{deep_it,density_it,n_elec},M.C_true_ind{deep_it,density_it},50,75);
            FPR_raw{deep_it,density_it,n_elec} = tmp.FPR_all;
            TPR_raw{deep_it,density_it,n_elec} = tmp.TPR_all;
            ACC_raw{deep_it,density_it,n_elec} = tmp.ACC_all;
        end
    end
end
RR.FPR = FPR_raw;
RR.TPR = TPR_raw;
RR.ACC = ACC_raw;
RR.bic_index = M.bic_index;
% save([savepath,'result_C_100realization_raw_accuracy'],'RR')

clear RR
density_it = 0;
for n_act = n_active_source
    density_it = density_it+1;
    deep_it = 0;
    for nd = n_deep
        deep_it = deep_it+1;
        n_elec = 0;
        for LL=Lname
            n_elec = n_elec+1;
            tmp = compare_C(M.C_nz_ind_bic{deep_it,density_it,n_elec},M.C_true_ind{deep_it,density_it},50,-1);
            RR.FPR_array{deep_it,density_it,n_elec} = tmp.FPR_all;
            RR.TPR_array{deep_it,density_it,n_elec} = tmp.TPR_all;
            RR.ACC_array{deep_it,density_it,n_elec} = tmp.ACC_all;
            RR.FPR_array_mean{deep_it,density_it,n_elec} = mean(tmp.FPR_all);
            RR.TPR_array_mean{deep_it,density_it,n_elec} = mean(tmp.TPR_all);
            RR.ACC_array_mean{deep_it,density_it,n_elec} = mean(tmp.ACC_all);
            RR.FPR_array_std{deep_it,density_it,n_elec} = std(tmp.FPR_all);
            RR.TPR_array_std{deep_it,density_it,n_elec} = std(tmp.TPR_all);
            RR.ACC_array_std{deep_it,density_it,n_elec} = std(tmp.ACC_all);
        end
    end
end
% save([savepath,'result_C_100realization_bic_accuracy.mat'],'RR')

%% Experiment: SIMEEG_SOURCERECON
clear
clc

savepath = './saved_results_toplot/';
modelpath = './input_data/source_recon_brainstorm/';
resultpath = './saved_experiment_results/simeeg_sourcerecon/';
nbatch = 1;
n_model = 100;
n_deep = 2;
n_active_source = 10;

C_ind_bic = cell(n_model,nbatch);
C_bic = cell(n_model,nbatch);
C_nz_ind_bic = cell(n_model,nbatch);
C = cell(n_model,nbatch);
C_nz_ind = cell(n_model,nbatch);
C_true = cell(n_model,1);
C_true_ind = cell(n_model,1);
for mm = 1:n_model
    postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
    load([modelpath,'model_',postfix])
    L = model.L10;
    disp(postfix)
    load([resultpath,'simeeg_sourcerecon_sysest_',postfix],'sys_est')
    r = size(L,1);
    for tt=1:nbatch
        C_ind_bic{mm,tt} =  sys_est(tt).C_out.ind_chosen_Lpq.bic;% use bic score
        C_bic{mm,tt} = sparse(sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind_bic{mm,tt})); %bic:1,
        C_nz_ind_bic{mm,tt} = sys_est(tt).C_out.nz_ind_C_Lpq{C_ind_bic{mm,tt}};
        C{mm,tt} =  sparse(sys_est(tt).C_out.C_Lpq_CLS);
        C_nz_ind{mm,tt} = sys_est(tt).C_out.nz_ind_C_Lpq;
    end
    C_true{mm} = model.source_model0.C;
    C_true_ind{mm} = model.PARAMETER.ind_active;
end
M.C_true{1} = C_true;
M.C_true_ind{1} = C_true_ind;
M.C_nz_ind_bic{1} = C_nz_ind_bic;
M.C_bic{1} = C_bic;
M.bic_index{1} = C_ind_bic;
M.C_nz_ind{1} = C_nz_ind;
M.C{1} = C;


% save([savepath,'simeeg_sourcerecon_C'],'M')

C_ind_bic = cell(n_model,nbatch);
C_bic = cell(n_model,nbatch);
C_nz_ind_bic = cell(n_model,nbatch);
C = cell(n_model,nbatch);
C_nz_ind = cell(n_model,nbatch);
C_true = cell(n_model,1);
C_true_ind = cell(n_model,1);
for mm = 1:n_model
    postfix = ['_act_',int2str(n_act),'_deep_',int2str(nd),'_',int2str(mm)];
    load([modelpath,'model_',postfix])
    L = model.L10;
    disp(postfix)
    load([resultpath,'simeeg_sourcerecon_sysest_',postfix,'_m50'],'sys_est')
    r = size(L,1);
    for tt=1:nbatch
        C_ind_bic{mm,tt} =  sys_est(tt).C_out.ind_chosen_Lpq.bic;% use bic score
        C_bic{mm,tt} = sparse(sys_est(tt).C_out.C_Lpq_CLS(:,:,C_ind_bic{mm,tt})); %bic:1,
        C_nz_ind_bic{mm,tt} = sys_est(tt).C_out.nz_ind_C_Lpq{C_ind_bic{mm,tt}};
        C{mm,tt} =  sparse(sys_est(tt).C_out.C_Lpq_CLS);
        C_nz_ind{mm,tt} = sys_est(tt).C_out.nz_ind_C_Lpq;
    end
    C_true{mm} = model.source_model0.C;
    C_true_ind{mm} = model.PARAMETER.ind_active;%find(norms(model.source_model0.C,2,2));
end
M.C_true{1} = C_true;
M.C_true_ind{1} = C_true_ind;
M.C_nz_ind_bic{1} = C_nz_ind_bic;
M.C_bic{1} = C_bic;
M.bic_index{1} = C_ind_bic;
M.C_nz_ind{1} = C_nz_ind;
M.C{1} = C;


% save([savepath,'result_C_augment_L_m50_100realization_with_bic_index'],'M')
% This result is not reported because it is the same with SIMEEG_VARYFACTOR
% experiment

load('simeeg_sourcerecon_C.mat')
tmp = compare_C(M.C_nz_ind{1},M.C_true_ind{1},100,75);
FPR_raw{1,1,1} = tmp.FPR_all;
TPR_raw{1,1,1} = tmp.TPR_all;
ACC_raw{1,1,1} = tmp.ACC_all;
RR.FPR = FPR_raw;
RR.TPR = TPR_raw;
RR.ACC = ACC_raw;
RR.bic_index = M.bic_index;
% save([savepath,'result_C_augment_L_100realization_raw_accuracy_fullcompare'],'RR')

gt_source_nz_ind = cell(n_model,GridSize);
GridSize = size(M.C_nz_ind{1},2); % number of regularization level
for ii=1:n_model
    for jj=1:GridSize
        tmp = M.C_nz_ind{1}{ii}{jj};
        tmp(tmp>50) = [];
        gt_source_nz_ind{ii}{jj} = tmp;
    end
    
end
tmp = compare_C(gt_source_nz_ind,M.C_true_ind{1},50,75);
FPR_raw{1} = tmp.FPR_all;
TPR_raw{1} = tmp.TPR_all;
ACC_raw{1} = tmp.ACC_all;
RR.FPR = FPR_raw;
RR.TPR = TPR_raw;
RR.ACC = ACC_raw;
RR.bic_index = M.bic_index;
% save([savepath,'result_C_augment_L_100realization_raw_accuracy_partial_compare'],'RR')
