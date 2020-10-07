
%% This script runs the experiment of performing GC estimation compared to Brainstorm
% The EEG data is generated under the setting: 50% deep source, 20% active source, 61 No. of electrodes
% Input: Estimated state-space parameters
% Output: Estimated GC matrices (mapped to Haufe ROI format)
% 
% Step:
% (1) calculate GC matrices from system estimated parameters. 
% This GC matrices are with the original reference of selected 4 ROIs (out of 8 ROIs)
% (2) Permute the calculated GC matrices to the uniform format compared to 8 ROIs and with the same ROI order as Haufe

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


clear
clc

datapath = './input_data/source_recon_brainstorm/eeg_simulated_data_deep50percent_active20percent';

inpath = './saved_experiment_results/simeeg_sourcerecon/';
outpath = './saved_experiment_results/simeeg_sourcerecon/';

n_model = 100;
n_deep = 2; % run only the case that there are 2 deep source ROIs
n_active_source = 10; % run the case that there are 10 active sources (from 50 sources)

load('./saved_experiment_results/simeeg_sourcerecon/simeeg_sourcerecon_source_perm_ind'); % load source cluster index and permutation mapping
groundtruth_source_cluster_index = sourceind.true_source_cluster_ind;
estimate_source_cluster_index = sourceind.est_source_cluster_ind;
perm_roi_ind = sourceind.perm_roi_ind;
perm_node_ind = sourceind.perm_node_ind;
%% Step (1): Calculate GC matrices (reference of selected 4 ROIs)
for mm = 1:n_model
    postfix = ['_act_',int2str(n_active_source),'_deep_',int2str(n_deep),'_',int2str(mm)];
    load([datapath,'model_',postfix])
    load([datapath,'eegdata_',postfix])
    E = reduce_L(eegdata.L); % Lhat, Lead field matrix of size (108 x 100)
    L = E.L10; % This experiment used 61 electrodes.L has size of 61 x 100
    disp(postfix)
    load([inpath,'simeeg_sourcerecon_sysest_',postfix],'sys_est')
    r = size(L,1);
    At = model.source_model0.A;
    Ct = model.source_model0.C;
    Wt = model.PARAMETER.sigma_w;
    Nt = model.PARAMETER.pinknoise_cov;
    n_all_source = size(Ct,1);
    % ground truth source cluster index is a cell array of size n_source_cluster
    % each cell contains array of indices of source in each cluster (for now, each cluster is each roi)
    F_roi_true = calgcss_block(At,Ct,Wt,Nt,[],groundtruth_source_cluster_index{mm});
    

    
    C_ind = sys_est(tt).C_out.ind_chosen_Lpq.bic;% use bic score
    A = sys_est.sys.A;
    C = sys_est.C_out.C_Lpq_CLS(:,:,C_ind);
    E = sys_est.C_out.E_Lpq_CLS(:,:,C_ind);
    W = sys_est.C_out.W_Lpq_CLS(:,:,C_ind);
    [N,V,obj] = noisecovest(E,L,'homo'); % require cvx
    [Ftmp,~,~] = calgcss(A,C,W,N,[]);
    F{mm,1} = sparse(Ftmp);
    F_nz_ind{mm,1} = find(Ftmp);
    Ftmp_roi = calgcss_block(A,C,W,N,[],estimate_source_cluster_index{mm});
    F_roi{mm,1} = sparse(Ftmp_roi);
    F_roi_nz_ind{mm,1} = find(Ftmp_roi);
    
    F_true{mm,1} = model.F0;
    F_true_ind{mm,1} = find(model.F0);
    F_true_roi{mm,1} = F_roi_true;
end
M.F_true{1} = F_true;
M.F_true_roi{1} = F_true_roi;
M.F_true_ind{1} = F_true_ind;
M.F_nz_ind{1} = F_nz_ind;
M.F{1} = F;
M.F_roi{1} = F_roi;
M.F_roi_nz_ind{1} = F_roi_nz_ind;

%% Step (2): Map the calculated GC matrices into the reference of 8 ROIs 
for mm=1:n_model
    postfix = ['_act_',int2str(n_active_source),'_deep_',int2str(n_deep),'_',int2str(mm)];
    load([datapath,'model_',postfix])
    
    F_true_roi_tmp = [M.F_true_roi{1}{mm} zeros(4,4);zeros(4,8)];
    F_true_roi_tmp = F_true_roi_tmp.*(F_true_roi_tmp>1e-6);
    E.F_true_roi{1}{mm} = F_true_roi_tmp(perm_roi_ind{mm},perm_roi_ind{mm});
    E.F_true_roi_ind{1}{mm} = find(M.F_true_roi{1}{mm});
    F_roi_tmp = M.F_roi{1}{mm};
    F_roi_tmp = F_roi_tmp.*(F_roi_tmp>1e-6); 
    % The value 1e-6 is used to threshold numerical error in zero causality
    % estimation. The magnitude of ground-truth GC is in order of 10^0
    E.F_roi{1}{mm} = F_roi_tmp(perm_roi_ind{mm},perm_roi_ind{mm});
    E.F_roi_nz_ind{1}{mm} = find(E.F_roi{1}{mm});
    
    F_true_tmp = [M.F_true{1}{mm} zeros(50,50);zeros(50,100)];
    E.F_true{1}{mm} = F_true_tmp(perm_node_ind{mm},perm_node_ind{mm});
    E.F{1}{mm} = M.F{1}{mm}(perm_node_ind{mm},perm_node_ind{mm});
    E.F_nz_ind{1}{mm} = {find(E.F{1}{mm})};
    E.perm_node_ind{mm,1} = perm_node_ind{mm};
    
    
    
end



save([outpath,'simeeg_sourcerecon_Fresult','M']);


