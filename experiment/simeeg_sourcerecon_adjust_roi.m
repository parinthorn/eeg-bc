%% This script prepare the indices of ground-truth/estimate sources 
% The EEG data is generated under the setting: 50% deep source, 20% active source, 61 No. of electrodes
% 
% Input: Ground-truth models
% Output: Ground-truth/Estimate sources indices (in clusters) with the reference of 4 ROIs and
% Permutation indices (to convert the indices to 8 ROIs)
% 
% Pre-knowledge
% - All sources (both active and inactive) in each generated model are chosen randomly in 4 ROIs.
% - When we generate Ltilde (for estimation), 
% the sources in consideration lie in the first ground-truth 4 ROIs and the remainning 4 ROIS.
% - In our experiment, all sources in consideration is twice the number of sources in ground-truth (mtilde = 2m).
% - Therefore, the indices of estimated sources assigned in the remaining 4 ROIs will be shifted by m
% 
% Step:
% (1) Read model (all sources in each model are chosen randomly in 4 ROIs)


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

datapath = './input_data/source_recon_brainstorm/eeg_simulated_data_deep50percent_active20percent/';

inpath = './saved_experiment_results/simeeg_sourcerecon/';
outpath = './saved_experiment_results/simeeg_sourcerecon/';

n_model = 100;
n_deep = 2; % read only the case that there are 2 deep source ROIs
n_active_source = 10; % read the case that there are 10 active sources (from 50 sources)

load([datapath,'model_act_10_deep_2_',int2str(1)])
 
%%  Construct ground-truth/estimate source cluster indices in the reference of 4 ROIs (differ from models to models)

    % ground truth/estimate source cluster index is a cell array of size n_source_cluster
    % each cell contains array of indices of source in each cluster (for now, each cluster is each roi)
    
    groundtruth_source_cluster_index = cell(model.PARAMETER.n_source_cluster,1); % correspond to 4 ROIs in our exp
    estimate_source_cluster_index = cell(2*model.PARAMETER.n_source_cluster,1); % correspond to 8 ROIs in our exp
    permind = cell(n_model,1);
    
    % stack all indices from all models
    sourceind.perm_node_ind = cell(n_model,1); % permutation mapping (node-wise)
    sourceind.true_source_cluster_ind = cell(n_model,1); % true source cluster index (4 ROIs)
    sourceind.est_source_cluster_ind = cell(n_model,1); % estimation source cluster index (4 ROIs)
    
for mm=1:n_model    
%     postfix = 
    load([datapath,'model_act_10_deep_2_',int2str(mm)])
    for ii=1:model.PARAMETER.n_source_cluster-1
        groundtruth_source_cluster_index{ii} = [model.PARAMETER.ind_cluster_source(ii):1:model.PARAMETER.ind_cluster_source(ii+1)-1];
    end
      
    n_all_source = model.PARAMETER.m; % parameter 'm' in the paper
    
    % in our estimation, the number of all sources in condiration is that
    % mtilde = 2*m. The first 'm' sources stay in the same ROIs as
    % ground-truth sources, the remaining 'm' sources are assigned to the
    % remaining 4 ROIs. Therefore, when we construct Ltilde, we shift the index of additional sources to
    % ground-truth by m
    groundtruth_source_cluster_index{end} = model.PARAMETER.ind_cluster_source(end):1:n_all_source;

    tmpf = @(x) x+ n_all_source;
        
    estimate_source_cluster_index = [groundtruth_source_cluster_index; ...
        cellfun(tmpf,groundtruth_source_cluster_index,'UniformOutput',false)];
    
    [~,permind{mm,1}]=sort([unique(model.truth.in_roi),unique(model.truth.augment_data.in_roi)]);
    
    permute_source_cluster_index = (estimate_source_cluster_index(permind{mm}));
    permute_source_cluster_index= [permute_source_cluster_index{:}];
    
    sourceind.perm_node_ind{mm} = permute_source_cluster_index; 
    sourceind.true_source_cluster_ind{mm} = groundtruth_source_cluster_index; 
    sourceind.est_source_cluster_ind{mm} = estimate_source_cluster_index;
    

end
    sourceind.perm_roi_ind = permind;
save([outpath,'simeeg_sourcerecon_source_perm_ind'],'sourceind');
