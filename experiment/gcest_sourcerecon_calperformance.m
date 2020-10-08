%   This function calculates performance indices for each threshold and
%   export to 'outpath/performance.mat'. The output file is in the
%   following format;
%
%   performance_F_[roi or MN or LCMV or sLORETA]
%       {between 1 to n_model}      =   performance of each model
%           .TPR
%           .FPR
%           .ACC
%			.F1
%			.TPR_avg
%           .FPR_avg
%           .ACC_avg
%			.F1_avg
%			.TPR_std
%           .FPR_std
%           .ACC_std
%			.F1_std
%           
%   INPUT   inpath_mvgc  =   input file path for F obtained from MVGC method
%			inpath_roi   =   input file path for F obtained from our proposed method
%   OUTPUT  outpath 	 =   output file path ('outpath/performance')
%
%   Written by ANAWAT NARTKULPAT, 2020
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
%========================================================================== 

clear

inpath_mvgc = './saved_experiment_results/source_recon_brainstorm/gcest_sourcerecon_F_mvgc';
inpath_roi = './saved_experiment_results/source_recon_brainstorm/F_result_L_augment_permuted';
outpath = './saved_results_toplot/gcest_sourcerecon_performance';

M_mvgc = load(inpath_mvgc);
M_roi = load(inpath_roi);
% load ground truth non-zero indices
F_true = M_roi.M.F_true_roi{1};
F_true_ind = cellfun(@(x) find(x > 0),F_true,'UniformOutput',false);
F_true_ind = F_true_ind';

% call size of matrix and data cell
[dim_F,~] = size(M_roi.M.F_roi{1,1,1}{1,1});
[n_ratio_source,n_active,n_electrode] = size(M_roi.M.F_roi_nz_ind);
[n_model,n_exp] = size(M_roi.M.F_nz_ind{1,1,1});

% predefine performance and non-zero index cell and other temporary
% variable.
% performance = cell(n_ratio_source,n_active,n_electrode);
% ind_nz = cell(n_model,n_exp);
ind_diag = [1:dim_F+1:dim_F^2]';

% create threshold level from 0 to 0.1 with logspace between 10^-10 to 10^-1 
% and the rest are evenly spaced with linspace
Num_thresh = 1000;
Nf = floor(0.5*Num_thresh);
% threshold_level = [0,logspace(-10,-1,Nf-1),linspace(0.1,1,Num_thresh-Nf)];
threshold_level = [0,logspace(-6,0,Num_thresh-1)];


for m = 1:Num_thresh
    % threshold every estimated F with (0-1)*(max(F)-min(F))
    ind_nz = cellfun(@(x) setdiff(find(x >= threshold_level(m)*(max(x(:))-min(x(x > 0)))),ind_diag),M_roi.M.F_roi{1},'UniformOutput',false);
    % obtain performance indices for each model
    perf_temp = compare_F(ind_nz',F_true_ind,dim_F);
    performance_F_roi.TPR(:,m) = perf_temp.TPR(:);
    performance_F_roi.FPR(:,m) = perf_temp.FPR(:);
    performance_F_roi.ACC(:,m) = perf_temp.ACC(:);
    performance_F_roi.TPR_avg(m) = mean(perf_temp.TPR(:));
    performance_F_roi.FPR_avg(m) = mean(perf_temp.FPR(:));
    performance_F_roi.ACC_avg(m) = mean(perf_temp.ACC(:));
    performance_F_roi.TPR_std(m) = std(perf_temp.TPR(:));
    performance_F_roi.FPR_std(m) = std(perf_temp.FPR(:));
    performance_F_roi.ACC_std(m) = std(perf_temp.ACC(:));
    performance_F_roi.F1(:,m) = perf_temp.F1(:);
    
    % check for zero and non-zero
    F = M_mvgc.M.F_sparse_MN(:,m);
    ind_nz = cellfun(@(x) setdiff(find(x > 0),ind_diag),F,'UniformOutput',false);
    perf_temp = compare_F(ind_nz,F_true_ind,dim_F);
    performance_F_MN.TPR(:,m) = perf_temp.TPR(:);
    performance_F_MN.FPR(:,m) = perf_temp.FPR(:);
    performance_F_MN.ACC(:,m) = perf_temp.ACC(:);
    performance_F_MN.TPR_avg(m) = mean(perf_temp.TPR(:));
    performance_F_MN.FPR_avg(m) = mean(perf_temp.FPR(:));
    performance_F_MN.ACC_avg(m) = mean(perf_temp.ACC(:));
    performance_F_MN.TPR_std(m) = std(perf_temp.TPR(:));
    performance_F_MN.FPR_std(m) = std(perf_temp.FPR(:));
    performance_F_MN.ACC_std(m) = std(perf_temp.ACC(:));
    performance_F_MN.F1(:,m) = perf_temp.F1(:);

    F = M_mvgc.M.F_sparse_LCMV(:,m);
    ind_nz = cellfun(@(x) setdiff(find(x > 0),ind_diag),F,'UniformOutput',false);
    perf_temp = compare_F(ind_nz,F_true_ind,dim_F);
    performance_F_LCMV.TPR(:,m) = perf_temp.TPR(:);
    performance_F_LCMV.FPR(:,m) = perf_temp.FPR(:);
    performance_F_LCMV.ACC(:,m) = perf_temp.ACC(:);
    performance_F_LCMV.TPR_avg(m) = mean(perf_temp.TPR(:));
    performance_F_LCMV.FPR_avg(m) = mean(perf_temp.FPR(:));
    performance_F_LCMV.ACC_avg(m) = mean(perf_temp.ACC(:));
    performance_F_LCMV.TPR_std(m) = std(perf_temp.TPR(:));
    performance_F_LCMV.FPR_std(m) = std(perf_temp.FPR(:));
    performance_F_LCMV.ACC_std(m) = std(perf_temp.ACC(:));
    performance_F_LCMV.F1(:,m) = perf_temp.F1(:);

    F = M_mvgc.M.F_sparse_sLORETA(:,m);
    ind_nz = cellfun(@(x) setdiff(find(x > 0),ind_diag),F,'UniformOutput',false);
    perf_temp = compare_F(ind_nz,F_true_ind,dim_F);
    performance_F_sLORETA.TPR(:,m) = perf_temp.TPR(:);
    performance_F_sLORETA.FPR(:,m) = perf_temp.FPR(:);
    performance_F_sLORETA.ACC(:,m) = perf_temp.ACC(:);
    performance_F_sLORETA.TPR_avg(m) = mean(perf_temp.TPR(:));
    performance_F_sLORETA.FPR_avg(m) = mean(perf_temp.FPR(:));
    performance_F_sLORETA.ACC_avg(m) = mean(perf_temp.ACC(:));
    performance_F_sLORETA.TPR_std(m) = std(perf_temp.TPR(:));
    performance_F_sLORETA.FPR_std(m) = std(perf_temp.FPR(:));
    performance_F_sLORETA.ACC_std(m) = std(perf_temp.ACC(:));
    performance_F_sLORETA.F1(:,m) = perf_temp.F1(:);
end
       
save('./saved_results_toplot/gcest_sourcerecon_performance','performance_F_roi',...
    'performance_F_MN',...
    'performance_F_LCMV',...
    'performance_F_sLORETA');
