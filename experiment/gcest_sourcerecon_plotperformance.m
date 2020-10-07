%   This function creates 7 plots relating to performances
%
%   INPUT   inpath_F    =   input estimated F file path
%           inpath_perf =   input performance file path
%   OUTPUT  outpath     =   output file path ('outpath/figx')
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

clear;
outpath = 'Experiment/figure';
load('Experiment/data/gcest_sourcerecon_performance');
load('Experiment/data/miscdata');
M_mvgc = load('Experiment/data/gcest_sourcerecon_F_mvgc');
M_roi = load('Experiment/data/F_result_L_augment_permuted');

%==========================================================================

% plot figure 1: example of ground truth F with estimated F from varying number of electrodes
figure(1);
[n_model,n_realization] = size(M_roi.M.F{1,1,1});
ground_truth_ind = 24;
thresh_chosen_ind = 2;
alpha_chosen = 10;
t = tiledlayout(2,3);
ax1 = nexttile; imagesc(M_roi.M.F_true_roi{1,1}{ground_truth_ind}); colormap(flipud(hot)); 
axis square; title('Ground-truth ROI GC');
set(gca,'xtick',1:8,'xticklabel',rois,'ytick',1:8,'yticklabel',rois);
xtickangle(ax1,-45);

TPR = performance_F_roi.TPR(ground_truth_ind,thresh_chosen_ind);
FPR = performance_F_roi.FPR(ground_truth_ind,thresh_chosen_ind);
ACC = performance_F_roi.ACC(ground_truth_ind,thresh_chosen_ind);
nexttile; imagesc(M_roi.M.F_roi{1,1,1}{ground_truth_ind}); colormap(flipud(hot)); axis square;
title('Proposed method');% caxis([0 1]);
xlabel(sprintf('(%.3f, %.3f, %.3f)',TPR,FPR,ACC),'Interpreter','None');
set(gca,'xticklabel',[],'ytick',1:8,'yticklabel',rois);

TPR = performance_F_MN.TPR(ground_truth_ind,alpha_chosen);
FPR = performance_F_MN.FPR(ground_truth_ind,alpha_chosen);
ACC = performance_F_MN.ACC(ground_truth_ind,alpha_chosen);
temp = nexttile;
nexttile; imagesc(M_mvgc.M.F_sparse_MN{ground_truth_ind,alpha_chosen}); colormap(flipud(hot)); axis square;
title('WMNE (VAR)');% caxis([0 1]);
xlabel(sprintf('(%.3f, %.3f, %.3f)',TPR,FPR,ACC),'Interpreter','None');
set(gca,'xticklabel',[],'ytick',1:8,'yticklabel',rois);

TPR = performance_F_LCMV.TPR(ground_truth_ind,alpha_chosen);
FPR = performance_F_LCMV.FPR(ground_truth_ind,alpha_chosen);
ACC = performance_F_LCMV.ACC(ground_truth_ind,alpha_chosen);
nexttile; imagesc(M_mvgc.M.F_sparse_LCMV{ground_truth_ind,alpha_chosen}); colormap(flipud(hot)); axis square;
title('LCMV (VAR)');% caxis([0 1]);
xlabel(sprintf('(%.3f, %.3f, %.3f)',TPR,FPR,ACC),'Interpreter','None');
set(gca,'xticklabel',[],'ytick',1:8,'yticklabel',rois);

TPR = performance_F_sLORETA.TPR(ground_truth_ind,alpha_chosen);
FPR = performance_F_sLORETA.FPR(ground_truth_ind,alpha_chosen);
ACC = performance_F_sLORETA.ACC(ground_truth_ind,alpha_chosen);
nexttile; imagesc(M_mvgc.M.F_sparse_sLORETA{ground_truth_ind,alpha_chosen}); colormap(flipud(hot)); axis square;
title('sLORETA (VAR)');% caxis([0 1]);
xlabel(sprintf('(%.3f, %.3f, %.3f)',TPR,FPR,ACC),'Interpreter','None');
set(gca,'xticklabel',[],'ytick',1:8,'yticklabel',rois);

delete(temp);

t.Padding = 'none';
t.TileSpacing = 'none';

saveas(gcf,[outpath,'/gcest_sourcerecon_roi_fhat24']); hold off;
print([outpath,'/gcest_sourcerecon_roi_fhat24'], '-depsc')
