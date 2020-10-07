%   This script calculates performance indices for each threshold and
%   export to 'outpath/performance.mat'. The output file is in the
%   following format;
%
%   performance{rois_deep,m_active,n_electrodes}
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
%           .TPR_q  (interquartiles)
%           .FPR_q (interquartiles)
%           .ACC_q (interquartiles)
%			.F1_q (interquartiles)
%           
%           
%   INPUT   inpath  =   input file path
%   OUTPUT  outpath =   output file path ('outpath/gcest_performance')
%
%   Written by ANAWAT NARTKULPAT
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
inpath = './saved_experiment_results/simeeg_varyfactor/simeeg_varyfactor_Fresult';
outpath = './saved_results_toplot/gcest_performance';

load(inpath);
% load ground truth non-zero indices
F_true_ind = M.F_true_ind;

% call size of matrix and data cell
[dim_F,~] = size(M.F{1,1,1}{1,1});
[n_ratio_source,n_active,n_electrode] = size(M.F_nz_ind);
[n_model,n_exp] = size(M.F_nz_ind{1,1,1});

% predefine performance and non-zero index cell and other temporary variable.
performance = cell(n_ratio_source,n_active,n_electrode);
ind_nz = cell(n_model,n_exp);
ind_diag = [1:dim_F+1:dim_F^2]';

% create threshold level from 0 to 0.1 with logspace between 10^-10 to 10^-1 
% and the rest are evenly spaced with linspace
Num_thresh = 1000;
Nf = floor(0.5*Num_thresh);
% threshold_level = [0,logspace(-10,-1,Nf-1),linspace(0.1,1,Num_thresh-Nf)];
threshold_level = [0,logspace(-6,0,Num_thresh-1)];

for i = 1:n_ratio_source
    disp(i);
    for j = 1:n_active
        disp(j);
        for k = 1:n_electrode
            disp(k);
            for m = 1:Num_thresh
                cnt = 0;
                % threshold every estimated F with (0-1)*(max(F)-min(F))
                ind_nz = cellfun(@(x) setdiff(find(x >= threshold_level(m)*(max(x(:))-min(x(x > 0)))),ind_diag),M.F{i,j,k},'UniformOutput',false);
                % obtain performance indices for each model
                perf_temp = compare_F(ind_nz,F_true_ind{i,j},dim_F);
                performance{i,j,k}.TPR(:,m) = perf_temp.TPR(:);
                performance{i,j,k}.FPR(:,m) = perf_temp.FPR(:);
                performance{i,j,k}.ACC(:,m) = perf_temp.ACC(:);
                performance{i,j,k}.F1(:,m) = perf_temp.F1(:);
                performance{i,j,k}.TPR_avg(m) = mean(perf_temp.TPR(:));
                performance{i,j,k}.FPR_avg(m) = mean(perf_temp.FPR(:));
                performance{i,j,k}.ACC_avg(m) = mean(perf_temp.ACC(:));
                performance{i,j,k}.F1_avg(m) = mean(perf_temp.F1(:));
                performance{i,j,k}.TPR_std(m) = std(perf_temp.TPR(:));
                performance{i,j,k}.FPR_std(m) = std(perf_temp.FPR(:));
                performance{i,j,k}.ACC_std(m) = std(perf_temp.ACC(:));
                performance{i,j,k}.F1_std(m) = std(perf_temp.F1(:));

                % calculate the 1st, 2nd and 3rd qurtile of TPR, FPR and ACC
                dat =  performance{i,j,k}.TPR(:,m);
                q123 = [prctile(dat,25); prctile(dat,50); prctile(dat,75)]; % 3 x 1
                performance{i,j,k}.TPR_q(:,m) = q123;
                dat =  performance{i,j,k}.FPR(:,m);
                q123 = [prctile(dat,25); prctile(dat,50); prctile(dat,75)];
                performance{i,j,k}.FPR_q(:,m) = q123;
                dat =  performance{i,j,k}.ACC(:,m);
                q123 = [prctile(dat,25); prctile(dat,50); prctile(dat,75)]; 
                performance{i,j,k}.ACC_q(:,m) = q123;
                dat =  performance{i,j,k}.F1(:,m);
                q123 = [prctile(dat,25); prctile(dat,50); prctile(dat,75)];
                performance{i,j,k}.F1_q(:,m) = q123;
            end
        end
    end
end
save(outpath,'gcest_performance');
