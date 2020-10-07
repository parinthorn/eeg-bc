%   This script estimate VAR model and calculate Granger causality of
%   reconstucted source data from Brainstorm. The causality is then 
%   statistically tested to classify zero or non-zero causality
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

clc
clear

inpath = './saved_experiment_results/simeeg_sourcerecon/brainstorm_avg_recon_source/'; 
outpath = './saved_experiment_results/source_recon_brainstorm/gcest_sourcerecon_F_mvgc';

n_model = 100;    % number of data files

%% Parameters

nx = 1; % number of target ("to") variables
ny = 1; % number of source ("from") variables
nz = 0; % number of conditioning variables (default: 0)

morder = 10; % maximum model order for model order estimation
regmode = 'LWR'; % VAR model estimation regression mode 'LWR'(default) or 'OLS'
acmaxlags = 1000;  % maximum auto covariance lags to calculate

tstat     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.01;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction

for i = 1:n_model
    tic;
    load([inpath,num2str(i)]);
    x_MN = x.MN;
    x_LCMV = x.LCMV;
    x_sLORETA = x.sLORETA;
    
    %% Calculate
    [numvar,numobserve,ntrials] = size(x_MN);
    ind_diag = 1:numvar+1:numvar^2;
    [aic_MN,bic_MN,moaic_MN,mobic_MN] = tsdata_to_infocrit(x_MN,morder); % calculate aic&bic from time series return best model order as moaic and mobic
    [aic_LCMV,bic_LCMV,moaic_LCMV,mobic_LCMV] = tsdata_to_infocrit(x_LCMV,morder);
    [aic_sLORETA,bic_sLORETA,moaic_sLORETA,mobic_sLORETA] = tsdata_to_infocrit(x_sLORETA,morder);
    
    [A_MN,SIG_MN,E_MN] = tsdata_to_var(x_MN,mobic_MN,regmode); % fit VAR model to time series, return A:VAR coefficients matrix
    [A_LCMV,SIG_LCMV,E_LCMV] = tsdata_to_var(x_LCMV,mobic_LCMV,regmode);
    [A_sLORETA,SIG_sLORETA,E_sLORETA] = tsdata_to_var(x_sLORETA,mobic_sLORETA,regmode);
    % SIG:residuals covariance matrix, E:residuals time series

    [G_MN,info_MN] = var_to_autocov(A_MN,SIG_MN,acmaxlags); % calculate G:autocovariance sequence for a VAR model
    [G_LCMV,info_LCMV] = var_to_autocov(A_LCMV,SIG_LCMV,acmaxlags);
    [G_sLORETA,info_sLORETA] = var_to_autocov(A_sLORETA,SIG_sLORETA,acmaxlags);

    M.F_MN{i,1} = autocov_to_pwcgc(G_MN); % calculate pairwise-conditional time-domain multivariate Granger causalities
    M.F_LCMV{i,1} = autocov_to_pwcgc(G_LCMV);
    M.F_sLORETA{i,1} = autocov_to_pwcgc(G_sLORETA);
    
    pval_MN{i,1} = mvgc_pval(M.F_MN{i,1},mobic_MN,numobserve,ntrials,nx,ny,nz,tstat); % p-values for sample MVGC based on theoretical asymptotic null distribution
    pval_LCMV{i,1} = mvgc_pval(M.F_LCMV{i,1},mobic_LCMV,numobserve,ntrials,nx,ny,nz,tstat);
    pval_sLORETA{i,1} = mvgc_pval(M.F_sLORETA{i,1},mobic_sLORETA,numobserve,ntrials,nx,ny,nz,tstat);
    
    alpha_range = linspace(0,1,1000);
    for j = 1:length(alpha_range)
        sig_MN = significance(pval_MN{i,1},alpha_range(j),mhtc); % statistical significance adjusted for multiple hypotheses
        sig_LCMV = significance(pval_LCMV{i,1},alpha_range(j),mhtc);
        sig_sLORETA = significance(pval_sLORETA{i,1},alpha_range(j),mhtc);
        M.F_sparse_MN{i,j} = M.F_MN{i,1}.*sig_MN;
        M.F_sparse_LCMV{i,j} = M.F_LCMV{i,1}.*sig_LCMV;
        M.F_sparse_sLORETA{i,j} = M.F_sLORETA{i,1}.*sig_sLORETA;
        
        M.F_sparse_MN{i,j}(ind_diag) = 0;
        M.F_sparse_LCMV{i,j}(ind_diag) = 0;
        M.F_sparse_sLORETA{i,j}(ind_diag) = 0;
    end
    
    M.F_MN{i,1}(ind_diag) = 0;
    M.F_LCMV{i,1}(ind_diag) = 0;
    M.F_sLORETA{i,1}(ind_diag) = 0;
    
    M.mobic_MN{i,1} = mobic_MN;
    M.mobic_LCMV{i,1} = mobic_LCMV;
    M.mobic_sLORETA{i,1} = mobic_sLORETA;
    toc;
end
save(outpath,'M');
