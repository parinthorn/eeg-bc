%% This script computes source amplitude distribution for each generated
% model in SIMEEG_SOURCERECON experiment
% Step 
% (1) We have all sources in considerations (from selected mtilde) for each model
% (2) From the location of each source, the amplitude distribution follows the Gaussian function
%
%  Written by Parinthorn Manomaisaowapak
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
addpath('./input_data/nyhead_model/');    
inpath = './input_data/source_recon_brainstorm/';
outpath = './saved_experiment_results/source_recon_brainstorm/';
if exist('sa','var')
    clearvars -except sa
else
    clear
    load('sa');
end
load('miscdata');

%% 
n_model = 100;
s_amp_augmented = cell{n_model,1};
for ii=1:n_model
    tic
load([inpath,'model_act_10_deep_2_',int2str(ii)])
load([inpath,'eegdata_act_10_deep_2_',int2str(ii)])
sigmas = model.truth.sigmas;
n_source_augment = 50;
cortex2K = [];
cortex2K.vc = sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K, :);
cortex2K.tri = sa.cortex2K.tri;

in_centers = model.truth.augment_data.in_centers;
in_roi = model.truth.augment_data.in_roi;

for i = 1:n_source_augment
[~, source_amp(:, i)] = graphrbf(cortex2K, sigmas(i), in_centers(i));
end
for i = 1:n_source_augment
    % index that is not in INDS_ROI_OUTER_2K
    ind_tmp = setdiff(1:size(source_amp,1), inds_roi_outer_2K{ in_roi(i) });
    source_amp( ind_tmp, i) = 0;
end
for k=1:n_source_augment
  TMP = norm(source_amp(:,k)); % since the source of some channel can be entirely zero
  if TMP > 1e-10
    source_amp(:, k) = source_amp(:, k) ./ TMP;
  end
end

full_source_amp = [model.truth.source_amp source_amp];
s_amp_augmented{ii,1} = sparse(full_source_amp);
toc
disp(ii)
end
% save([outpath,'source_amplitude','s_amp_augmented']);

