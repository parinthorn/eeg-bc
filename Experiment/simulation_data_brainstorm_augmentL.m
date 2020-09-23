%% Simulation data experiment: 50% deep source
%%                             20% active source
%%                             61 No. of electrodes
%% This script generate source amplitude for generating Leadfield matrix
%% Written by: PARINTHORN MANOMAISAOWAPAK
clc
if exist('sa','var')
    clearvars -except sa
else
    clear
    load('G:\Shared drives\MASTER_DRIVE\Journal\DATA_GENERATION\BBCB_code\data\sa')
end
load('G:\Shared drives\MASTER_DRIVE\Journal\DATA_GENERATION\BBCB_code\data\miscdata')

for ii=1:100
    tic
load(['G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\model__act_10_deep_2_',int2str(ii)])
load(['G:\Shared drives\MASTER_DRIVE\Journal\DATA\experiment_paper_L_augment\eegdata__act_10_deep_2_',int2str(ii)])
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
