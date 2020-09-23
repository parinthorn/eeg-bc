function [L_augment,truth_augment] = upsample_L(n_source,ind_roi_cluster,ind_cluster_source,n_source_cluster,sa,rois_in,sigmas,inds_roi_inner_2K,inds_roi_outer_2K)
% assign ROI to each source. Sources in the same cluster are assigned to the same ROI
ind_roi_cluster_augment = setdiff(1:8,ind_roi_cluster);
truth_augment.in_roi = zeros(1,n_source);
for i = 1:n_source_cluster-1
truth_augment.in_roi(ind_cluster_source(i):ind_cluster_source(i+1)-1) = ind_roi_cluster_augment(i);
end
truth_augment.in_roi(ind_cluster_source(n_source_cluster):n_source) = ind_roi_cluster_augment(n_source_cluster);

truth_augment.rois = rois_in(truth_augment.in_roi);

% sample source extents

% sample source centers within rois
for i = 1:n_source
  truth_augment.in_centers(i) = inds_roi_inner_2K{truth_augment.in_roi(i)}(ceil(length(inds_roi_inner_2K{truth_augment.in_roi(i)})*rand(1)));
end
truth_augment.centers = sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K(truth_augment.in_centers), :);

% calculate source amplitude distribution
cortex2K = [];
cortex2K.vc = sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K, :);
cortex2K.tri = sa.cortex2K.tri;

% generate Gaussian distributions on cortical manifold % loop 1:n_sources
for i = 1:n_source
[~, truth_augment.source_amp(:, i)] = graphrbf(cortex2K, sigmas(i), truth_augment.in_centers(i));
end

% set activity outside of source octant to zero %%%==========================================
for i = 1:n_source
    % index that is not in INDS_ROI_OUTER_2K
    ind_tmp = setdiff(1:size(truth_augment.source_amp,1), inds_roi_outer_2K{ truth_augment.in_roi(i) });
    truth_augment.source_amp( ind_tmp, i) = 0;
end

for k=1:n_source
  TMP = norm(truth_augment.source_amp(:,k)); % since the source of some channel can be entirely zero
  if TMP > 1e-10
    truth_augment.source_amp(:, k) = truth_augment.source_amp(:, k) ./ TMP;
  end
end

% obtain the lead field matrix
L_augment = sa.cortex75K.EEG_V_fem_normal(:, sa.cortex2K.in_from_cortex75K)*truth_augment.source_amp; % (n_electrode x 2K) x (2K x n_source)

end
