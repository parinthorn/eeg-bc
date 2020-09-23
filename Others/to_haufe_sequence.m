function [F,F_block] = to_haufe_sequence(F_in,F_block_in,permute_roi_index,m,mtilde,n_source_cluster,ind_cluster_source)
% creating source cluster index in cell format
% order will be {'LPS','LPI','LAS','LAI','RPS','RPI','RAS','RAI'}
source_cluster_index = cell(n_source_cluster,1);
for ii=1:n_source_cluster-1
    source_cluster_index{ii} = [ind_cluster_source(ii):1:ind_cluster_source(ii+1)-1];
end
source_cluster_index{end} = ind_cluster_source(end):1:50;
tmpf = @(x) x+m;
source_cluster_index = [source_cluster_index; ...
    cellfun(tmpf,source_cluster_index, ...
    'UniformOutput',false)];
permute_source_index = (source_cluster_index(permute_roi_index));
if size(F_in,1)==m
    F = [F_in zeros(mtilde-m,mtilde-m);zeros(mtilde-m,mtilde)];
    F = F([permute_source_index{:}],[permute_source_index{:}]);
    F_block = [F_block_in zeros(4,4);zeros(4,8)]; % hard code part
    F_block = F_block(permute_roi_index,permute_roi_index);
else
    F = F_in(permute_source_index{:},permute_source_index{:});
    F_block = F_block_in(permute_roi_index,permute_roi_index);
end
end
