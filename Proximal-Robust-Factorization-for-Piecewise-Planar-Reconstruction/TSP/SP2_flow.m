function [utotal,vtotal] = SP2_flow(sp_labels)

xdim = size(sp_labels,1);
ydim = size(sp_labels,2);

K = double(max(sp_labels(:)))+1;

[means1] = SP_mean_locations(K, uint32(sp_labels(:,:,1)));
[means2] = SP_mean_locations(K, uint32(sp_labels(:,:,2)));

mask = means1(:,1)>=0 & means2(:,1)>=0;

v = means2(mask,1) - means1(mask,1);
u = means2(mask,2) - means1(mask,2);

[indices, ~] = populate_indices(double(K), uint64(sp_labels(:,:,2)));
indices(~mask) = [];

K = numel(indices);

[utotal, vtotal] = flow_total(K,xdim,ydim,indices,u,v);