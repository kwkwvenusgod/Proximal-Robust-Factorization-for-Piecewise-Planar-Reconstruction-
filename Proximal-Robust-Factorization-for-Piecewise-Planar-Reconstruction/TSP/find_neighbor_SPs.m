function [neighbors] = find_neighbor_SPs(label, k)

mask = (label==k);
mask = imdilate(mask, [0 1 0; 1 1 1; 0 1 0]);
neighbors = label(mask);

neighbors = unique(neighbors)';
neighbors = neighbors(neighbors~=k);