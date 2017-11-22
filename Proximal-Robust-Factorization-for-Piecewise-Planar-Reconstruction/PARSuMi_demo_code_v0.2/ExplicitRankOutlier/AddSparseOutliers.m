function [Wcap_outlier,E_gt,E_mask] = AddSparseOutliers(Wcap,sample_mask,fraction)
%ADDSPARSEOUTLIERS Summary of this function goes here
%   Detailed explanation goes here

amplitude = 5.0;

[m,n] = size(sample_mask);
sample_index = find(sample_mask(:));
permuted_sample_index = randperm(length(sample_index));
n_outlier = round(length(sample_index)*fraction);
outlier_index = sample_index(permuted_sample_index(1:n_outlier));
E_mask_vec = zeros(m*n,1);
E_mask_vec(outlier_index) = 1;
E_mask = reshape(E_mask_vec,m,[]);
G = randn(m,n)* amplitude;
E_gt = E_mask .* G ;
Wcap_outlier = Wcap + E_gt;

end

