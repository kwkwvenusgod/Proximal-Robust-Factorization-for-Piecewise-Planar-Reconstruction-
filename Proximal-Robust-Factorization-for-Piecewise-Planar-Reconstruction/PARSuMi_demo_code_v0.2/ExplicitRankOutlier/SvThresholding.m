function W_threshold = SvThresholding(W,gamma)
%SVT Summary of this function goes here
%   Detailed explanation goes here

[U,S,V] = svd(W,'econ');
s = diag(S);
s_threshold = max(s-gamma,0);
S_threshold = diag(s_threshold);
W_threshold = U*S_threshold*V';

end

