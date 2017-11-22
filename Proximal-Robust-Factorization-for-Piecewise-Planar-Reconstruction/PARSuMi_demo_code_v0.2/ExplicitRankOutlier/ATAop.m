function M = ATAop(M,sample_mask)
%ATAOPERATOR Summary of this function goes here
%   Detailed explanation goes here

M(~sample_mask) = 0;

end

