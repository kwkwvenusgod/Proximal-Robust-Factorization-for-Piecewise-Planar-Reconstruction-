function [ RMSE ] = OracleBound( sigma,m,n, r,p )
%ORACLEBOUND Summary of this function goes here
%   This function calculated the oracle bound as if prior knowledge of the
%   linear subspace.
% sigma is the standard deviation 

RMSE=sigma*sqrt((m+n-r)*r./p);

end
