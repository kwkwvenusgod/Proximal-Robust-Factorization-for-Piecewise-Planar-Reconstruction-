function [ normal ] = Calc_Normal( position_matrix,f )
%CALC_NORMAL Summary of this function goes here
%   Detailed explanation goes here

Z0=Calc_Plane_Offset(position_matrix,f);
ratio_xz=position_matrix(1,3:6:size(position_matrix,2)).*Z0';
ratio_yz=position_matrix(1,5:6:size(position_matrix,2)).*Z0';
normalz=-1./sqrt(ratio_xz.*ratio_xz+ratio_yz.*ratio_yz+ones(size(ratio_yz)));
normaly=-ratio_yz.*normalz;
normalx=-ratio_xz.*normalz;
normal=[normalx',normaly',normalz'];
end

