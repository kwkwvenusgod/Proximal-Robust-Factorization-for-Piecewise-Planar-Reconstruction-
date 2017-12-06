function [ normal ] = Calc_Normal( position_matrix,f )
%CALC_NORMAL Summary of this function goes here
%   Detailed explanation goes here
depth_inv=position_matrix(1,1:6:size(position_matrix,2))/f;
depth=-(1./depth_inv);
ratio_xz=position_matrix(1,3:6:size(position_matrix,2)).*depth;
ratio_yz=position_matrix(1,5:6:size(position_matrix,2)).*depth;
normalz=-1./sqrt(ratio_xz.*ratio_xz+ratio_yz.*ratio_yz+ones(size(ratio_yz)));
normaly=-ratio_yz.*normalz;
normalx=-ratio_xz.*normalz;
normal=[normalx',normaly',normalz'];
end

