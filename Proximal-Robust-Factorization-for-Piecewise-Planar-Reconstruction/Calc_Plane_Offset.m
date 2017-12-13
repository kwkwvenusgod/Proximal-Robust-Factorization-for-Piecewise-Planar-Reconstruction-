function [ plane_Z0 ] = Calc_Plane_Offset(structure_matrix,focal_length )
%CAL_DEPTH Summary of this function goes here
%   Detailed explanation goes here
plane_Z0_inv=structure_matrix(1,1:6:size(structure_matrix,2))/focal_length;
plane_Z0=-(1./plane_Z0_inv);
plane_Z0=plane_Z0';
end