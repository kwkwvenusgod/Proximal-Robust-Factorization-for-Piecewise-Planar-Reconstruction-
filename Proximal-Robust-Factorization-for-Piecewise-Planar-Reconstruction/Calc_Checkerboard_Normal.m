function [ checkerboard_normal ] = Calc_Checkerboard_Normal( normal_map,checkerboard_mask, label_map_consistency, label_set )
%CALC_CHECKERBOARD_NORMAL Summary of this function goes here
%   Detailed explanation goes here
graymask = rgb2gray(checkerboard_mask);
ind = find(graymask==255);
normal_label = label_map_consistency(ind);
normal_label = unique(normal_label);
tmp = find(normal_label==inf);
normal_label(tmp)=[];
normal_label = intersect(normal_label,label_set);
checkerboard_normal = zeros(1,3);
for i=1:length(normal_label)
    ll = normal_label(i);
    ind_tmp = find(label_set==ll);
    checkerboard_normal = checkerboard_normal + normal_map(ind_tmp,:);
end


end

