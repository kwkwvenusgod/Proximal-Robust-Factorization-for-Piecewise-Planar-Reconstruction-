function [ Depth_map ] = Depth_Completion( label_map, label_map_consistent, structure_matrix, focal_length  )
%DEPTH_COMPLETION Summary of this function goes here
%   Detailed explanation goes here
label_set = unique(label_map(:));
label_set(label_set == inf) = [];
label_set_consistent = unique(label_map_consistent(:));
label_set_consistent(label_set_consistent == inf) = [];
[ depth, nagetive_label,label_map_consistent_new,label_consistent_new]  = Calc_Depth(structure_matrix, focal_length,label_set_consistent, label_map_consistent);
neighbour_structure  = Calc_Neighbour_Structure(label_map);
[values, ind1, ind2] = intersect(label_set, label_set_consistent);
missing_label = label_set;
missing_label(ind1) = [];

missing_label = [missing_label; nagetive_label];

[values, ind1, ind2] = intersect( label_set_consistent,label_consistent_new);
Z0 = Calc_Plane_Offset(structure_matrix, focal_length);
Normal = Calc_Normal(structure_matrix, focal_length);
label_set_completion = label_consistent_new;
Z0_completion = Z0(ind1,:);
Normal_completion = Normal(ind1,:);

while isempty(missing_label) == 0
    ml=missing_label(1);
    indhost=find(label_set==ml);
    %find neighbour
    nn=neighbour_structure(:,indhost);
    nn(indhost)=0;
    nn_index = find(nn);
    if isempty(nn)==1
       missing_label(1)=[]; 
    else
       count=0;
       Normal_Missiong_Labels = zeros(1,3);
       Z0_Missing_Labels = 0;
       find_missing = 0;
       for k=1:length(nn_index)
           nl=label_set(nn_index(k));
           indtmp=find(label_set_completion==nl);
           if isempty(indtmp)==0
               count=count+1;
               Normal_Missiong_Labels = Normal_Missiong_Labels + Normal_completion(indtmp,:);
               Z0_Missing_Labels = Z0_Missing_Labels + Z0_completion(indtmp,:);
               
               find_missing = 1;
           end
       end
       if  find_missing == 1
           normal_missing = Normal_Missiong_Labels/ norm(Normal_Missiong_Labels);
           Z0_missing = Z0_Missing_Labels/count;
           
           originallength = length(Z0_completion);
           Z0_completion(originallength+1,:) = Z0_missing;
           Normal_completion(originallength+1,:) = normal_missing;
           label_set_completion(originallength+1,:) = ml;
           missing_label(1)=[];
       else
           missing_label(1)=[];
           missing_label(length(missing_label)+1)=ml;
       end
    end
end

% calculate depth
Depth_map = inf*ones(size(label_map));
[H, W]= size(label_map);
for i = 1:length(label_set_completion)
    label_tmp = label_set_completion(i);
    normal_tmp = Normal_completion(i,:);
    Zx = -normal_tmp(1)/(normal_tmp(3)+eps);
    Zy = -normal_tmp(2)/(normal_tmp(3)+eps); 
    Z0_tmp = Z0_completion(i);
    ind_tmp = find(label_map == label_tmp);
    [y,x] = find(label_map == label_tmp);
    y=H/2*ones(length(y),1)-y;
    x=W/2*ones(length(x),1)-x;
    
    Depth_map(ind_tmp) = Z0_tmp./(1-x*Zx/focal_length-y*Zy/focal_length);
end


end

