function [ depth, nagetivelabel,new_label_map,label_set_new] = Calc_Depth(structure_matrix,focal_len,label_set, label_map  )
%CAL_DEPTH Summary of this function goes here
%   Detailed explanation goes here


plane_Z0 = Calc_Plane_Offset(structure_matrix,focal_len );
normal = Calc_Normal(structure_matrix,focal_len );

depth=zeros(size(label_map,1),size(label_map,2));
H=size(label_map,1);
W=size(label_map,2);
for i =1:length(label_set)
        ll=label_set(i);
        normal_tmp=normal(i,:);
        Zx=-normal_tmp(1)/normal_tmp(3);
        Zy=-normal_tmp(2)/normal_tmp(3);
        ind= label_map==ll;
        [y,x]=find(label_map==ll);
        y=H/2*ones(length(y),1)-y;
        x=W/2*ones(length(x),1)-x;
        depth(ind)=plane_Z0(i)./(1-x*Zx/focal_len-y*Zy/focal_len);
end
ind= depth<0;
nagetivelabel=unique(label_map(ind));
new_label_map=label_map;
label_set_new=label_set;
if isempty(nagetivelabel)==0
    for i=1:length(nagetivelabel)
        new_label_map(new_label_map==nagetivelabel(i))=inf;
    end
    label_set_new=unique(new_label_map(:));
    label_set_new(label_set_new==inf)=[];
end
end

