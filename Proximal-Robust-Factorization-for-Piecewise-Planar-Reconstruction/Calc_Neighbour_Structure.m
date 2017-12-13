function [ neighbour_structure ] = Calc_Neighbour_Structure( label_map )
%SCENE_CONSTRAINT Summary of this function goes here
%   Detailed explanation goes here
    [H,W]=size(label_map);
    label=unique(label_map(:));
    label(find(label==inf))=[];
    se=strel('square',2);
    neighbour_structure=zeros(length(label));
    for i=1:length(label)
        ll_ind=find(label_map==label(i));
        mask=zeros(size(label_map));
        mask(ll_ind)=1;
        mask_dilate=imdilate(mask,se);
        tmp=mask_dilate-mask;
        indtmp=find(tmp);
        nn=label_map(indtmp);
        ind=find(nn==inf);
        if isempty(ind)==0
            nn(ind)=[];
        end
        neighbour=unique(nn);
        for k=1:length(neighbour)
            rrr=find(label==neighbour(k));
            neighbour_structure(rrr,i)=-1;
        end
    end
end