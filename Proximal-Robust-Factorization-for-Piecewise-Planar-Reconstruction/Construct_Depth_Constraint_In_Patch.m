function [ depth_cons ] = Construct_Depth_Constraint_In_Patch( label_map, uv, uvback ,patch_similarity,imageref, label_map_full )
%DEPTH_CONSTRAINT_CONSTRUCT Summary of this function goes here
%   Detailed explanation goes here
% process scenConstraint to make less similar patch more close in depth of boundary pixels

scenConstraint_ = patch_similarity;
[nRows,nCols] = size(scenConstraint_);
scenConstraint_(1:(nRows+1):nRows*nCols) = 0;
scenConstraint_ = abs(scenConstraint_);
max_scenConstraint_ = max(scenConstraint_);
max_scenConstraint_(find(max_scenConstraint_==0)) = eps;
scenConstraint_ = scenConstraint_./repmat(max_scenConstraint_,length(max_scenConstraint_),1);
scenConstraint_(find(scenConstraint_==0)) = 1;

scenConstraint_ = 1./scenConstraint_;
scenConstraint_ = (scenConstraint_+scenConstraint_')/2;

label_set = unique(label_map(:));
label_set(find(label_set==inf)) = [];
[H, W] = size(label_map);

%find borders
ll=zeros(H,W,'int32');
ll(1:end,1:end)=label_map_full;
border=is_border_valsIMPORT(double(reshape(ll+1, [H,W])));
border_labels = border.*label_map;
border_labels(isnan(border_labels))=0;
border_labels(isinf(border_labels))=0;
% get forward-backward flow check
[overoccluded]=Check_Occlussions(uv,uvback);

se=strel('square',3);
neighbour_structure=zeros(length(label_set));
global_discard_neighour_num=0;
for i=1:length(label_set)
    ll_ind=find(label_map_full==label_set(i));
    mask=zeros(size(label_map_full));
    mask(ll_ind)=1;
    mask_dilate=imdilate(mask,se);
    tmp=mask_dilate-mask;
    indtmp=find(tmp);
    nn=label_map_full(indtmp);
    ind=find(nn==inf);
    if isempty(ind)==0
        nn(ind)=[];
    end
    neighbour=unique(nn);
    if length(neighbour)>4
        c=zeros(length(neighbour),1);
        for k=1:length(neighbour)
            c(k)=length(find(nn==neighbour(k)));
        end
        [freq, indtmp]=sort(c,'descend');
        for tt=1:length(indtmp(5:end))
            global_discard_neighour_num=global_discard_neighour_num+1;
            discard_neighbour_pair(global_discard_neighour_num,1)=label_set(i);
            discard_neighbour_pair(global_discard_neighour_num,2)=neighbour(indtmp(tt+4));
        end
        
    end
    for k=1:length(neighbour)
        rrr=find(label_set==neighbour(k));
        neighbour_structure(rrr,i)=-1;
    end
  
end
if exist('discard_neighbour_pair','var')
    for i=1:size(discard_neighbour_pair,1)
        r=find(label_set==discard_neighbour_pair(i,1));
        c=find(label_set==discard_neighbour_pair(i,2));
        neighbour_structure(r,c)=0;
        neighbour_structure(c,r)=0;
    end
end
borders=is_border_valsIMPORT(double(reshape(label_map_full, [H,W])));
se=strel('square',3);
borders=imdilate(borders,se);
se=strel('square',5);
se1=strel('square',6);
neighbour_structure_occlusion=ones(size(neighbour_structure));
for i=1:length(label_set)
    nn=neighbour_structure(i,:);
    
    tmp=find(nn);
    for j=1:length(tmp)
        mask=zeros(size(label_map_full));
        mask(find(label_map_full==label_set(i)))=1;
        mask(find(label_map_full==label_set(tmp(j))))=1;
        mask=imerode(mask,se);
        tmpind=find(mask==1);
        borderline=borders(tmpind);
        borderline_ind=find(borderline==1);
        tmpborderline=zeros(size(borders));
        tmpborderline(tmpind(borderline_ind))=1;
        tmpborderline=imdilate(tmpborderline,se1);
        occlusion_in_borderline=overoccluded(find(tmpborderline==1));
        occluded_ind=find(occlusion_in_borderline~=0);
        ratio=length(occluded_ind)/length(borderline_ind);
               
        if ratio>0.4
            neighbour_structure_occlusion(i,tmp(j))=0;
            neighbour_structure_occlusion(tmp(j),i)=0;
        end
    end
end

% Constraint from color texture map
im=double(imageref)/255;
texture=rgb2gray(imageref);
texture=entropyfilt(texture);
texture=mat2gray(texture);
im1=im(:,:,1);
im2=im(:,:,2);
im3=im(:,:,3);

%color weights
coeff1=zeros(length(label_set),1);
coeff2=zeros(length(label_set),1);
coeff3=zeros(length(label_set),1);
%texture weights
coeff4=zeros(length(label_set),1);
for i=1:length(label_set)
    ll=label_set(i);
    ind=find(label_map_full==ll);
    mu1=mean(im1(ind));
    mu2=mean(im2(ind));
    mu3=mean(im3(ind));
    mu4=mean(texture(ind));
    coeff1(i)=mu1;
    coeff2(i)=mu2;
    coeff3(i)=mu3;
    coeff4(i)=mu4;
end
aff1=repmat(coeff1,1,length(label_set))-repmat(coeff1',length(label_set),1);
aff2=repmat(coeff2,1,length(label_set))-repmat(coeff2',length(label_set),1);
aff3=repmat(coeff3,1,length(label_set))-repmat(coeff3',length(label_set),1);
aff4=repmat(coeff4,1,length(label_set))-repmat(coeff4',length(label_set),1);

simlarity1=aff1.*aff1+aff2.*aff2+aff3.*aff3;
simlarity2=aff4.*aff4;
sigma1=0.04;
simlarity1=-simlarity1/(sigma1^2);
sigma2=0.25;
simlarity2=-simlarity2/(sigma2^2);
simlarity=exp(simlarity1+simlarity2);
simlarity(find(simlarity<10e-2)) = 0;
neighour_structure = neighbour_structure_occlusion.*neighbour_structure.*scenConstraint_;

count = 1;
count_pair = 1;
pair_flag = false;
for i=1:length(label_set)
    host_label = label_set(i);
    sceCons = neighour_structure(i,:);
    sceCons(i) = 0;
    neighbor_ind = find(sceCons);
    
    [host_y, host_x] = find(border_labels == host_label);
    host_y = H/2-host_y;
    host_x = W/2-host_x;
    for j= 1: length(neighbor_ind)
         neighbour_label = label_set(neighbor_ind(j));
         neighbour_weight = abs(sceCons(neighbor_ind(j)));
        % check whether the pair processed or not
        if exist('processed_clique_pair','var')
            for k=1:size(processed_clique_pair,1)
                if all([host_label,neighbour_label]==processed_clique_pair(k,:))||all([neighbour_label,host_label]==processed_clique_pair(k,:))
                    pair_flag = true;
                    break;
                end
            end
            count_pair = count_pair+1;
            processed_clique_pair(count_pair,:)=[host_label,neighbour_label];
        else
            processed_clique_pair(count_pair,:)=[host_label,neighbour_label];
        end
        
        if pair_flag==true
            break;
        end
       
        [neighbour_y, neighbour_x] = find(border_labels == neighbour_label);
        neighbour_y = H/2-neighbour_y;
        neighbour_x = W/2-neighbour_x;
        dis_y = (repmat(host_y,1,size(neighbour_y,1)) - repmat(neighbour_y',size(host_y,1),1)).^2;
        dis_x = (repmat(host_x,1,size(neighbour_x,1)) - repmat(neighbour_x',size(host_x,1),1)).^2;
        dis = sqrt(dis_y+dis_x);
        [dis_ind_row, dis_ind_col] = find(dis<2); 
        if ~isempty(dis_ind_row)
            tmp = [neighbour_weight*host_y(dis_ind_row),neighbour_weight*host_x(dis_ind_row),repmat(host_label,length(dis_ind_row),1),...
                -neighbour_weight*neighbour_y(dis_ind_col),-neighbour_weight*neighbour_x(dis_ind_col),repmat(neighbour_label,length(dis_ind_row),1)];
            depth_cons_coord (count:count+length(dis_ind_row)-1,:)= tmp;
            count = count+length(dis_ind_row);
        end
    end
end

%construct sparse depth constraint matrix by sparse(i,j,s,m,n)
m = size(depth_cons_coord, 1);
n = length(label_set)*36;
ind_r = zeros(size(depth_cons_coord, 1)*6,1);
ind_c = zeros(size(depth_cons_coord, 1)*6,1);
value = zeros(size(depth_cons_coord, 1)*6,1);
for i = 1:6: 6*m
    row = floor(i/6)+1;
    tmp = depth_cons_coord(row,:);
    host_label = tmp(3); neighbor_label = tmp(6);
    ind_r(i) = row; ind_c(i) = find(label_set == host_label)*36 -35; value(i) = 1; 
    ind_r(i+1) = row; ind_c(i+1) = find(label_set == host_label)*36 -23; value(i+1) = tmp(1);
    ind_r(i+2) = row; ind_c(i+2) = find(label_set == host_label)*36 -11; value(i+2) = tmp(2);
    
    ind_r(i+3) = row; ind_c(i+3) = find(label_set == neighbor_label)*36 - 35; value(i+3) = -1;
    ind_r(i+4) = row; ind_c(i+4) = find(label_set == neighbor_label)*36 - 23; value(i+4) = tmp(4);
    ind_r(i+5) = row; ind_c(i+5) = find(label_set == neighbor_label)*36 - 11; value(i+5) = tmp(5);
end
depth_cons = sparse(ind_r,ind_c,value,m,n);

end