function [ patch_similarity,patch_clique ] = Construct_Patch_Similarity( flowf,flowb, correspondence , imageref)
%SCENE_CONSTRAINT Summary of this function goes here
%   Detailed explanation goes here
[H,W]=size(correspondence);
label=unique(correspondence(:));
label(find(label==inf))=[];
se=strel('square',3);
neighbour_structure=zeros(length(label));
global_discard_neighour_num=0;
for i=1:length(label)
    ll_ind=find(correspondence==label(i));
    mask=zeros(size(correspondence));
    mask(ll_ind)=1;
    mask_dilate=imdilate(mask,se);
    tmp=mask_dilate-mask;
    indtmp=find(tmp);
    nn=correspondence(indtmp);
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
            discard_neighbour_pair(global_discard_neighour_num,1)=label(i);
            discard_neighbour_pair(global_discard_neighour_num,2)=neighbour(indtmp(tt+4));
        end
        
    end
    for k=1:length(neighbour)
        rrr=find(label==neighbour(k));
        neighbour_structure(rrr,i)=-1;
    end
end
if exist('discard_neighbour_pair','var')
    for i=1:size(discard_neighbour_pair,1)
        r=find(label==discard_neighbour_pair(i,1));
        c=find(label==discard_neighbour_pair(i,2));
        neighbour_structure(r,c)=0;
        neighbour_structure(c,r)=0;
    end
end

neighbour_structure_enlarge=neighbour_structure;
for i=1:length(label)
    Nind=find(neighbour_structure_enlarge(:,i)==-1);
    if isempty(Nind)==0
        for k=1:length(Nind)
            NNind=find(neighbour_structure(:,Nind(k))==-1);
            newN=neighbour_structure_enlarge(:,i);
            newN(NNind)=-1;
            neighbour_structure_enlarge(:,i)=newN;
            neighbour_structure_enlarge(i,i)=0;
        end
    end
end

im=double(imageref)/255;
texture=rgb2gray(imageref);
texture=entropyfilt(texture);
texture=mat2gray(texture);
im1=im(:,:,1);
im2=im(:,:,2);
im3=im(:,:,3);

%color weights
coeff1=zeros(length(label),1);
coeff2=zeros(length(label),1);
coeff3=zeros(length(label),1);
%texture weights
coeff4=zeros(length(label),1);
for i=1:length(label)
    ll=label(i);
    ind=find(correspondence==ll);
    mu1=mean(im1(ind));
    mu2=mean(im2(ind));
    mu3=mean(im3(ind));
    mu4=mean(texture(ind));
    coeff1(i)=mu1;
    coeff2(i)=mu2;
    coeff3(i)=mu3;
    coeff4(i)=mu4;
end
aff1=repmat(coeff1,1,length(label))-repmat(coeff1',length(label),1);
aff2=repmat(coeff2,1,length(label))-repmat(coeff2',length(label),1);
aff3=repmat(coeff3,1,length(label))-repmat(coeff3',length(label),1);
aff4=repmat(coeff4,1,length(label))-repmat(coeff4',length(label),1);

simlarity1=aff1.*aff1+aff2.*aff2+aff3.*aff3;
simlarity2=aff4.*aff4;
sigma1=0.04;
simlarity1=-simlarity1/(sigma1^2);
sigma2=0.25;
simlarity2=-simlarity2/(sigma2^2);
simlarity=exp(simlarity1+simlarity2);
simlarity(find(simlarity<0.05*1e-1))=0;

overoccluded=Check_Occlussions(flowf,flowb);

borders=is_border_valsIMPORT(double(reshape(correspondence, [H,W])));
se=strel('square',3);
borders=imdilate(borders,se);

se=strel('square',5);
se1=strel('square',6);
neighbour_structure_occlusion=ones(size(neighbour_structure));
for i=1:length(label)
    nn=neighbour_structure(i,:);
    
    tmp=find(nn);
    for j=1:length(tmp)
        mask=zeros(size(correspondence));
        mask(find(correspondence==label(i)))=1;
        mask(find(correspondence==label(tmp(j))))=1;
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

for i=1:length(label)
    Nind=find(neighbour_structure(:,i)==-1);
    Noccluded=neighbour_structure_occlusion(:,i);
    if isempty(Nind)==0&&isempty(Noccluded)==0
        for k=1:length(Nind)
            NNind=find(neighbour_structure(:,Nind(k))==-1);   
            if Noccluded(Nind(k))==0
                breakN=neighbour_structure_enlarge(:,i);
                for mm=1:length(NNind)
                    if isempty(find(neighbour_structure(NNind(mm),i)==-1))==0
                        NNind(mm)=0;
                    end
                end
                NNind(find(NNind==0))=[];
                breakN(NNind)=0;
                neighbour_structure_enlarge(:,i)=breakN;
                neighbour_structure_enlarge(i,:)=breakN;
                neighbour_structure_enlarge(i,i)=0;
            end          
        end
    end
end

simlarity=simlarity.*neighbour_structure_occlusion.*neighbour_structure_enlarge;
patch_similarity=simlarity;
ind_diag=sub2ind([length(label),length(label)],[1:length(label)]',[1:length(label)]');
patch_similarity(ind_diag)=-sum(simlarity,2);
%treat the neighbour equal importance
patch_clique=zeros(size(patch_similarity));
for i=1:size(patch_similarity,1)
    indtmp=find(simlarity(:,i));
    if isempty(indtmp)==0
        tmpvec=zeros(size(patch_similarity,1),1);
        tmpvec(indtmp)=-1;
        patch_clique(:,i)=tmpvec;
        patch_clique(i,i)=length(indtmp);
    end
end
end
