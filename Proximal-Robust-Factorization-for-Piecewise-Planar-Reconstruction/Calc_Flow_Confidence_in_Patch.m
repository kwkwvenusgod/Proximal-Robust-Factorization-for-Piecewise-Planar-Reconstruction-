function [ patch_confidence_map ] = Calc_Flow_Confidence_in_Patch( flowl2r,flowr2l,image_ref,image_set,label_ref,fovx,fovy )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
imref=double(image_ref)/255;

H=size(image_ref,1);
W=size(image_ref,2);
numofflow=size(flowl2r,2);
DuCheck=zeros(H,W,numofflow);
Gradientmap=zeros(H,W,numofflow);
for i=1:numofflow
    DuCheck(:,:,i)=checkocclusiton(flowl2r{i},flowr2l{i});
    im=double(image_set{i+1})/255;
    Gradientmap(:,:,i)=calcimgradient(imref,im);
end
masklabel=unique(label_ref(:));
masklabel(find(masklabel==inf))=[];
patch_confidence_map=zeros(numofflow,length(masklabel));
distweight=ones(1,length(masklabel));
distthreshx=fovx/2;
distthreshy=fovy/2;
alpha=5;
for i=1:length(masklabel)
    ind=find(label_ref==masklabel(i));
    [y,x]=find(label_ref==masklabel(i));
    y=mean(y);
    x=mean(x);
    x=abs(W/2-x);
    y=abs(H/2-y);
    if x>=distthreshx||y>=distthreshy
        if x>=y
            tmp=abs(x)/W;        
            distweight(i)=0.5*exp(-alpha*tmp);
        else
            tmp=abs(y)/H;        
            distweight(i)=0.5*exp(-alpha*tmp);
        end

%        distweight(i)=0.01;
    end
    for j=1:numofflow
        tmpgradient=Gradientmap(:,:,j);
        tmppatchgradient=tmpgradient(ind);
        tmppatchgradient=sort(tmppatchgradient,'descend');
        selectednum=ceil(0.20*length(tmppatchgradient));
        tmppatchgradient=tmppatchgradient(1:selectednum);
        tmpconfgradient=mean(tmppatchgradient);
        
        DuChecktmp=DuCheck(:,:,j);
        checktmp=DuChecktmp(ind);
        checktmp=sort(checktmp);
        maxnum=ceil(0.8*length(checktmp));
        checktmp=mean(checktmp(1:maxnum));
        
        tmppatchConf=checktmp*tmpconfgradient;
        patch_confidence_map(j,i)=tmppatchConf;
    end
end
patch_confidence_map=mean(patch_confidence_map);
patch_confidence_map=patch_confidence_map.*distweight;
end

function [check]=checkocclusiton(flowf,flowb)
Uf=flowf(:,:,1);
Vf=flowf(:,:,2);
H=size(Uf,1);
W=size(Uf,2);
[coorIm1to2X,coorIm1to2Y]=meshgrid(1:W,1:H);
coorIm1to2X=round(coorIm1to2X+Uf);
coorIm1to2Y=round(coorIm1to2Y+Vf);
coorIm1to2X(find(coorIm1to2X<=1))=1;
coorIm1to2X(find(coorIm1to2X>=W))=W;
coorIm1to2Y(find(coorIm1to2Y<=1))=1;
coorIm1to2Y(find(coorIm1to2Y>=H))=H;
Ub=flowb(:,:,1);
Vb=flowb(:,:,2);
ind=sub2ind(size(coorIm1to2X),coorIm1to2Y(:),coorIm1to2X(:));
Ub=Ub(ind);
Vb=Vb(ind);
Ub=reshape(Ub,H,W);
Vb=reshape(Vb,H,W);

checkU=Uf+Ub;
checkV=Vf+Vb;

check=sqrt(checkU.^2+checkV.^2);
%  sigma=std(check(:));
sigma=1.5;
check=(check/sigma);
check=exp(-(check).^2);
end
function  [gradientmap]=calcimgradient(imref,im)
% imrefIx=zeros(size(imref,1),size(imref,2));
% imrefIy=zeros(size(imref,1),size(imref,2));
gradientmapIx=zeros(size(imref,1),size(imref,2));
% imIx=zeros(size(im,1),size(im,2));
% imIy=zeros(size(im,1),size(im,2));
gradientmapIy=zeros(size(im,1),size(im,2));
h = fspecial('gaussian', 5, 5);
imref=imfilter(imref,h);
im=imfilter(im,h);
for i=1:size(imref,1)-2
    tmp1=2*imref(i,:,1)-imref(i+1,:,1)-imref(i+2,:,1);
    tmp2=2*imref(i,:,2)-imref(i+1,:,2)-imref(i+2,:,2);
    tmp3=2*imref(i,:,3)-imref(i+1,:,3)-imref(i+2,:,3);
    tmp4=2*im(i,:,1)-im(i+1,:,1)-im(i+2,:,1);
    tmp5=2*im(i,:,2)-im(i+1,:,2)-im(i+2,:,2);
    tmp6=2*im(i,:,3)-im(i+1,:,3)-im(i+2,:,3);
    gradientmapIy(i,:)=(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6)/4;
end

for i=1:size(imref,2)-2
    tmp1=2*imref(:,i,1)-imref(:,i+1,1)-imref(:,i+2,1);
    tmp2=2*imref(:,i,2)-imref(:,i,2)-imref(:,i+2,2);
    tmp3=2*imref(:,i,3)-imref(:,i,3)-imref(:,i+2,3);
    tmp4=2*im(:,i,1)-im(:,i+1,1)-im(:,i+2,1);
    tmp5=2*im(:,i,2)-im(:,i+1,2)-im(:,i+2,2);
    tmp6=2*im(:,i,3)-im(:,i+1,3)-im(:,i+2,3);
    gradientmapIx(:,i)=(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6)/4;
end
epson=0.001;
gradientmap=sqrt(gradientmapIy.^2+gradientmapIx.^2+epson^2);
rho=std(gradientmap(:));
alpha=4;
gradientmap=exp(alpha*rho*gradientmap);
% sigma=1;
% gradientmap=gradientmap/sigma;
% gradientmap=exp(-)
end
