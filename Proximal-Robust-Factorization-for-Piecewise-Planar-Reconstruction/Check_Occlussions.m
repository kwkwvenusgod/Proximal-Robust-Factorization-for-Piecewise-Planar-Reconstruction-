function [ overoccluded ] = Check_Occlussions( flowf,flowb )
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
check = (check- min(check(:)))/max(check(:));

mu = mean(check(:));
variance= var(check(:));
[B,index] = sort(check(:),'descend');

% select top per as occlusions
per = 0.2;
num = ceil(length(check(:))* per);
index = index(num+1:end);

check(index)=0;

overoccluded=zeros(size(check));
overoccluded(find(check))=1;
tmp=overoccluded;
se=strel('square',8);
tmp=imerode(tmp,se);
overoccluded(find(tmp==1))=0;
end

