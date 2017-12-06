function [ an_error,error_rate,an_error_tr, an_error_tr_rate,angular_error_rotation, epe] = Evaluation( groundtruthfile ,normal,motion,superpixel,masklabel,normal_mask)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rotations_gt=load([groundtruthfile,'\Rotation_Para']);
translations_gt=load([groundtruthfile,'\Translation_Para']);
rotations_gt=struct2cell(rotations_gt);
translations_gt=struct2cell(translations_gt);
for i=1:size(translations_gt,1)
    tgt(i,:)=translations_gt{i}';
end

r1=rotations_gt{1};
dcm_anglegt=zeros(size(rotations_gt,1)-1,3);
for i=2:size(rotations_gt,1)
    ri=rotations_gt{i};
    rr=ri/(r1);
    dcm_anglegt(i-1,:)=dcm2angle(rr);
    tmp=rr*tgt(1,:)'-tgt(i,:)';
    tr_r_gt(i-1,:)=tmp';
    rr=rr-eye(3);
    dcm_anglegt(i-1,1)=-rr(2,3);
    dcm_anglegt(i-1,2)=rr(1,3);
    dcm_anglegt(i-1,3)=-rr(1,2);
end
mag_relative_tgt=sqrt(sum(tr_r_gt.*tr_r_gt,2));
directiongt=tr_r_gt./repmat(mag_relative_tgt,1,3);
directiongt(:,1:2)=-directiongt(:,1:2);
% directiongt(:,1)=-directiongt(:,1);
% directiongt(:,2)=-directiongt(:,2);

normal1_gt=r1*[0;0;1];
normal1_gt=[-normal1_gt(1);-normal1_gt(2);normal1_gt(3)];


%evaluation of normal

graymask=rgb2gray(normal_mask);
ind=find(graymask==255);
normal_label=superpixel(ind);
normal_label=unique(normal_label);
tmp=find(normal_label==inf);
normal_label(tmp)=[];
normal_label=intersect(normal_label,masklabel);
chess_normal=zeros(length(normal_label),3);

for i=1:length(normal_label)
    ll=normal_label(i);
    ind=find(masklabel==ll);
    chess_normal(i,:)=normal(ind,:);
end
chess_normal=mean(chess_normal);
chess_normal=chess_normal/sqrt(chess_normal*chess_normal');
ab_error=chess_normal*normal1_gt;
an_error=acos(ab_error);
error_rate=180*an_error/pi;
% evaluation of translation

translation=motion(:,1:3);
magnitute=translation.*translation;
magnitute=sqrt(sum(magnitute,2));
direction=translation./repmat(magnitute,1,3);

err_tr=sum(directiongt(1:end,:).*direction,2);
an_error_tr=acos(err_tr);
an_error_tr_rate=180*an_error_tr/pi;

%angular error of  rotation
rotation=motion(:,4:6);
r_angle=rotation./repmat(sqrt(sum(rotation.*rotation,2)),1,3);
dcm_anglegt(1:2,:)=-dcm_anglegt(1:2,:);
rgt_angle=dcm_anglegt./repmat(sqrt(sum(dcm_anglegt.*dcm_anglegt,2)),1,3);
angular_error_rotation=sum(r_angle.*rgt_angle(1:end,:),2);
angular_error_rotation=acos(angular_error_rotation);
angular_error_rotation=180*angular_error_rotation/pi;

%EPE of rotation

epe = rotation - dcm_anglegt;
epe = epe.*epe;
epe = sum(epe,2);
epe = sqrt(epe);
end

