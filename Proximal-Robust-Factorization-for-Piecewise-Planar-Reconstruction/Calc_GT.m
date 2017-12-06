function [ motion, normal ] = Calc_GT( translations_gt,rotations_gt)
%CALC_GT Summary of this function goes here
%   Detailed explanation goes here
tgt = zeros(size(translations_gt,1), 3);
for i=1:size(translations_gt,1)
    tgt(i,:)=translations_gt{i}';
end

r1=rotations_gt{1};
dcm_anglegt = zeros(size(rotations_gt,1)-1,3);
tr_r_gt = zeros(size(translations_gt,1)-1,3);
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
tr_r_gt(:,1:2) = -tr_r_gt(:,1:2);
dcm_anglegt(:,1:2) = -dcm_anglegt(:,1:2);
motion = [tr_r_gt, dcm_anglegt];
normal=r1*[0;0;1];
normal=[-normal(1);-normal(2);normal(3)];

end

