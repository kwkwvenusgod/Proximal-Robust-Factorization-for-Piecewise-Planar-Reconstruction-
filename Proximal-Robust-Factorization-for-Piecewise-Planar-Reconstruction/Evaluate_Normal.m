function [ angular_error ] = Evaluate_Normal( normal_gt, normal_est)
%EVALUATE_NORMAL Summary of this function goes here
%   Detailed explanation goes here
normal_est_normalize = normal_est/norm(normal_est);
cos_similarity = normal_est_normalize * normal_gt;

angular_error = 180*acos(cos_similarity)/pi;
end

