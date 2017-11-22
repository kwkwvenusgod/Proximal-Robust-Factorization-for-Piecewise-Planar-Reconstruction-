function [Ai_list]=compute_Ai(N,mask)
% this function compute the corresponding pseudo inverse of N
% at every column, the output is a cell list of all pinv(N_i)
[m,n]=size(mask);
Ai_list=cell(n);
parfor i=1:n
    Ai_list{i}=pinv(N(mask(:,i),:));
end
end