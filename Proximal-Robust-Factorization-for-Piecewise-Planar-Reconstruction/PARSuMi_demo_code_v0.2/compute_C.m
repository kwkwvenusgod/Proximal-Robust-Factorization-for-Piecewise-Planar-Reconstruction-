function [C]=compute_C(N,M_list,mask)
% this function compute the corresponding pseudo inverse of N
% at every column, the output is a cell list of all pinv(N_i)
[m,n]=size(mask);
C=zeros(size(N,2),n);
    for i=1:n
        y=M_list{i};
        Ni=N(mask(:,i),:);
        C(:,i)=Ni\y;
    end
end