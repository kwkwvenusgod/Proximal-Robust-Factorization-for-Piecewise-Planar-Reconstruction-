function [Ai_list,C]=compute_Ai_C(N,M_list,mask)
% this function compute the corresponding pseudo inverse of N
% at every column, the output is a cell list of all pinv(N_i)
[m,n]=size(mask);
Ai_list=cell(n,1);
C=zeros(size(N,2),n);
    parfor i=1:n
        y=M_list{i};
        Ni=N(mask(:,i),:);
        pinvNi=pinv(Ni);
        C(:,i)=pinvNi*y;
        Ai_list{i}=pinvNi;
    end
end