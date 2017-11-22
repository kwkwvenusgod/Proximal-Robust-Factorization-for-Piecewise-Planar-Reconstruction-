function [C]=compute_C_H(N,B,H)
% this function compute the corresponding pseudo inverse of N
% at every column, the output is a cell list of all pinv(N_i)
[m,n]=size(H);
C=zeros(size(N,2),n);
    for i=1:n
        y=B(:,i);
        Ni=diag(H(:,i))*N;
        C(:,i)=Ni\y;
    end
end