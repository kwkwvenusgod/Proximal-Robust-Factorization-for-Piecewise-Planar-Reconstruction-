function [Ek_new, maskE, normEk]=Prox_ErrorCorrection(a,b,num,beta,thresh,flowconfidencW)
%implement algorithm 2 for l0 constrained error correction
D=beta/(1+beta)*(a-b).^2-a.*a-beta*(b.*b);
[sorted, idx]=sort(D,'ascend');
maskE=false(size(a));
maskE(idx(1:num))=1;

Ek_new=zeros(size(a));
Ek_new(maskE)=(a(maskE)+beta*b(maskE))./(flowconfidencW(maskE)+beta);
Ek_new(abs(Ek_new)<thresh)=0;
normEk=norm(Ek_new(maskE));
end