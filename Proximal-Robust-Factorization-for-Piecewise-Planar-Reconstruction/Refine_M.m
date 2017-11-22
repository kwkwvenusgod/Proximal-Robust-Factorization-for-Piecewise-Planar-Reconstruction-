function [ Wk_new ,  Nk, C, normWk ] = Refine_M( M, Wk,r,matWk, mask, beta1, U, V, beta2,flowconfidencW )
%REFINE_M Summary of this function goes here
%   Detailed explanation goes here
W_uv=U*V;

W_uv=W_uv(mask);
Bk=(M+beta1*Wk+beta2*W_uv)./(flowconfidencW+beta1+beta2);


Gk=(1+beta1)*(Wk-Bk);
Uk=sqrt(1+beta1)*matWk;
Uk(mask)=Uk(mask)-(1/sqrt(1+beta1))*Gk;

[Nk, S, V_s]=svds(Uk,r);
C=(1/sqrt(1+beta1))*S*V_s';

Wk_new=Nk*C;
Wk_new=Wk_new(mask);
normWk=norm(C(:));
end

