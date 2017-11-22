function [Wk_new, Nk,C,normWk]=PM_step_H(M, Wk,r,matWk, mask,beta,lmbd)
% Compute the analytical solution for the majorized solution
%   Detailed explanation goes here
% M is the Measurement - E in vector form
% Wk is in matrix form


[m,n]=size(mask);
H=zeros(m,n);%the \bar{H} in the paper
H(mask)=sqrt(1+beta);
H(~mask)=sqrt(lmbd*(1+beta));

p=max(H,[],2);q=max(H,[],1);

Bk=zeros(m,n);
Bk(mask)=(M+beta*Wk)/(1+beta);
Bk(~mask)=lmbd*beta/(lmbd*(1+beta)+eps)*matWk(~mask);

Gk=H.^2.*(matWk-Bk);

Uk=diag(p.^(1/2))*matWk*diag(q.^(1/2));
% in this case multiplying with p and q are just constant
Uk=Uk-diag(p.^(-1/2))*Gk*diag(q.^(-1/2));

[U, S, V]=svds(Uk,r);
PrUk=U*S*V';
Wk_new=diag(p.^(-1/2))*PrUk*diag(q.^(-1/2));
[Nk, S, V]=svds(Wk_new,r);
C=S*V';
Wk_new=Wk_new(mask);
normWk=norm(diag(S));


end

