function [Wk_new, Nk,C,normWk]=PM_step(M, Wk,r,matWk, mask,beta)
% Compute the analytical solution for the majorized solution
%   Detailed explanation goes here
% M is the Measurement - E in vector form
% Wk is in matrix form


[m,n]=size(mask);

p=sqrt(1+beta)*ones(m,1);q=sqrt(1+beta)*ones(n,1);

Bk=(M+beta*Wk)/(1+beta); % vector form
Gk=(1+beta)*(Wk-Bk);%vector form

Uk=sqrt(1+beta)*matWk;
% in this case multiplying with p and q are just constant
Uk(mask)=Uk(mask)-(1/sqrt(1+beta))*Gk;

[Nk, S, V]=svds(Uk,r);
C=(1/sqrt(1+beta))*S*V';

Wk_new=Nk*C;
Wk_new=Wk_new(mask);
normWk=norm(C(:));


end

