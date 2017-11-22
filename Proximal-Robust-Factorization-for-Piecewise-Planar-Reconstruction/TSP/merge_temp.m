
D = 1;

x1 = randn(20)*10 + 0;
x2 = randn(20)*10 + 200;
x1 = x1(:);
x2 = x2(:);

kappa = 1;
nu = 100;
theta = 100;
Delta = 10;



% merged together
x = [x1;x2];
pkappa = kappa + numel(x);
pnu = nu + numel(x);
ptheta = 1/pkappa * (kappa*theta + sum(x));
pDelta = 1/pnu * (nu*Delta + sum(x'*x) + kappa*theta'*theta - pkappa*ptheta'*ptheta);

cov = pDelta * pnu / (pnu+D+1);

j = 1:D;
lprob = 0.5*pnu*log(pnu*det(pDelta)) - 0.5*pnu*D*log(2) - (D*(D-1)/4*log(pi) + sum(gammaln(pnu/2 + (1-j/2)))) + ...
    -0.5*(pnu+D+1)*log(det(cov)) - 0.5*D*(D+pnu+1);
%lprob = lprob -0.5*D*log(2*pi) - 0.5*log(det(cov));






pkappa = kappa + numel(x1);
pnu = nu + numel(x1);
ptheta = 1/pkappa * (kappa*theta + sum(x1));
pDelta = 1/pnu * (nu*Delta + sum(x1'*x1) + kappa*theta'*theta - pkappa*ptheta'*ptheta);

cov = pDelta * pnu / (pnu+D+1);

j = 1:D;
lprob2 = 0.5*pnu*log(pnu*det(pDelta)) - 0.5*pnu*D*log(2) - (D*(D-1)/4*log(pi) + sum(gammaln(pnu/2 + (1-j/2)))) + ...
    -0.5*(pnu+D+1)*log(det(cov)) - 0.5*D*(D+pnu+1);
%lprob2 = lprob2 -0.5*D*log(2*pi) - 0.5*log(det(cov));


pkappa = kappa + numel(x2);
pnu = nu + numel(x2);
ptheta = 1/pkappa * (kappa*theta + sum(x2));
pDelta = 1/pnu * (nu*Delta + sum(x2'*x2) + kappa*theta'*theta - pkappa*ptheta'*ptheta);

cov = pDelta * pnu / (pnu+D+1);

j = 1:D;
lprob2 = lprob2 + 0.5*pnu*log(pnu*det(pDelta)) - 0.5*pnu*D*log(2) - (D*(D-1)/4*log(pi) + sum(gammaln(pnu/2 + (1-j/2)))) + ...
    -0.5*(pnu+D+1)*log(det(cov)) - 0.5*D*(D+pnu+1);
%lprob2 = lprob2 -0.5*D*log(2*pi) - 0.5*log(det(cov));

lprob
lprob2