function [Wk_new C ]=huber_fitting(M,mask,N,sigma)
[m,n]=size(mask);
W=zeros(size(mask));
W(mask)=M;
C=zeros(size(N,2),n);
lambda=sigma;
for i=1:n
    mi=W(:,i);
    mm=mask(:,i);
    y=mi(mm);
%    lambda=1/sqrt(ni);
    Ni=N(mm,:);
    C(:,i)=huber_fit(Ni/lambda,y/lambda,2.0,1.0);
    W(:,i)=N*C(:,i);
end
Wk_new=W(mask);
end