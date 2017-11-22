function [Wk_new C ]=huber_fitting_iterative(M,mask,N,sigma, lmax)
[m,n]=size(mask);
W=zeros(size(mask));
W(mask)=M;
C=zeros(size(N,2),n);
lambda=sigma;
%epsilon=0.1*sigma;
for i=1:n
    mi=W(:,i);
    mm=mask(:,i);
    y=mi(mm);
%    lambda=1/sqrt(ni);
    Ni=N(mm,:);
    weight=ones(length(y),1);
    for j=1:lmax
        ci=huber_fit(bsxfun(@rdivide,Ni,weight)/lambda,y./weight/lambda,2.0,1.0);
        res=abs(y-Ni*ci);
        weight(res>lambda)=max(1,(res(res>lambda)/lambda));
%         
%         plot(res,'*');
%         pause(0.1);
    end
    C(:,i)=ci;
    W(:,i)=N*C(:,i);
end
Wk_new=W(mask);
end