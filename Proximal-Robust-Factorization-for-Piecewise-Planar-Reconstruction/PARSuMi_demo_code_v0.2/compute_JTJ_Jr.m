function [JTJ Jr]=compute_JTJ_Jr(N,M_list,mask)
[m,n]=size(mask);
r=size(N,2);
Jr=zeros(m*r,1);
JTJ=zeros(m*r);
% if n>5*m*r
%     nlist=randsample(n,5*m*r);
% else
%    nlist=(1:n)';
% end
for i=1:n;
    y=M_list{i};
    mm=mask(:,i);
    Ni=N(mm,:);
    I=eye(m);
    L=I(mm,:);
    D=kron(eye(r),L);
    mi=size(Ni,1);
    Ai=pinv(Ni);
    ci=Ai*y;
    C=communication(length(y),r);
    tmp=eye(mi)-Ni*Ai;
    Ji=(kron(ci',tmp)+kron((tmp*y)',Ai')*C)*D;
    
    Jr=Jr+Ji'*(y-Ni*ci);
    JTJ=JTJ+Ji'*Ji;
end


    