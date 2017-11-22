function [JTJ Jr]=compute_JTJ_Jr_H(N,B,H)
[m,n]=size(H);
r=size(N,2);
Jr=zeros(m*r,1);
JTJ=zeros(m*r);

parfor i=1:n
    y=B(:,i);
    mm=H(:,i);
    Ni=diag(mm)*N;
    
    Ai=pinv(Ni);
    Ci=Ai*y;
    
    C=communication(length(y),r);
    tmp=(eye(m)-Ni*Ai);
    Ji=(kron(Ci',tmp*diag(mm))+kron((diag(mm)*(tmp*y))',Ai')*C);
    
    Jr=Jr+Ji'*(y-Ni*Ci);
    JTJ=JTJ+Ji'*Ji;
end


    