function [ V_new ] = MF_CVX( M,cons_A, value_a, cons_B,value_b,plane_similarity_constraint,depth_similarity_constraint,Uk,patch_flowConf)

F=size(M,1);
P=size(M,2);
P=P/6;

b=reshape(M,6*F*P,1);
WeightedM=zeros(size(M));
%construct sparse matrix A
row1=1:6*F*P;
row1=repmat(row1,6,1);
row1=reshape(row1,36*P*F,1);
col1=1:36*P;
col1=reshape(col1,6,6*P);
col1=repmat(col1,F,1);
col1=reshape(col1,36*F*P,1);
AH=6*F*P;
AW=36*P;
Uk_weighted=[];
row_w=[];
col_w=[];
for i=1:P
    patch_flowConftmp=patch_flowConf(i);
    Uk_weightedtmp=Uk.*repmat(patch_flowConftmp,F,6);
    Uk_weightedtmp=repmat(Uk_weightedtmp,1,6);
    Uk_weightedtmp=Uk_weightedtmp(:);
    Uk_weighted=[Uk_weighted;Uk_weightedtmp];
   
    
    row1=6*F*i-6*F+1:6*F*i;
    row1=reshape(row1,F,6);
    row1=repmat(row1,6,1);
    row1=reshape(row1,36*F,1);
    row_w=[row_w;row1];
    col1=36*i-35:36*i;  
    col1=repmat(col1,F,1); 
    col1=reshape(col1,36*F,1);
    col_w=[col_w;col1];
    WeightedM(:,i*6-5:i*6)=repmat(patch_flowConftmp,F,6).*M(:,i*6-5:i*6);
end

coeff_A=sparse(row_w,col_w,Uk_weighted,6*F*P,36*P);
b=reshape(WeightedM,6*F*P,1);

cvx_begin
    variable x(36*P)
    minimize( norm( coeff_A * x -b, 2 ) +CoeffSceneCons*norm(plane_similarity_constraint*x,1)+CoeffSceneCons*norm(depth_similarity_constraint*x,2))
    subject to
        cons_A * x == value_a
        cons_B * x == value_b
cvx_end
V_new=reshape(x,6,6*P);
end