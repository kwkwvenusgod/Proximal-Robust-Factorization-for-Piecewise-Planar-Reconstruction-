function [ Wr,Er,N, Ur, Vr] = Plane_Factorization( M, mask, r, N0,W_ini, E_ini ,focal, label_map, label_set , label_map_full ,imageref ,patch_flowConf,uv,uvback )

%%
% INITIALIZE PARAMETERs

[m,n]=size(M);
p=sum(mask(:));
iter=1;
maxiter=2000;
beta1=1e-3/sqrt(max(m,n));
beta2=1e-3*1e-3/sqrt(max(m,n));
beta3=1e-3/sqrt(max(m,n));
beta4=1e-3/sqrt(max(m,n));
tol=1e-6;
gamma=1/sqrt(max(m,n));
tau=1;
thresh=0;

vecM=M(mask);
Ek=E_ini(mask);
[N_ini, S, V]=svds(W_ini,r);
Nk=N_ini;
C=S*V';
Wk=Nk*C;Wk=Wk(mask);

U=N_ini*sqrt(S);
V=sqrt(S)*V';
converged=0;
%%
lmax=3;
lambda=1e-6;
Wk_tmp=reshape(Wk,m,n);

depth=Calc_depth(V,focal,label_set,label_map);

s_depth=depth./abs(depth+eps);
s_depth(find(s_depth<0))=0;
vote=sum(sum(double(s_depth)));
vote=vote/(size(s_depth,1)*size(s_depth,2));
if vote<0.5
    U(:,1:3)=-U(:,1:3); V(1:3,:)=-V(1:3,:);
end
[cons_A,value_a, cons_B, value_b] = Construct_Equivalent_Constraint(length(unique(label_set)),focal);
patch_similarity = Construct_Patch_Similarity(uv,uvback,label_map, imageref);
[col,row]=find(patch_similarity);
ind=find(patch_similarity);
col=36*col-35*ones(length(col),1);
row=36*row-35*ones(length(row),1);
smoothness_Z0=sparse(row,col,patch_similarity(ind),36*length(label_set),36*length(label_set));
col=col+12*ones(length(col),1);
row=row+12*ones(length(row),1);
smoothness_normalxz=sparse(row,col,patch_similarity(ind),36*length(label_set),36*length(label_set));
col=col+12*ones(length(col),1);
row=row+12*ones(length(row),1);
smoothness_normalyz=sparse(row,col,patch_similarity(ind),36*length(label_set),36*length(label_set));
plane_similarity_constraint = [smoothness_Z0/focal^1.4;smoothness_normalxz;smoothness_normalyz];

%calc depth constraint
depth_similarity_constraint =  Construct_Depth_Constraint_In_Patch( label_map, uv, uvback ,patch_similarity,imageref, label_map_full )/focal^1.3;


count=1;
flowconfidencW=repmat(patch_flowConf,m,1);
flowconfidencW=repmat(flowconfidencW',1,6);
flowconfidencW=flowconfidencW';
flowconfidencW=reshape(flowconfidencW(:),m,n);


while ~converged
    
    Wk_tmp=reshape(Wk,m,n);
    count=count+1;
    
    [V_new ] = MF_CVX( Wk_tmp,cons_A,value_a, cons_B, value_b,plane_similarity_constraint,depth_similarity_constraint,U,patch_flowConf);
    
    %refine step
    beta2=beta2*1.02;
    WeightME=(M(mask)-Ek).*flowconfidencW(mask);
    [ Wk_new1 ,  Nk1, C1, normWk1 ] = refine_M( WeightME, Wk,r,Nk*C, mask, beta1, U, V_new, beta2, flowconfidencW(:));
    
    %update U
    tmpW=reshape(Wk_new1,m,n);
    U_new=tmpW*V_new'/(V_new*V_new');
    
    WeightEM=(vecM-Wk_new1).*flowconfidencW(mask);
    
    [Ek_new1, maskE1, normEk1]=Prox_ErrorCorrection(WeightEM,Ek,N0,beta4,thresh,flowconfidencW(:));
    
    e1=norm(Wk_new1(:)-Wk(:))/normWk1;
    
    tmp=U_new*V_new;
    tmp=tmp(mask);
    e2=norm(Wk_new1(:)-tmp(:))/norm(tmp(:));
    e3=norm(Wk(:)-tmp(:))/norm(tmp(:));
    V_new_measure=V_new(1:3,:);
    V_measure=V(1:3,:);
    e4=norm(V_new_measure(:)-V_measure(:))/norm(V_new_measure(:));
    e5=norm(U_new(:)-U(:))/norm(U_new(:));
    
    
    
    Nk=Nk1;
    C=C1;
    Wk=Wk_new1;
    Ek=Ek_new1;
    
    U=U_new;
    V=V_new;
    iter=iter+1;
    depth=Calc_depth(V,focal,label_set,label_map);
    normal=Calc_Normal(V,focal);
    
    s_depth=depth./abs(depth+eps);
    s_depth(find(s_depth<0))=0;
    vote=sum(sum(double(s_depth)));
    vote=vote/length(find(s_depth~=0));
    
    if vote<0.5
        U(:,1:3)=-U(:,1:3); V(1:3,:)=-V(1:3,:);
    end
    fprintf('The convergence rate of each element  e1 %f e2 %f e3 %f \n e4 %f e5 %f iter %i \n',e1,e2,e3,e4,e5,iter);
    if iter>maxiter
        converged=1;
        fprintf('Max iterations reached and converge may fail; please check final convergence. Quit.\n');
    end
    if e1<1e-5&&e5<0.5*1e-4&&e2<0.5*1e-4&&e4<0.5*1e-4
        fprintf('The convergence has been complete\n ');
        break
    end
end

Wr=reshape(Wk,size(M));
Er=reshape(Ek,size(M));
N=Nk;
Ur=U;
Vr=V;

end

