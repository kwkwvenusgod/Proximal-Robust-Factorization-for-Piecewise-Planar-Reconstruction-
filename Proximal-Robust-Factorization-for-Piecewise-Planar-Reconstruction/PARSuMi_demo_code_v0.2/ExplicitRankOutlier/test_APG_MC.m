%
clear all;
clc;
%%
%generate random data for testing
m=500;n=500;r=5;
X=orth(randn(m,r))*randn(r,n);

rate=0.8;
d=floor(rate*m*n);
idx=randperm(m*n);
idx=idx(1:d);

mask=zeros(m,n);mask(idx)=1;


%%
gamma_W=1e-2;
[Wr,residue_X,n_iter] = APG_MC(X,mask,zeros(size(X)),gamma_W);

%     gamma_W=0.2; %regularization weight for nuclear norm approximation
%     gamma_E = 1e5; %regularization weight for sparse error
%     fprintf('Getting a good initial estimate by using convex nuclear norm function in place of rank\n');
%     [Wr,Er,residue_nuc,iter] = NuclearNormAPG(X,mask,zeros(size(X)),zeros(size(X)),gamma_W,gamma_E);
    
    norm(Wr-X,'fro')
