% script to test RMC algorithms

%generate data
m=100;n=100; r=4; rate=0.3;orate=0.1;
sigma=0; mag=2;

U=randn(m,r);
V=randn(n,r);
W=U*V';

p=ceil(m*n*rate);
idx=randperm(m*n);
idx=idx(1:p);
mask=false(m,n);
mask(idx)=1;

po=ceil(p*0.1);
idxo=idx(1:po);
E=zeros(m,n);
E(idxo)=mag*(rand(po,1)-0.5);
masko=false(m,n);
masko(idxo)=1;


Z=sigma*randn(m,n);

M=W+E+Z;
M(~mask)=0;

%% run robust matrix completion
fprintf('Start RMC...\n')
[Wr, Er, t, conv]=RMC_GRASTA(M,mask,r);

%%

RMSE=norm(W-Wr,'fro')/sqrt(m*n);
figure(1);plot(E(:),Er(:),'x');
Obj=norm(W(mask)-Wr(mask),'fro');
fprintf('Done.\n')
fprintf('RMSE=%f, ObjVal=%f, conv=%d.\n',RMSE,Obj,conv);
