%test script
clc;
addpath ExplicitRankOutlier
%%
m=100;n=100;
r=4;
rate=0.5;
W=rand(m,r)*rand(r,n);
mask=rand(m,n)<rate;
sigma=0.01;
Z=sigma*randn(m,n);
E=2*(rand(m,n)-0.5);
orate=0.1;
E(rand(m,n)>orate*rate)=0;
E(~mask)=0;
%E=zeros(m,n);
M=mask.*(W+E+Z);

p=sum(mask(:));
%%
E0=zeros(size(E));
W0=M-E0;
Nini=randn(m,r);


% Initialization by nuclear norm relaxation
    gamma_W = 0.5; %regularization weight for nuclear norm approximation
    gamma_E = 1/sqrt(n); %regularization weight for sparse error
    fprintf('Getting a good initial estimate by using convex nuclear norm function in place of rank\n');
    [W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(M,mask,W0,E0,gamma_W,gamma_E);
   
%        W_nuc=W0;
%        E_nuc=E0;
%% Initialization by nuclear norm relaxation

%W_0=rand(m,n);
                gamma_W = 5e-1; %regularization weight for nuclear norm approximation
                gamma_E = 1/sqrt(n); %regularization weight for sparse error
                %r=4;
                fprintf('Getting a good initial estimate by using convex nuclear norm function in place of rank\n');
                %[W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(Wcap,sample_mask,W_0,E_0,gamma_W,gamma_E);
            %    [U,S,V] = svd(W_nuc,'econ');
            %    S_sqrt = sqrt(S(1:r,1:r));
            %    U0 = U(:,1:4)* S_sqrt;
            %    V0 = V(:,1:4)* S_sqrt;

            top=1;bot=0;%define range we want to tune gamma_W
            [W_nuc,E_nuc,residue_nuc,iter,gamma_W]=optimal_gamma_W_initialization(M,mask,W0,E0,r,top,bot);

%%
tic;
[Nini S V]=svds(W_nuc,r);
E0=E_nuc;
%[Wr,Er,N]=convergentARSuMi(M,mask,r,ceil(0.15*p),Nini,E0,max(3*sigma,0.001));
%[Wr, Er, t, conv]=RMC_GRASTA(M,mask,r); [N S V]=svds(Wr,r);
[Wr,Er,N]=pmPARSuMi(M,mask,r,ceil(0.15*p),W_nuc,E0,max(3*sigma,0.001));
%lmbd=1e-3*sigma/m/n*sum(mask(:));
%[Wr,Er,N]=pmPARSuMi1(M,mask,r,ceil(0.15*p),W_nuc,E0,lmbd,max(3*sigma,0.001));
toc;
%%
Err=norm(W-Wr,'fro');
Obj=0.5*norm(M(mask)-Wr(mask)-Er(mask))^2;
fprintf('Obj=%f, Err=%f\n',Obj,Err)
Oracle=sqrt((m+n-4)*r/p/(1-orate));
fprintf('Denoising ratio=%f, oracle ratio=%f\n',norm(W(mask)-Wr(mask))/norm(Z(mask)),Oracle)
figure(1)
plot(E(mask),Er(mask),'xr',E(mask),E_nuc(mask),'*b')
hold off
figure(2);plot(E(mask),'.c-')
hold on;plot(Er(mask),'ro')
hold on;plot(E_nuc(mask),'b.')
figure(3);imagesc(abs(Wr-W));
colorbar;
%%
hfig=figure(1);
set(hfig,'position',[100,100,400,400])
plot(E(mask),Er(mask),'xr',E(mask),E_nuc(mask),'*b','markersize',10)
hold on; plot([-1,1],[-1,1],'-c')
hlg=legend('ARSuMi E vs. True E','Nuclear Norm E vs. True E');
set(hlg,'fontsize',12,'location','northwest');
set(gca,'position',[0.13,0.1,0.85,.85])
xlabel('Injected Corruptions (ground truth)','fontsize',12)
ylabel('Estimated Corruptions','fontsize',12)
xlim([-1,1]);ylim([-1,1])
grid on


