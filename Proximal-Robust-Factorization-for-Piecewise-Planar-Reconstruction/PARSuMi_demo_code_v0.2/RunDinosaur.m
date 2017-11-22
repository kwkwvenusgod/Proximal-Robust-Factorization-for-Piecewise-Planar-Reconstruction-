% This script is written to run ChenPei's LM_X algorithm on dinosaur
% sequence and simulation
clear all;
current_directory=pwd;
path=sprintf('%s\\ExplicitRankOutlier',pwd);
addpath(path);%for nuclear norm initialization
path=sprintf('%s\\abmDataset',pwd);
addpath(path);%add data path
rng('default') 
% rng(100)

%% Read data matrix dinosaur and block segmentation from file
%information about the images
N=36;%number of frames
width=720;%1360;
height=576;%1024;

S=load('rotdino.mat','M_inliers');
data_matrix=full(S.M_inliers);

Wcap=data_matrix;

[m, n]= size(Wcap);
r=4;%targeted rank
%image(Wcap.*sample_mask*50);

fprintf('Size of matrix: %d x %d\n',m,n);
n_frames = m/2;
N=n_frames;

sample_mask=logical(Wcap);
[row,col] = find(sample_mask == true);
n_sample = size(row,1);

%normalization of data to [0,1]
Wcap(1:2:m,:)=(Wcap(1:2:m,:))/width;
Wcap(2:2:m,:)=(Wcap(2:2:m,:))/height;

W=Wcap;

Wcap=Wcap.*sample_mask;


stats1=sum(sample_mask,1);
stats2=sum(sample_mask,2);

fprintf('Sampling rate = %.2f\n',n_sample/(m*n));
fprintf('Sample per degree of freedom = %.2f\n',n_sample/(m+n-r)/r);
fprintf('Number of frames = %d\n',n_frames);
fprintf('minimum samples per column and per row are %d and %d\n',min(stats1),min(stats2));

num=1;


%Wcap_N=zeros(m,n,num);
UiniALL=zeros(m,r,num);
ViniALL=zeros(n,r,num);

%report variables: 1 for Wiberg  2 for ExplicitRank 3 for LM_GN
Z1=zeros(num,1);%fitting residual
T1=zeros(num,1);%run time
Wr1=zeros(m,n,num);%Recovered data matrix
C1=zeros(num,1);% convergence residual
%ExplicitRank
Z2=zeros(num,1);
T2=zeros(num,1);
Wr2=zeros(m,n,num);
C2=zeros(num,1);
%LM_GN
Z3=zeros(num,1);
E3=zeros(m,n,num);
T3=zeros(num,1);
Wr3=zeros(m,n,num);
C3=zeros(num,1);
%LM_S
Z4=zeros(num,1);
T4=zeros(num,1);
Wr4=zeros(m,n,num);
C4=zeros(num,1);
%LM_M
Z5=zeros(num,1);
T5=zeros(num,1);
Wr5=zeros(m,n,num);
C5=zeros(num,1);

%recording initialization performance
E_nuc_base=zeros(m,n,num);
gamma_W_rec=zeros(num,1);
E_rec=zeros(m,n,num);
E_nuc_rec=zeros(m,n,num);
E_diff=zeros(num,1);
E_diff_base=zeros(num,1);

i=1;
%% run tests
for i=1:num
    
    
       % Initial values
        Vini = rand(n, r);
        ViniALL(:,:,i)=Vini;
        Uini = rand(m, r);
        UiniALL(:,:,i)=Uini;
        
        
        
%% generating error
d=floor(0.01*m*n);
idx=randperm(m*n);
E=zeros(m,n);
%E(idx(1:d))=0.01*randn(d,1);

%E(idx(1:d))=rand(1,d)-Wcap(idx(1:d));
E(idx(1:d))=4*(rand(d,1)-0.5);
%E=zeros(m,n);
%E=E+0.005*randn(m,n);
E=E.*sample_mask;

E_rec(:,:,i)=E;
%E(3,9)=1;
%E(4,9)=1;
%  Wcap(3,9)=1;
%  Wcap(4,9)=1;
Wcap=Wcap.*sample_mask;
Wcap=Wcap+E;
        
 %% Initialization of the nuclear norm initialization

% Change flag value to use different initialization
flag=1;
E_0 = zeros(m,n);

% Initialization
if flag==0
    W_0 = zeros(m,n);
elseif flag ==1;
    % Initialization using svd of image center
    %W_0 = zeros(m,n);
    W_0=Wcap;
    W_0(~sample_mask)=0.5;
    %[U S V] =svds(W_0,r);
    %W_0=U*S*V';
    
elseif flag==2
    %initialization using the same W as WibergL2
    W_0=Uini*Vini';
else
    % initialization of data matrix using interpolation and extrapolation
    W_0=Wcap;
    for j=1:n

        trk=W_0(:,j);    
        idx = find(trk>0);
        trk1=trk(idx);%get the sampled trk1
        sz=size(trk1,1);
        %obj=fit(trk1(1:2:sz),trk1(2:2:sz),'poly1'); %L2 Linear regression
        %P=[obj.p1 -1 obj.p2];
        obj=L1LinearRegression(trk1(1:2:sz),trk1(2:2:sz)); %L1 linear regression
        
        P=[obj(2) -1 obj(1)];
        xx=vec2mat(trk1',2)';
        xx1=[xx;ones(1,sz/2)];
        lambda=-P*xx1/(P(1)^2+1);
        xx2=xx+[P(1) ;-1]*lambda;
        step=(xx2(1,sz/2)-xx2(1,1))/(idx(sz)-idx(1))*2;    
        x0=xx2(1,1)-((idx(1)+1)/2-1)*step;
        trk_x=x0+(1:N)*step;
        trk_y=trk_x*P(1)+P(3);
        idx2=find(W_0(:,j)==0);
        sz=size(idx2,1);
        for i=idx2'
            if mod(i,2)==1
                W_0(i,j)=trk_x((i+1)/2);
            else
                W_0(i,j)=trk_y((i)/2);
            end         
        end

    end

end
       
   
%% Initialization by nuclear norm relaxation

    gamma_W=0.2; %regularization weight for nuclear norm approximation
    gamma_E = 1/sqrt(n); %regularization weight for sparse error
    fprintf('Getting a good initial estimate by using convex nuclear norm function in place of rank\n');
   [W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(Wcap,sample_mask,W_0,E_0,gamma_W,gamma_E);
   
    E_nuc_base(:,:,i)=E_nuc;
    E_diff_base(i)=norm(E_nuc(:)-E(:));
   
    top=1;bot=0;%define range we want to tune gamma_W
    
% uncomment the following line if want to use the rank guided tuning
% heuristic (performance is not guaranteed)

%[W_nuc,E_nuc,residue_nuc,iter,gamma_W]=optimal_gamma_W_initialization(Wcap,sample_mask,W_0,E_0,r,top,bot);



gamma_W_rec(i)=gamma_W;
E_nuc_rec(:,:,i)=E_nuc;
E_diff(i)=norm(E_nuc(:)-E(:));

%% Outlier preliminary removal heuristic
% doesn't matter much in most cases, may help in really harsh condition
E1norm=abs(E_nuc);

index=find(abs(E_nuc)>0.1);

mask_clean=sample_mask;
mask_clean(index)=0;
p=sum(sum(sample_mask));

%% run PARSuMi
    
% Transpose data matrices into tall matrix (suitable for LM_X)
    M=Wcap;
    mask=mask_clean;
    M_ini=W_nuc;
    E_ini=E_nuc.*mask;
    
    tic;
    %[Wr, Er, t, conv]=RMC_GRASTA(M,mask,r); [N S V]=svds(Wr,r);
    [Wr,Er,N]=pmPARSuMi(M,mask,r,ceil(0.1*p),W_nuc,E_ini,0.001);
    %lmbd=1e-6*0.001/m/n*sum(mask(:))
    %[Wr,Er,N]=pmPARSuMi1(M,mask,r,ceil(0.1*p),W_nuc,E_ini,lmbd,0.001);
    timeElapsed=toc;
    
    T3(i)= timeElapsed;
    Z3(i)=norm(M(mask)-Wr(mask)-Er(mask))^2;
    E3(:,:,i)=Er;
    %M=N*S
    %now only return the factor N, Based on N, calculate M
    Wr3(:,:,i)=Wr;


end
%% Display feature trajectory
for i=1:num
plot_traj( Wr3(:,:,i) );
title('recovered feature trajectory')
[hours,minutes,seconds] = Sec2HMS(T3(i)); 
fprintf('This is %dth case!\n',i);
fprintf('Time elapsed = %d hours : %d minutes : %d seconds \n',hours,minutes,seconds);
fprintf('Residual is %f. Frobenious norm error is %.6f. Sparse Error term is %.6f. \n', C3(i),Z3(i),norm(vec(E3(:,:,i)),1));
fprintf('Singular values are \n')
svds(Wr3(:,:,i),10)
% waitforbuttonpress 
end

figure;
subplot(2,1,1)
plot(Er(:))
title('recovered outliers')
subplot(2,1,2)
plot(E(:).*mask_clean(:))
title('ground truth outliers')
figure;
plot((Er(:)-E(:)).*mask_clean(:))
title('difference between recovered and ground truth outliers')
figure;plot(vec(mask_clean.*(Wr3(:,:,i)-Wcap+Er)))
title('recovered gaussian noise')

%%
mask=sample_mask;
Err=norm(W-Wr,'fro');
Obj=0.5*norm(M(mask)-Wr(mask)-Er(mask))^2;
fprintf('Obj=%f, Err=%f\n',Obj,Err)

figure(2)
plot(E(mask),Er(mask),'xb',E(mask),E_nuc(mask),'*r')
hold off
%
hfig=figure(3);
set(hfig,'position',[0 0 800 400])
hold off;
plot(E(mask),'.c-')
hold on;plot(Er(mask),'ro')
hold on;plot(E_nuc(mask),'b.')
hleg=legend('Injected E (Ground Truth)','ARSuMi Recovered E','Nuclear Norm Init. E');
set(hleg,'fontsize',13,'orientation','horizontal','location','southoutside')
set(gca,'position',[0.05,0.15,0.9,0.82])
xlim([0,5302])

%
hfig=figure(4);
set(hfig,'position',[0 0 800 400])
hold off;
plot((E_nuc(mask)-E(mask))*sqrt(width*height),'.-b')
hold on;
plot((Er(mask)-E(mask))*sqrt(width*height),'.-r')
hleg=legend('Initial E - GroundTruth E','ARSuMi E - GroundTruth E');
set(hleg,'fontsize',13,'orientation','vertical','location','southeast')
set(gca,'position',[0.05,0.1,0.9,0.85])
xlim([0,5302])
ylim([-100,100])

