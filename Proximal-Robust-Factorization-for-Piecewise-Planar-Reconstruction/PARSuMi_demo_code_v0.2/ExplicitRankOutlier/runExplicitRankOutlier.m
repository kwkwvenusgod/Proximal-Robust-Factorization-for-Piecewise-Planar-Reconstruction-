clear;clc;

%% Read dinosaur sequence from file
% This is a data matrix with missing data, but no outlier error.
load viff.xy;

i = find(viff== -1); % find where ends of tracks marked - indicator var
viff(i) = nan;       % changes these to nan

%plot feature trajectory
plot(viff(1,1:2:72),viff(1,2:2:72)) % plots track for point 1
plot(viff(1:end,1:2:72)',viff(1:end,2:2:72)') % plots all tracks
axis ij  % y in -ve direction
title('Feature point trajectories of the dinosaur sequence');

viff(i) = 0; % change it to 0 for our purpose.

%information about the images
N=36;%number of frames
width=720;%1360;
height=576;%1024;


%% Removing short feature tracks
data_matrix=viff';
idx=find(sum(data_matrix~=0)>12); %keep only tracks longer than 12
%i.e., keep only those features appear in at least 7 images
data_matrix=data_matrix(:,idx);
%The matrix has 336 feature tracks in total, this is what we deal with here

%display the data matrix
figure;image(data_matrix);
title('Data matrix');colorbar;

%% Read data matrix and block segmentation from file

% to test on a smaller submatrix for better sample rate
Wcap=data_matrix(1:42,1:130);
[m, n]= size(Wcap);
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

Wcap=Wcap.*sample_mask;

r=4;%targeted rank

fprintf('Sampling rate = %.2f\n',n_sample/(m*n));
fprintf('Sample per degree of freedom = %.2f\n',n_sample/(m+n+r)/r);
fprintf('Number of frames = %d\n',n_frames);

%% Add sparse outliers
fraction = 4e-2;
[Wcap_outlier,E_gt,E_mask] = AddSparseOutliers(Wcap,sample_mask,fraction);

%% Initialization of CM's explicit rank 4

% Change flag value to use different initialization
flag=1;
E_0 = zeros(m,n);

% Initialization
if flag==0
    W_0 = zeros(m,n);
elseif flag ==1;
    % Initialization using svd of image center
    W_0 = zeros(m,n);
    W_0(~sample_mask)=0.5;
    [U S V] =svd(W_0);
    for i=5:m
        S(i,i)=0;
    end
    W_0 = U(:,1:r)*S(1:r,1:r)*(V(:,1:r)');
    
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

%% Run CM's explicit Rank 4
fprintf('Entering main Rosper loop:\n');
fprintf('**********************************\n');
tic;
gamma_W = 1e-1; %regularization weight for nuclear norm approximation
gamma_E = 1e-2; %regularization weight for sparse error
% gamma_E = 1e1; %regularization weight for sparse error
r=4;
fprintf('Getting a good initial estimate by using convex nuclear norm function in place of rank\n');
[W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(Wcap_outlier,sample_mask,W_0,E_0,gamma_W,gamma_E);
% Wr = W_nuc;
% Er = E_nuc;
fprintf('Refine solution by using explicit rank constraint\n');
[Wr,Er,residue_xrank,n_iter] = ExplicitRank(Wcap_outlier,sample_mask,r,W_nuc,E_nuc,gamma_E);
timeElapsed = toc;

%%Compute the recovery accuracy and other statistics
[hours,minutes,seconds] = Sec2HMS(timeElapsed);
fprintf('Time elapsed = %d hours : %d minutes : %d seconds \n',hours,minutes,seconds);
fprintf('No. of iterations for nuclear norm initialization: %d\n',iter); 
fprintf('No. of iterations for explicit rank constraint: %d\n',n_iter); 

res2=norm(sample_mask.*(Wr-Wcap),'fro');
fprintf('Size of data matrix: %d x %d\n',m,n);
fprintf('\nFrobenious norm of error = %.6f \n\n',res2);

E2_mask=sample_mask.*Er;
%ppp=Wr+Er-Wcap;

% fprintf('Rosper: The L1 norm of sampled error is %.4f \n',sum(sum(abs(E2_mask))));
% fprintf('Rosper: The L2 norm of sampled error is %.4f \n',norm(E2_mask,'fro'));
% fprintf('Rosper: Mean L1 error for each sample is %.4f \n',sum(sum(abs(E2_mask)))/n_sample);
% fprintf('Rosper: RMS error for each sample is %.4f \n',norm(E2_mask,'fro')/sqrt(n_sample));


%% Visualization of theresults

figure(3)
subplot(2,2,1)
image(Wcap.*sample_mask*50);
title('Sampled Data Matrix ')
subplot(2,2,2)
image(Wr*50);
title('Recovered Rank-4 Data matrix')
subplot(2,2,3)
image(abs(Er)*50);
title('Recovered sparse noise term')
subplot(2,2,4)
image(abs(Wr-Wcap+Er).*sample_mask*100);
title('Recovered dense gaussian noise term')


% sv = svd(Wcap-Er);
% figure;plot(sv,'*');
% title('Singular values of Wcap - E');

%% display feature trajectory

figure(4)
% ExplicitRank
hold off;
for j=1:n
    trk=Wr(:,j);    
    idx = find(trk~=0);
    trk1=trk(idx);%get the sampled trk1
    sz=size(trk1,1);
    plot(trk1(1:2:sz),trk1(2:2:sz));
    title('Explicit rank');
    hold on;
end
hold off;