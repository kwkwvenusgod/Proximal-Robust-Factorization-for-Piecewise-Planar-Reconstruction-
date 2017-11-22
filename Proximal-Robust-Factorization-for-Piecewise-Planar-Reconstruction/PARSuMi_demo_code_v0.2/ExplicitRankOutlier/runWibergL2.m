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
Wcap=data_matrix;%(1:34,1:80);
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

%% Run Wiberg L2
%H is sample_mask, Y is Wcap
%m = 50;
%n = 30;
%r = 3;
cur_folder=pwd;
addpath(strcat(cur_folder,'\Wiberg L2\'));
sigma = 0.01;

H=sample_mask;
% Initial values
Vini = rand(n, r);
Uini = rand(m, r);

fprintf('Start Wiberg L2:\n');
fprintf('Checking the completeness of H...\n');
data_ok = datacheck_nomean(H, r);
fprintf('done.\n');
if data_ok == 1,
    fprintf('Starting wiberg_nomean...\n');
    %profile on;
    
    tic;
    
    [U2,V2] = wiberg_nomean(Wcap, H, r, Vini);
    
    timeElapsed_WibergL2=toc;
    V2=V2';
    %[U2,V2] = wiberg_nomean_naive(Wcap, H, r, Uini, Vini, 1e-6, 30);
    fprintf('complete.\n');
end

[hours,minutes,seconds] = Sec2HMS(timeElapsed_WibergL2);
Y = U2*V2;
fprintf('Time elapsed = %d hours : %d minutes : %d seconds \n',hours,minutes,seconds);
res=norm(sample_mask.*(Y-Wcap),'fro');
fprintf('\nFrobenious norm of error = %f \n\n',res);


%% display feature trajectory

figure(1)
%for Wiberg L2
hold off;
for j=1:n
    trk=Y(:,j);    
    idx = find(trk~=0);
    trk1=trk(idx);%get the sampled trk1
    sz=size(trk1,1);
    plot(trk1(1:2:sz),trk1(2:2:sz));
    title('Wiberg L2');
    hold on;
end
hold off;
