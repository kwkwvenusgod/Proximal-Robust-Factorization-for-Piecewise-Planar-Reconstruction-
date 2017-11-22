%     [100,1000,10000],...
%     [100,1000,10000],...
%     [400],...
%     [-15,-5],...
%     [-3,-1],...
%     [0.1,0.33,1]);

params.cov_var_p = 1000;
params.cov_var_a = 1000;
params.area_var = 400;
params.Kpercent = 0.8;
params.K = 600;
params.alpha = -10;
params.beta = -1;
params.deltap_scale = 0.5;


% params.cov_var_p = 1000;
% params.cov_var_a = 1000;
% params.area_var = 400;
% params.Kpercent = 0.8;
% params.K = 600;
% params.alpha = -10;
% params.beta = -9;
% params.deltap_scale = 0.1;


params.cov_var_p = 100;
params.cov_var_a = 10000;
params.area_var = 400;
params.Kpercent = 0.8;
params.K = 200;
params.alpha = -15;
params.beta = -2;
params.deltap_scale = 1;

params.cov_var_p = 1000;
params.cov_var_a = 100;
params.area_var = 400;
params.Kpercent = 0.8;
params.K = 100;
params.alpha = -15;
params.beta = -10;
params.deltap_scale = 1e-3;
params.deltaa_scale = 100;



%close force all;
reset(RandStream.getGlobalStream());
setenv('GUROBI_HOME', '/csail/fisher3/projects/Vseg_sp/C/gurobi500/linux64');
username = getenv('USER');
[~,hostname] = system('hostname');
hostname(numel(hostname)) = [];
setenv('GRB_LICENSE_FILE', ['/csail/fisher3/licenses/gurobi/' username '/' hostname '/gurobi.lic']);

%mex test_mex.c
%make

addpath('/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_backup/MM_SP/flow/');
addpath('/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_backup/MM_SP/flow/optical_flow_celiu');
addpath('/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_backup/MM_SP/flow/optical_flow_celiu/mex');


%monkeydog
name = 'parachute';
name = 'girl';
% name = 'garden';
% name = 'monkeydog';
root = ['/csail/fisher5/datasets/SegTrack_201111_clean/' name '/'];
%root = ['/csail/fisher3/datasets/SegTrack_201111/' name '/'];
%root = '/csail/fisher3/datasets/SegTrack_201111/parachute/';
%root = '/csail/fisher3/datasets/SegTrack_201111/penguin/';

% name = 'dog';
% name = 'parachute';
% name = 'garden';
% name = 'sample1';
% root = ['/csail/fisher5/datasets/LayerSegmentation_clean/Data/Seg/' name '/'];

name = 'Urban2';
root = ['/csail/fisher5/datasets/SV/img/Middlebury/' name '/'];

files = dir([root '*.png']);
files = dir([root '*.bmp']);
files = dir([root '*.ppm']);
% files = dir([root '*.jpg']);


% root = ['/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_/table/frames/'];
% files = dir([root '*.jpg']);


% root = ['/csail/fisher3/datasets/garden/'];
% files = dir([root '*.ppm']);

frame_it = 0;
frame_start = 1;
frames = frame_start:1:numel(files);
frames = frame_start:1:10;
%frames = frame_start;
%frames = frame_start:1:numel(files);
%frames = repmat(frame_start, [1 20]);
% frames = 25:35;
% frames = [8 10 10 10 10 10];



% frames = 1:3;

% % middlebury
% root = '/csail/fisher3/datasets/Middlebury/other-data/';
% folders = {'RubberWhale', 'Hydrangea', 'Dimetrodon', 'Grove2', 'Grove3', 'Urban2', 'Urban3'};
% root = [root folders{1} '/'];
% files = dir([root '*.png']);
% frame_it = 0;
% frame_start = 1;
% frames = frame_start:1:numel(files);

% root = '/csail/fisher3/datasets/BSDS300/images/test/';
% name = '126007';
% name = '175043';
% name = '219090';
% % name = '87046';
% files = dir([root name '.jpg']);
% files = dir([root name '.jpg']);
% files = dir([root name '.jpg']);
% files = dir([root name '.jpg']);
% frame_it = 0;
% frame_start = 1;
% frames = 1;


%K = 4;

oim = imread([root files(frame_start).name]);
oim = imresize(oim, 1);
%oim = oim(1:40,1:40,:);
%oim = oim(1:100,1:100,:);
%oim(:) = randn(numel(oim),1)*10;
%oim = oim(1:10,1:40,:);

sp_labels = zeros(size(oim,1), size(oim,2), 'uint64');
sp_flow = cell(numel(frames), 1);
sp_flow{end} = zeros(size(oim));


% root = 'test/';
% files = dir([root '*.png']);
%
% oim = imread([root files(1).name]);
% sp_labels = zeros(size(oim,1), size(oim,2), 'uint64');
% frame_it = 0;
% frame_start = 1;
% frames = frame_start:1:2;

dispOn = true;

for f=frames
    frame_it = frame_it + 1;
    disp(['Frame ' num2str(f) ' ----------------------------']);
    oim1 = imread([root files(f).name]);
    oim1 = imresize(oim1, 1);
    
%     oim1 = oim1 + uint8(randn(size(oim1))*10);
    %oim1 = oim1(1:100,1:100,:);
    %oim1 = oim1(1:40,1:40,:);
    %oim1(:) = randn(numel(oim1),1)*10;
    %oim1 = oim1(1:10,1:40,:);


    if (frame_it==1)
        %IMG = IMG_init(rgb2lab(oim1),K);
        IMG = IMG_init(oim1, params);
%        IMG.log_alpha = -(IMG.xdim*IMG.ydim/K)^2/(12.5123*500);
%        IMG.log_alpha = -30;
    else

%         lab0 = rgb2lab(oim);
%         lab1 = rgb2lab(oim1);
%         lab0(:,:,1) = lab0(:,:,1) / IMG.hyper.a_Sigma(1);
%         lab0(:,:,2) = lab0(:,:,2) / IMG.hyper.a_Sigma(2);
%         lab0(:,:,3) = lab0(:,:,3) / IMG.hyper.a_Sigma(3);
%         lab1(:,:,1) = lab1(:,:,1) / IMG.hyper.a_Sigma(1);
%         lab1(:,:,2) = lab1(:,:,2) / IMG.hyper.a_Sigma(2);
%         lab1(:,:,3) = lab1(:,:,3) / IMG.hyper.a_Sigma(3);
%
%         [u, v, aep] = SP_likelihood_flow_IMG(IMG, lab0, lab1, true, [], [], 5, 5);
%
%         % optical flow returns actual x and y flow... flip it
%         vx = zeros(size(oim,1), size(oim,2));
%         vy = zeros(size(oim,1), size(oim,2));
        [vx vy] = coarse_flow(oim, oim1,0);
%        [vy vx] = coarse_flow(oim, oim1,10);
%         vx(:) = 0;
%         vy(:) = 0;


%         load(['/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp/table/flows/sample ' num2str(f, '%02d') '.mat']);

%         [u, v, ~] = SP_likelihood_flow_IMG(IMG, oim, oim1, true, [], [], 50, 50);


%         vx = zeros(size(oim1,1),size(oim1,2));
%         vy = zeros(size(oim1,1),size(oim1,2));

        IMG = IMG_prop(oim1,vy,vx,IMG);

        %newIMG = IMG_prop(oim1,vx,vy,IMG);
        %IMG = newIMG;
        %[u, v] = SP_likelihood_flow(IMG.label, oim, oim1, [], [], true, [], [], 5, 25, [], []);
%         mask = imread([root 'gt.png']);
%         u = double(mask)/255*10;
%         v = double(mask)/255*2;
%         IMG = IMG_prop(oim1,u,v,IMG);
    end

%     if (f==2)
%         IMG.alpha = 1.25;
%     end

    oim = oim1;
    if (frame_it>1)
        IMG.prev_oim = IMG.oim;
    end
    IMG.oim = oim;

    %sfigure(3);
    %imagesc(IMG.label);


    %[IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID, E(end+1)] = local_move(IMG,1000);

    E = [];
    it = 0;
    IMG.alive_dead_changed = true;
    IMG.SxySyy = [];
    IMG.Sxy = [];
    IMG.Syy = [];
    converged = false;
    disp('Split/Merge');
    while (~converged && it<5 && true && frame_it==1)
        it = it + 1;

        oldK = IMG.K;
        IMG.SP_changed(:) = true;
        [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = split_move(IMG,1);
        E(end+1) = newE;
        converged = IMG.K - oldK < 2;

        if (dispOn)
            sfigure(1);
            subplot(1,1,1);
            im = zeros(size(oim));
            %im(:) = c(IMG.label+1,:);
            imagesc(IMG.label);
            %image(im,'parent',gca);
            title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);

            sfigure(2);
            subplot(1,1,1);
            im = zeros(size(oim,1)+2*IMG.w, size(oim,2)+2*IMG.w, 3);
            im(IMG.w+1:end-IMG.w, IMG.w+1:end-IMG.w, :) = double(oim)/255;
            borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [IMG.xdim IMG.ydim])));
            im = setPixelColors(im, find(borders), [1 0 0]);
            image(im,'parent',gca);
            drawnow;
        end
    end



    it = 0;
    converged = false;
    if (frame_it>1)
        IMG.SP_changed(:) = true;
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = merge_move(IMG,1);
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = split_move(IMG,10);
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = switch_move(IMG,1);
        [IMG.K, IMG.label, IMG.SP, ~, IMG.max_UID, ~, ~, ~] = localonly_move(IMG,1000);
    end
    IMG.SP_changed(:) = true;
    IMG.alive_dead_changed = true;
    disp('Local');drawnow;
    while (~converged && it<20)
        it = it + 1;
        disp(['Iteration: ' num2str(it)]);drawnow;
        times = zeros(1,5);
        tic;[IMG.K, IMG.label, IMG.SP, SP_changed0, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = local_move(IMG,1000);times(1)=toc;
        tic;[IMG.K, IMG.label, IMG.SP, SP_changed1, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = localonly_move(IMG,500);times(2)=toc;
        if (frame_it>1 && it<5)
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed2, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = merge_move(IMG,1);times(3)=toc;
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed3, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = split_move(IMG,1);times(4)=toc;
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed4, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = switch_move(IMG,1);times(5)=toc;
            IMG.SP_changed = SP_changed0 | SP_changed1 | SP_changed2 | SP_changed3 | SP_changed4;
        else
            IMG.SP_changed = SP_changed0 | SP_changed1;
        end
        disp(times);
        E(end+1) = newE;
        converged = ~any(~arrayfun(@(x)(isempty(x{1})), {IMG.SP(:).N}) & IMG.SP_changed(1:IMG.K));

        if (dispOn)
            sfigure(1);
            im = zeros(size(oim));
            %im(:) = c(IMG.label+1,:);
            imagesc(IMG.label);
            %image(im,'parent',gca);
            title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);

            sfigure(2);
            im = zeros(size(oim,1)+2*IMG.w, size(oim,2)+2*IMG.w, 3);
            im(IMG.w+1:end-IMG.w, IMG.w+1:end-IMG.w, :) = double(oim)/255;
            borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [IMG.xdim IMG.ydim])));
            im = setPixelColors(im, find(borders), [1 0 0]);
            image(im,'parent',gca);
            drawnow;
            %waitforbuttonpress;
            if (frame_it>1)
                %flow = reshape([IMG.SP(1:IMG.prev_K).v], [2 IMG.prev_K]);
                flow = reshape([IMG.SP(1:IMG.prev_K).p_mu] - [IMG.SP(1:IMG.prev_K).p_theta], [2, IMG.prev_K]);
                flow = reshape([IMG.SP(1:IMG.prev_K).v],[2, IMG.prev_K])*sqrt(IMG.hyper.op_Sigma(1));
                u = flow(2,:);
                v = flow(1,:);
                [utotal, vtotal] = flow_total(IMG.prev_K,IMG.xdim,IMG.ydim,IMG.prev_indices,u,v);
                a = flowToColor(cat(3,utotal, vtotal));
                sfigure(3);
                imagesc(a);
                sfigure(4);
                subplot(2,1,1);
                imagesc(utotal);colorbar;
                subplot(2,1,2);
                imagesc(vtotal);colorbar;
            end

            sfigure(5);
            plot(E);
        end
    end
    disp('Done');


    
    SP_UID = {IMG.SP(:).UID};
    mask = arrayfun(@(x)(isempty(x{1})), SP_UID);
    for m = find(mask)
        SP_UID{m} = -1;
    end
    sp_labels(:,:,frame_it) = reshape([SP_UID{IMG.label(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w) +1}], size(oim,1), size(oim,2));
    if (frame_it>1)
        sp_flow{frame_it-1} = flowToColor(cat(3,utotal(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w), vtotal(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w)));
    end

%    save(['results/' name '_' num2str(f, '%02d') '.mat']);
end

% [utotal,vtotal] = SP2_flow(sp_labels);
% a = flowToColor(cat(3,utotal,vtotal));

%%
%view_SP(sp_labels, root, files, frames)
view_sp_gui(sp_labels, root, files, frames,sp_flow)