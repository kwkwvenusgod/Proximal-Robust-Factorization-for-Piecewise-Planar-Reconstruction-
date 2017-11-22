function [] = track_SP_segmentation(name, orig_root)


params.cov_var_p = 1000;
params.cov_var_a = 1000;
params.area_var = 400;
params.Kpercent = 0.8;
params.K = 800;
params.alpha = -5;
params.beta = -3;
params.deltap_scale = 0.1;
params.deltaa_scale = 100;% 10000;

[deltap_scale, deltaa_scale] = meshgrid([0.1 0.3 1], [1 33 100]);

params.deltap_scale = 0.1;
params.deltaa_scale = 100;
%params.deltap_scale = deltap_scale(paramNum);
%params.deltaa_scale = deltaa_scale(paramNum);


paramname = [num2str(params.deltap_scale) '_' num2str(params.deltaa_scale)];

dispOn = false;
saveOn = true;


%close force all;
% reset(RandStream.getGlobalStream());
% setenv('GUROBI_HOME', '/csail/fisher3/projects/Vseg_sp/C/gurobi500/linux64');
% username = getenv('USER');
% [~,hostname] = system('hostname');
% hostname(numel(hostname)) = [];
% setenv('GRB_LICENSE_FILE', ['/csail/fisher3/licenses/gurobi/' username '/' hostname '/gurobi.lic']);

%mex test_mex.c
%make

addpath('/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_backup/MM_SP/flow/');
addpath('/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_backup/MM_SP/flow/optical_flow_celiu');
addpath('/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_backup/MM_SP/flow/optical_flow_celiu/mex');


%monkeydog
% name = 'parachute';
% name = 'penguin';
% name = 'birdfall2';
% name = 'cheetah';
% name = 'girl';
% name = 'monkeydog';
% name = 'garden';
root = [orig_root name '/'];
%root = ['/csail/fisher3/datasets/SegTrack_201111/' name '/'];
%root = '/csail/fisher3/datasets/SegTrack_201111/parachute/';
%root = '/csail/fisher3/datasets/SegTrack_201111/penguin/';

% name = 'dog';
% name = 'parachute';
% name = 'garden';
% name = 'sample1';
% root = ['/csail/fisher5/datasets/LayerSegmentation_clean/Data/Seg/' name '/'];

% files = dir([root '*.png']);
% files = dir([root '*.bmp']);
files = dir_images(root);
% files = dir([root '*.ppm']);
% files = dir([root '*.jpg']);


% root = ['/afs/csail.mit.edu/u/j/jchang7/code/experiments/Vseg_sp_/table/frames/'];
% files = dir([root '*.jpg']);


% root = ['/csail/fisher3/datasets/garden/'];
% files = dir([root '*.ppm']);

frame_it = 0;
frame_start = 1;
frames = frame_start:1:min(30,numel(files));
% frames = frame_start:1:5;
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

sp_means = -inf(10000, numel(frames), 9);
sp_indices = cell(10000, numel(frames));


flow_files = dir([root '*_flow.mat']);

for f=frames
    frame_it = frame_it + 1;
    disp(['Frame ' num2str(f) ' ----------------------------']);
    oim1 = imread([root files(f).name]);
    oim1 = imresize(oim1, 1);
    

    if (frame_it==1)
        IMG = IMG_init(oim1, params);
        sp_labels_all = zeros(IMG.xdim, IMG.ydim, 'int32');
        sp_flows_all = cell(numel(frames), 1);
        sp_flows_all{end} = zeros(IMG.xdim, IMG.ydim);
    else
        load([root flow_files(f-1).name]);
        %[folder flow_files(f-1).name]
        vx = -flow.bvx;
        vy = -flow.bvy;
        %[vx vy] = coarse_flow(oim, oim1,0);
%         vx(:) = 0;
%         vy(:) = 0;
        IMG = IMG_prop(oim1,vy,vx,IMG);
    end

    oim = oim1;


    if (frame_it>1)
        [indices, ~] = populate_indices(double(IMG.prev_K), IMG.prev_label);
    end

    E = [];
    it = 0;
    IMG.alive_dead_changed = true;
    IMG.SxySyy = [];
    converged = false;
    %disp('Split/Merge');
    while (~converged && it<5 && true && frame_it==1 )
        it = it + 1;
        %disp(['Iteration: ' num2str(it)]);

        oldK = IMG.K;
        IMG.SP_changed(:) = true;
        [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy, newE] = split_move(IMG,1);
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


%     save(['iteration_' num2str(frame_it) '.mat']);
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
    %disp('Local');drawnow;
    while (~converged && it<20)
        it = it + 1;
        %disp(['Iteration: ' num2str(it)]);drawnow;

        times = zeros(1,5);
        tic;[IMG.K, IMG.label, IMG.SP, SP_changed0, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy, newE] = local_move(IMG,1000);times(1)=toc;
        tic;[IMG.K, IMG.label, IMG.SP, SP_changed1, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy, newE] = localonly_move(IMG,500);times(2)=toc;
        if (frame_it>1 && it<5)
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed2, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy, newE] = merge_move(IMG,1);times(3)=toc;
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed3, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy, newE] = split_move(IMG,1);times(4)=toc;
            tic;[IMG.K, IMG.label, IMG.SP, SP_changed4, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy, newE] = switch_move(IMG,1);times(5)=toc;
            IMG.SP_changed = SP_changed0 | SP_changed1 | SP_changed2 | SP_changed3 | SP_changed4;
        else
            IMG.SP_changed = SP_changed0 | SP_changed1;
        end
        %disp(times);
        E(end+1) = newE;

        converged = ~any(~arrayfun(@(x)(isempty(x{1})), {IMG.SP(:).N}) & IMG.SP_changed(1:IMG.K));


        if (true)
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
%             borders = imdilate(borders, [0 1 0; 1 1 1; 0 1 0]);
            im = setPixelColors(im, find(borders), [1 0 0]);
            image(im,'parent',gca);
            drawnow;
            end
            %waitforbuttonpress;
            if (frame_it>1)
                %flow = reshape([IMG.SP(1:IMG.prev_K).v], [2 IMG.prev_K]);
                %flow = reshape([IMG.SP(1:IMG.prev_K).p_mu] - [IMG.SP(1:IMG.prev_K).p_theta], [2, IMG.prev_K]);
                %flow = reshape([IMG.SP(1:IMG.prev_K).p_mu],[2, IMG.prev_K]) - [IMG.prev_pos_mean];
                flow = reshape([IMG.SP(1:IMG.prev_K).v],[2, IMG.prev_K])*sqrt(IMG.hyper.op_Sigma(1));
                u = flow(2,:);
                v = flow(1,:);
                if (any(u(:)==0))
                    a = 1;
                end
                [utotal, vtotal] = flow_total(IMG.prev_K,IMG.xdim,IMG.ydim,indices,u,v);
                if (dispOn)
                a = flowToColor(cat(3,utotal, vtotal));
                sfigure(3);
                imagesc(a);
                sfigure(4);
                subplot(2,1,1);
                imagesc(utotal);colorbar;
                subplot(2,1,2);
                imagesc(vtotal);colorbar;
                
                sfigure(5);
                imagesc(IMG.label>IMG.prev_K)
                end
            end
        end
    end
    disp('Done');
    
    
    SP_UID = {IMG.SP(:).UID};
    mask = arrayfun(@(x)(isempty(x{1})), SP_UID);
    for m = find(mask)
        SP_UID{m} = -1;
    end
    SP_UID_all = [-1 int32([IMG.SP(:).UID])];
    sp_labels(:,:,frame_it) = reshape([SP_UID{IMG.label(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w) +1}], size(oim,1), size(oim,2));
    sp_labels_all(:,:,frame_it) = reshape([SP_UID_all(IMG.label +2)], size(IMG.label,1), size(IMG.label,2));
    if (frame_it>1)
        sp_flow{frame_it-1} = flowToColor(cat(3,utotal(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w), vtotal(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w)));
        sp_flows_all{frame_it-1} = flowToColor(cat(3,utotal, vtotal));
    end

    % get rid of SPs that are only on the outside
    a_mu = reshape([IMG.SP(:).a_mu], 3, numel(IMG.SP))';
    dead = any(a_mu~=a_mu,2);
    
    sp_inds = [IMG.SP(:).UID]+1;
    sp_means(sp_inds, frame_it, 1:3) = reshape([IMG.SP(:).a_mu], 3, numel(IMG.SP))';
    sp_means(sp_inds, frame_it, 4:5) = reshape([IMG.SP(:).p_mu], 2, numel(IMG.SP))';
    
    if (frame_it>1)
        old_SPs = [IMG.SP(:).old];
        old_sp_inds = [IMG.SP([IMG.SP(:).old]).UID]+1;
        sp_means(old_sp_inds, frame_it-1, 6:7) = reshape([IMG.SP(old_SPs).v], 2, numel(IMG.SP(old_SPs)))';
        
        mask = IMG.prev_label>=0;
        utotal = extend_image_simpleIMPORT(double(utotal), mask, IMG.xdim+IMG.ydim);
        vtotal = extend_image_simpleIMPORT(double(vtotal), mask, IMG.xdim+IMG.ydim);

        mean2_utotal = imfilter(utotal.^2, fspecial('gaussian', round(sqrt(IMG.area))*4+1, 0.5*sqrt(IMG.area)));
        mean2_vtotal = imfilter(vtotal.^2, fspecial('gaussian', round(sqrt(IMG.area))*4+1, 0.5*sqrt(IMG.area)));

        mean_utotal = imfilter(utotal, fspecial('gaussian', round(sqrt(IMG.area))*4+1, 0.5*sqrt(IMG.area)));
        mean_vtotal = imfilter(vtotal, fspecial('gaussian', round(sqrt(IMG.area))*4+1, 0.5*sqrt(IMG.area)));

        var_utotal = mean2_utotal - mean_utotal.^2;
        var_vtotal = mean2_vtotal - mean_vtotal.^2;

        varu = find_meansIMPORT(IMG.prev_label, double(max(IMG.prev_label(:))+1), var_utotal);
        varv = find_meansIMPORT(IMG.prev_label, double(max(IMG.prev_label(:))+1), var_vtotal);
        
        sp_means(old_sp_inds, frame_it-1, 8) = varu(:);
        sp_means(old_sp_inds, frame_it-1, 9) = varv(:);
    end
    
    
    [cur_indices, neighbors] = populate_indices(double(IMG.K), IMG.label(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w));
    sp_indices(sp_inds, frame_it) = cellfun(@uint32, {cur_indices.all}, 'UniformOutput', false);
    neighbors(sp_inds, frame_it) = cellfun(@int32, neighbors, 'UniformOutput', false);
    
    if (frame_it>1)
        E(end) = E(end) + logmvnpdf(u, zeros(1,size(IMG.prev_covariance,1)), IMG.prev_covariance) + logmvnpdf(v, zeros(1,size(IMG.prev_covariance,1)), IMG.prev_covariance);
    end
    
%    save(['results/' name '_' num2str(f, '%02d') '.mat']);
end


sp_means(:,:,1) = sp_means(:,:,1) * sqrt(IMG.hyper.oa_Sigma(1));
sp_means(:,:,2) = sp_means(:,:,2) * sqrt(IMG.hyper.oa_Sigma(2));
sp_means(:,:,3) = sp_means(:,:,3) * sqrt(IMG.hyper.oa_Sigma(3));
sp_means(:,:,4) = sp_means(:,:,4) * sqrt(IMG.hyper.op_Sigma(1));
sp_means(:,:,5) = sp_means(:,:,5) * sqrt(IMG.hyper.op_Sigma(2));

last_index = find(any(any(sp_means>=0,3),2), 1, 'last');
sp_means(last_index+1:end, :, :) = [];
sp_indices(last_index+1:end, :) = [];
neighbors(last_index+1:end, :) = [];
% if (saveOn)
%     save(['data_' name '_' paramname '.mat'], 'sp_means', 'sp_indices', 'sp_labels', 'sp_flow', 'neighbors');
%     load(['data_' name '_' paramname '.mat']);
% end
    
empty_indices = all(cellfun(@isempty, sp_indices)==1,2);
labelmap = -ones(size(sp_indices,1),1);
labelmap(~empty_indices) = 0:nnz(~empty_indices)-1;

labelmapint32 = int32(labelmap);
sp_labels = labelmapint32(sp_labels+1);
neighbors(empty_indices,:) = [];
sp_indices(empty_indices,:) = [];
sp_means(empty_indices,:,:) = [];
% for i=1:numel(neighbors)
%     if ~isempty(neighbors{i})
%         neighbors{i} = labelmapint32(neighbors{i}+1);
%     end
% end


T = size(sp_means,2);
neighbors = cell(size(sp_means,1),T);
for t=1:T
    [~, tneighbors] = populate_indices(find(~cellfun(@isempty,sp_indices(:,t)),1,'last'), int32(sp_labels(:,:,t)));
    neighbors(1:numel(tneighbors),t) = tneighbors;
end
neighbors = cellfun(@int32,neighbors,'UniformOutput',false);


neighbors(size(sp_indices,1)+1:end, :) = [];
if (saveOn)
    %save(['data_fixed_' name '_' paramname '.mat'], 'sp_means', 'sp_indices', 'sp_labels', 'sp_flow', 'neighbors');
    save(['output_seg/data_fixed_' name '.mat'], 'sp_means', 'sp_indices', 'sp_labels', 'sp_flow', 'neighbors');
end






% [utotal,vtotal] = SP2_flow(sp_labels);
% a = flowToColor(cat(3,utotal,vtotal));

%%
%view_SP(sp_labels, root, files, frames)
if (dispOn)
    view_sp_gui(sp_labels, root, files, frames,sp_flow)
end
