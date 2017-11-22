function [sp_labels, means, variances] = SP_segmentation(oim, dispOn)
% function [sp_labels] = SP_segmentation(oim)

addpath('flow');
addpath('flow/optical_flow_celiu');
addpath('flow/optical_flow_celiu/mex');

K = round(size(oim,1)*size(oim,2)/60);
K = 100;

IMG = IMG_init(oim,K);
IMG.SP_changed = true(1,IMG.K);
it = 0;
converged = false;

% first split
while (~converged && it<10)
    it = it + 1;
    disp(['Iteration: ' num2str(it)]);

    oldK = IMG.K;
    [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID] = split_move(IMG,1);
    converged = IMG.K - oldK < 2;
    %converged = ~any([IMG.SP.N]>0 & IMG.SP_changed(1:IMG.K));

    if (dispOn)
        sfigure(1);
        subplot(1,1,1);
        im = zeros(size(oim));
        %im(:) = c(IMG.label+1,:);
        imagesc(IMG.label,'parent',gca);
        %image(im,'parent',gca);
        title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);

        sfigure(2);
        subplot(1,1,1);
        im = double(oim)/255;
        borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [size(oim,1), size(oim,2)])));
        im = setPixelColors(im, find(borders), [1 0 0]);
        image(im,'parent',gca);
        drawnow;
        %waitforbuttonpress;
    end
end

IMG.SP_changed(:) = true;
converged = false;
it = 0;
tic;
while (~converged && it<4)
    it = it + 1;
    disp(['Iteration: ' num2str(it)]);

    %[IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID] = split_move(IMG,1);
    %[IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID] = merge_move(IMG,1);
    [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID] = local_move(IMG,50);
    converged = ~any(~arrayfun(@(x)(isempty(x{1})), {IMG.SP(:).N}) & IMG.SP_changed(1:IMG.K));
    %converged = ~any([IMG.SP.N]>0 & IMG.SP_changed(1:IMG.K));


    if (dispOn)
        sfigure(1);
        subplot(1,1,1);
        im = zeros(size(oim));
        %im(:) = c(IMG.label+1,:);
        imagesc(IMG.label,'parent',gca);
        %image(im,'parent',gca);
        title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);

        sfigure(2);
        subplot(1,1,1);
        im = double(oim)/255;
        borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [size(oim,1), size(oim,2)])));
        im = setPixelColors(im, find(borders), [1 0 0]);
        image(im,'parent',gca);
        drawnow;
        %waitforbuttonpress;
    end
end
toc;

emptySPs = [IMG.SP.N]==0;
means = reshape([IMG.SP.a_mu], [3 IMG.K]);
variances = reshape([IMG.SP.a_Sigma], [3 3 IMG.K]);
means(:, emptySPs) = [];
variances(:, :,emptySPs) = [];

sp_labels = prepareSegmentationLabels(double(IMG.label));