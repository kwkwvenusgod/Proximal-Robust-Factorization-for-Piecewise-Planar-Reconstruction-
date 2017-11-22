function [params] = track_sp_fnc_kp_JC(k_num, param_num,name, root, opticalflow)
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 1));

if (strcmp(name,'cheetah') && k_num==6)
    RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 2));
end

addpath('eval_code/');
% param_num = 1:576

% reset(RandStream.getGlobalStream());
[cov_var_p, cov_var_a, area_var, alpha, beta, deltap_scale] = ndgrid(...
    [1000],...
    [100],...
    [400],...
    [-30,-15,-5],...
    [0.5,1,2,5],...
    [0.001,0.01,0.1,0.33,1]);
    %[0.25,0.5,1,2]);

params.cov_var_p = cov_var_p(param_num);
params.cov_var_a = cov_var_a(param_num);
params.area_var = area_var(param_num);
params.alpha = alpha(param_num);
params.beta = params.alpha+beta(param_num);
params.deltap_scale = deltap_scale(param_num);
params.deltaa_scale = 100;
%kset = [100 200 500 1000 2000 5000 10000, 1e5,5e5,1e6];
%params.kkk=kset(k_num);
kset = [200,700,1200,1700,2200,2700,3200,3700,4200,1000];
kset = [200,300,450,700,1000,1500,2000,2500,3000,3500];
kset = [100,200,300,450,700,1000,1500,2000,2500,3000,3500];
params.K = kset(k_num);
params.Kpercent = 0.8;
return
% params.cov_var_p = 100;
% params.cov_var_a = 10;
% params.area_var = 400;
% params.Kpercent = 0.8;
% params.K = 500;
% params.alpha = -600;
% params.beta = -200;



% params.cov_var_p = 10000;
% params.cov_var_a = 1000;
% params.area_var = 400;
% params.Kpercent = 0.8;
% params.K = 300;
% params.alpha = -600;
% params.beta = -100;
fff= ['outputs_kp/' num2str(k_num*1e5+param_num, '%07d') '/'];
%fff= ['outputs_kp_MiddleburyTest/' num2str((k_num-1)*1e5+param_num, '%07d') '/'];
fff= ['outputs_kp_of/' num2str((k_num-1)*1e5+param_num, '%07d') '/'];
if exist(fff)==0
    mkdir(fff)
end
if exist([fff 'seg_' name '.mat'], 'file')
    disp('Already done!');
    return;
end

%{
params
return;
%}
dispOn = true;
dispOn = false;


        folder = [root name '/'];
        files = [dir([folder '*.jpg']); dir([folder '*.JPG']); dir([folder '*.png']); dir([folder '*.bmp']);dir([folder '*.ppm']);];
        flow_files = dir([folder '*_flow.mat']);
        frames = 1:numel(files);
        oim = imread([folder files(1).name]);
        sp_labels = zeros(size(oim,1), size(oim,2), numel(frames), 'uint64');
        frame_it = 0;

        for f=frames
%         for f=1:5
            disp([name ' : ' num2str(f) ' / ' num2str(numel(frames))]);
            
            if (f==10)
                a=1;
            end

            frame_it = frame_it + 1;
            %disp(['Frame ' num2str(f) ' ----------------------------']);
            oim1 = imread([folder files(f).name]);
            %[folder files(f).name]

            if (frame_it==1)
                IMG = IMG_init(oim1, params);
            else
                % optical flow returns actual x and y flow... flip it
                vx = zeros(size(oim,1), size(oim,2));
                vy = zeros(size(oim,1), size(oim,2));
                % load the optical flow
                load([folder flow_files(f-1).name]);
                %[folder flow_files(f-1).name]
                vx = -flow.bvx;
                vy = -flow.bvy;
                %{
                %}
    %             [vx vy] = coarse_flow(oim, oim1,0);
                IMG = IMG_prop(oim1,vy,vx,IMG);
            end

            oim = oim1;

            if (frame_it>1 && dispOn)
%                 [indices, ~] = populate_indices(double(IMG.prev_K), uint64(IMG.prev_label));
            end

            E = [];
            it = 0;
            IMG.alive_dead_changed = true;
            IMG.SxySyy = [];
            IMG.Sxy = [];
            IMG.Syy = [];
            converged = false;
%             disp('Split/Merge');
            while (~converged && it<5 && true && frame_it==1)
                it = it + 1;

                oldK = IMG.K;
                IMG.SP_changed(:) = true;
                [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy, newE] = split_move(IMG,1);
%                 [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = split_move(IMG,1);
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
%             disp('Local');drawnow;
            while (~converged && it<20)
                it = it + 1;
%                 disp(['Iteration: ' num2str(it)]);drawnow;
                times = zeros(1,5);

                if (opticalflow)
                    tic;[IMG.K, IMG.label, IMG.SP, SP_changed1, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = localonly_move(IMG,1500);times(2)=toc;
                    SP_changed0 = SP_changed1;
                else
                    tic;[IMG.K, IMG.label, IMG.SP, SP_changed0, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = local_move(IMG,1000);times(1)=toc;
                    tic;[IMG.K, IMG.label, IMG.SP, SP_changed1, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = localonly_move(IMG,500);times(2)=toc;
                end
                if (frame_it>1 && it<5)
                    tic;[IMG.K, IMG.label, IMG.SP, SP_changed2, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = merge_move(IMG,1);times(3)=toc;
                    tic;[IMG.K, IMG.label, IMG.SP, SP_changed3, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = split_move(IMG,1);times(4)=toc;
                    tic;[IMG.K, IMG.label, IMG.SP, SP_changed4, IMG.max_UID, IMG.alive_dead_changed, IMG.Sxy,IMG.Syy,IMG.SxySyy,newE] = switch_move(IMG,1);times(5)=toc;
                    IMG.SP_changed = SP_changed0 | SP_changed1 | SP_changed2 | SP_changed3 | SP_changed4;
                else
                    IMG.SP_changed = SP_changed0 | SP_changed1;
                end

%                 tic;[IMG.K, IMG.label, IMG.SP, SP_changed0, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = local_move(IMG,1000);times(1)=toc;
%                 tic;[IMG.K, IMG.label, IMG.SP, SP_changed1, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = localonly_move(IMG,500);times(2)=toc;
%                 if (frame_it>1 && it<5)
%                     tic;[IMG.K, IMG.label, IMG.SP, SP_changed2, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = merge_move(IMG,1);times(3)=toc;
%                     tic;[IMG.K, IMG.label, IMG.SP, SP_changed3, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = split_move(IMG,1);times(4)=toc;
%                     tic;[IMG.K, IMG.label, IMG.SP, SP_changed4, IMG.max_UID, IMG.alive_dead_changed, IMG.SxySyy,IMG.Sxy,IMG.Syy, newE] = switch_move(IMG,1);times(5)=toc;
%                     IMG.SP_changed = SP_changed0 | SP_changed1 | SP_changed2 | SP_changed3 | SP_changed4;
%                 else
%                     IMG.SP_changed = SP_changed0 | SP_changed1;
%                 end
                if (dispOn)
                    disp(times);
                end
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
                    if (frame_it>1 && false)
                        %flow = reshape([IMG.SP(1:IMG.prev_K).v], [2 IMG.prev_K]);
                        flow = reshape([IMG.SP(1:IMG.prev_K).p_mu] - [IMG.SP(1:IMG.prev_K).p_theta], [2, IMG.prev_K]);
                        flow = reshape([IMG.SP(1:IMG.prev_K).v],[2, IMG.prev_K])*sqrt(IMG.hyper.op_Sigma(1));
                        u = flow(2,:);
                        v = flow(1,:);
                        [utotal, vtotal] = flow_total(IMG.prev_K,IMG.xdim,IMG.ydim,indices,u,v);
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
%             disp('Done');

            SP_UID = {IMG.SP(:).UID};
            mask = arrayfun(@(x)(isempty(x{1})), SP_UID);
            for m = find(mask)
                SP_UID{m} = -1;
            end
%            sp_labels(:,:,frame_it) = reshape([SP_UID{IMG.label+1}], size(oim,1), size(oim,2));
            sp_labels(:,:,frame_it) = reshape([SP_UID{IMG.label(IMG.w+1:end-IMG.w,IMG.w+1:end-IMG.w) +1}], size(oim,1), size(oim,2));


%             tempim = zeros(IMG.xdim, IMG.ydim, 3, 'uint8');
%             tempim(:,:,1) = mod(sp_labels(:,:,frame_it), 256);
%             tempim(:,:,2) = mod(floor(sp_labels(:,:,frame_it)/256), 256);
%             tempim(:,:,3) = mod(floor(sp_labels(:,:,frame_it)/65536), 256);
%             imwrite(tempim,['images/' num2str(param_num, '%06d') '/' name '/' num2str(frame_it, '%06d') '.png'])
        end

    %     zip(['images/' num2str(param_num, '%06d') '/' name '.zip'], ['images/' num2str(param_num, '%06d') '/' name]);
    %     tar(['images/' num2str(param_num, '%06d') '/' name '.tar'], ['images/' num2str(param_num, '%06d') '/' name]);
    %
    %     writerObj = VideoWriter('images/temp.avi', 'Archival');
    %     open(writerObj);
    %     for a=1:size(sp_labels,3)
    %         tempim = zeros(IMG.xdim, IMG.ydim, 3, 'uint8');
    %         tempim(:,:,1) = mod(sp_labels(:,:,a), 256);
    %         tempim(:,:,2) = mod(floor(sp_labels(:,:,a)/256), 256);
    %         tempim(:,:,3) = mod(floor(sp_labels(:,:,a)/65536), 256);
    %
    %         frame.cdata = tempim;
    %         frame.colormap = [];
    %         writeVideo(writerObj,frame);
    %     end
    %     close(writerObj);



        sp_labels = uint32(sp_labels);
        save([fff 'seg_' name '.mat'], 'sp_labels', 'folder','params');
    end

