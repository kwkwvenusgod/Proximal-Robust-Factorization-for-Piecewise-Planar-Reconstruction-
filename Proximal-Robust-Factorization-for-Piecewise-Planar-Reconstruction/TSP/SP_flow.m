
if (false)
    %s = RandStream('mt19937ar','Seed','shuffle');
    %RandStream.setGlobalStream(s);
    reset(RandStream.getGlobalStream());

    root = '/csail/fisher3/datasets/SegTrack_201111/monkeydog/';
    files = dir([root '*.bmp']);
    root = '/csail/fisher3/datasets/Brox_motion/';
    files = dir([root '*.jpg']);
    %root = '/csail/fisher3/datasets/SegTrack_201111/penguin/';
    

    K = 1000;
    start = 28;

    oim = imread([root files(start).name]);
    sp_labels = zeros(size(oim,1), size(oim,2), 'uint64');

    IMG = IMG_init(oim,K);

    it = 0;
    converged = false;
    while (~converged && it<20)
        it = it + 1;
        disp(['Iteration: ' num2str(it)]);

        %[IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID] = split_move(IMG,1);
        %[IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID] = merge_move(IMG,1);
        [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID] = local_move(IMG,1);
        converged = ~any(~arrayfun(@(x)(isempty(x{1})), {IMG.SP(:).N}) & IMG.SP_changed(1:IMG.K));
        %converged = ~any([IMG.SP.N]>0 & IMG.SP_changed(1:IMG.K));



        sfigure(1);
        im = zeros(size(oim));
        %im(:) = c(IMG.label+1,:);
        imagesc(IMG.label,'parent',gca);
        %image(im,'parent',gca);
        title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);

        sfigure(2);
        im = double(oim)/255;
        borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [size(oim,1), size(oim,2)])));
        im = setPixelColors(im, find(borders), [1 0 0]);
        image(im,'parent',gca);
        drawnow;
        %waitforbuttonpress;
    end
    
    label = IMG.label;
    oim1 = imread([root files(start).name]);
    oim2 = imread([root files(start+1).name]);
    
    save('tennis', 'label', 'oim1', 'oim2');
end

clear

load monkey_label.mat
load tennis.mat
oim1 = double(oim1);
oim2 = double(oim2);

%oim2 = oim1(10:109,1:100,:);
%oim1 = oim1(1:200,1:200,:);
%oim2 = oim2(1:200,1:200,:);
%label = label(1:200,1:200);

label = prepareSegmentationLabels(label);


% oim1 = zeros(100)+randn(100)*10;
% oim1(20:40,20:40) = 128;
% oim2 = zeros(100)+randn(100)*10;
% oim2(25:45,30:50) = 128;
% f = fspecial('gaussian', 20,5);
% oim1 = imfilter(oim1, f);
% oim2 = imfilter(oim2, f);
% figure(1);
% imagesc(oim1/255);
% figure(2);
% imagesc(oim2/255);

X = size(oim1,1);
Y = size(oim1,2);
C = size(oim1,3);

% [yg xg] = meshgrid(1:X, 1:Y);
% 
% label = floor((xg-1)/20) + floor((yg-1)/20)*X/20;
% label = prepareSegmentationLabels(label);

%% stuff stored in label and oim now
unique_vals = unique(label)';

K = max(unique_vals);





%% find the mean stats
mu1 = zeros(C,K);
mu2 = zeros(X,Y,C, K);
xi = cell(K,1);
yi = cell(K,1);
Ni = zeros(K,1);
parfor k=unique_vals
    mask = (label==k);
    [xi{k},yi{k}] = find(mask);
    Ni(k) = numel(xi{k});
    indices = find(mask);
    for c=1:C
        mu1(c,k) = mean(oim1(indices));
        indices = indices + X*Y;
    end
    %mu2(:,:,:,k) = imfilter(oim2, mask/sum(mask(:)));
end

%% find the neighboring super pixels
neighbors = cell(K,1);
for k=unique_vals
    neighbors{k} = find_neighbor_SPs(label, k);
end



%% main grad descent loop
lambda = 10;
lambda = 100;
u = zeros(K,1);
v = zeros(K,1);

xoff = ceil(X/2);
yoff = ceil(Y/2);

for it=1:10
    
    tmu2 = nan(C,K);
    
    dmu_du = nan(C,K);
    dmu_dv = nan(C,K);
    dmu_du2 = nan(C,K);
    dmu_dv2 = nan(C,K);
    dmu_dudv = nan(C,K);
    
    G = zeros(2*K,1);
    H = zeros(2*K, 2*K);
    
    for k=1:K
    %for k=randi(K)
        x = xi{k} + round(u(k));
        y = yi{k} + round(v(k));
        mask = (x>0 & x<=X & y>0 & y<=Y);
        
        if (any(mask))
            x = x(mask);
            y = y(mask);
            indices = sub2ind([X,Y],x,y);

            indicesf = indices(x<X)+1;
            indicesb = indices(x>1)-1;
            indicesd = indices(y<Y)+X;
            indicesu = indices(y>1)-X;
            indicesfd = indices(x<X & y<Y)+X+1;
            indicesfu = indices(x<X & y>1)-X+1;
            indicesbd = indices(x>1 & y<Y)+X-1;
            indicesbu = indices(x>1 & y>1)-X-1;

            for c=1:C
                offset = X*Y*(c-1);

                mu = mean(oim2(indices+offset));
                muf = mean(oim2(indicesf+offset));
                mub = mean(oim2(indicesb+offset));
                mud = mean(oim2(indicesd+offset));
                muu = mean(oim2(indicesu+offset));
                mufd = mean(oim2(indicesfd+offset));
                mufu = mean(oim2(indicesfu+offset));
                mubd = mean(oim2(indicesbd+offset));
                mubu = mean(oim2(indicesbu+offset));

                tmu2(c,k) = mu;

                % 1st order central difference
                dmu_du(c,k) = (muf-mub) / 2;
                dmu_dv(c,k) = (mud-muu) / 2;

                % 2nd order central difference
                dmu_du2(c,k) = (muf-2*mu+mub) / 2;
                dmu_dv2(c,k) = (mud-2*mu+muu) / 2;

                % 2nd order 2var central difference
                dmu_dudv(c,k) = (mufd + mubu - mufu - mubd) / 4;
            end
        end
        
        
        
        
        
%         x = xoff + round(u(k));
%         y = yoff + round(v(k));
%         
%         tmu2(:,k) = reshape(mu2(x,y,:,k), [C,1]);
%         
%         % 1st order central difference
%         dmu_du(:,k) = (mu2(x+1, y, :,k) - mu2(x-1, y, :,k)) / 2;
%         dmu_dv(:,k) = (mu2(x, y+1, :,k) - mu2(x, y-1, :,k)) / 2;
%         
%         % 2nd order central difference
%         dmu_du2(:,k) = (mu2(x+1,y,:,k) - 2*mu2(x,y,:,k) + mu2(x-1,y,:,k)) / 2;
%         dmu_dv2(:,k) = (mu2(x,y+1,:,k) - 2*mu2(x,y,:,k) + mu2(x,y-1,:,k)) / 2;
%         
%         % 2nd order 2var central difference
%         dmu_dudv(:,k) = (mu2(x+1,y+1,:,k) + mu2(x-1,y-1,:,k) - mu2(x+1,y-1,:,k) - mu2(x-1,y+1,:,k)) / 4;
        
        % populate the gradient vector. stores u and then v
        G(k)   = 2*Ni(k)* sum( ((mu1(:,k) - tmu2(:,k)) .* (-dmu_du(:,k))) ) + 2*lambda*sum(u(k)-u(neighbors{k}));
        G(k+K) = 2*Ni(k)* sum( ((mu1(:,k) - tmu2(:,k)) .* (-dmu_dv(:,k))) ) + 2*lambda*sum(v(k)-v(neighbors{k}));
        
        if (isnan(G(k)));   G(k)   = 2*lambda*sum(u(k) - u(neighbors{k})); end;
        if (isnan(G(k+K))); G(k+K) = 2*lambda*sum(v(k) - v(neighbors{k})); end;
        
        % populate the hessian matrix. stores u and then v
        H(k,k)     = 2*Ni(k)*sum( ((mu1(:,k) - tmu2(:,k)) .* (-dmu_du2(:,k))) + (dmu_du(:,k)).^2 ) + 2*lambda*numel(neighbors{k});
        H(k+K,k+K) = 2*Ni(k)*sum( ((mu1(:,k) - tmu2(:,k)) .* (-dmu_dv2(:,k))) + (dmu_dv(:,k)).^2 ) + 2*lambda*numel(neighbors{k});
        H(k,k+K)   = 2*Ni(k)*sum( ((mu1(:,k) - tmu2(:,k)) .* (-dmu_dudv(:,k))) + (dmu_dv(:,k).*dmu_du(:,k)) );
        
        if (isnan(H(k,k)))
            H(k,k) = 2*lambda*numel(neighbors{k});
        end
        if (isnan(H(k+K,k+K)))
            H(k+K,k+K) = 2*lambda*numel(neighbors{k});
        end
        if (isnan(H(k,k+K)))
            H(k,k+K) = 0;
        end
        
        H(k+K,k)   = H(k,k+K);
        
        
        
        for n=neighbors{k}
            H(k,n) = -2*lambda;
            H(k+K,n+K) = -2*lambda;
        end
    end

    %new_uv = [u;v] - max(min(H^-1 * G,20),-20);
    %new_uv = [u;v] - 0.1*H^-1 * G;
    %new_uv = [u;v] - 0.0001*G;
    
    delta = 0.001*G;
    %delta = 0.5*H^-1 * G;
    max_delta = max(abs(delta(:)));
    if max_delta>10
        delta = delta / max_delta * 10;
    end
    new_uv = [u;v] - delta;
    
    u = new_uv(1:K);
    v = new_uv(K+1:end);
    
    utotal = zeros(X,Y);
    vtotal = zeros(X,Y);
    for k=1:K
        mask = label==k;
        utotal(mask) = u(k);
        vtotal(mask) = v(k);
    end
%     [wim] = myWarpFL_flip(oim1, cat(3,vtotal,utotal));
    
    [wim] = myWarpFL(oim2, cat(3,vtotal,utotal));
    
    sfigure(1);
    imagesc(wim/255);
    title(['Warped ' num2str(it)]);
    %waitforbuttonpress
    sfigure(2);
    imagesc(oim1/255);
    title(['True ' num2str(it)]);
    %waitforbuttonpress;
    drawnow;
end

disp('DONE!');
figure(3);
imagesc(wim/255);
waitforbuttonpress;
imagesc(oim1/255);
waitforbuttonpress;
imagesc(wim/255);
waitforbuttonpress
drawnow;

a=flowToColor(cat(3,utotal, vtotal));