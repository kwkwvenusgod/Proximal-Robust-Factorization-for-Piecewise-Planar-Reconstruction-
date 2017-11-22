% cluster the means
% sp_means
% sp_indices
% sp_labels

load('seg_try.mat', 'sp_means', 'sp_indices', 'sp_labels');
sp_labels = sp_labels + 1;

K = size(sp_means,1);
X = size(sp_labels,1);
Y = size(sp_labels,2);

mask = false(K,1);

uvals = unique(sp_labels);
mapping = zeros(K,1);
mapping(uvals) = 1:numel(uvals);
mask(uvals) = true;
mask(all(isinf(sp_means(:,:,1)) | isnan(sp_means(:,:,1)), 2)) = false;

sp_means(~mask,:,:) = [];
sp_indices(~mask,:) = [];
sp_labels = mapping(sp_labels);


K = size(sp_means,1);
T = size(sp_means,2);
d = zeros(K);
for i=1:K
    disp(i)
    overlap_ti = ~isinf(sp_means(i,:,4));
    for j=i+1:K
        overlap_t = overlap_ti & ~isinf(sp_means(j,:,4));
        if (any(overlap_t))
            dist = -inf;
            for t=1:T-1
                if (overlap_t(t))
                    possible = overlap_t(t : (min(t+5, T)));
                    dt = find(possible,1,'last') - 1;
                    dsp = mean(sqrt(sum((sp_means(i,t:t+dt,4:5) - sp_means(j,t:t+dt,4:5)).^2,3)), 2);
                    dv = sum((mean(sp_means(i,t:t+dt,6:7),2) - mean(sp_means(j,t:t+dt,6:7),2)).^2,3);
                    varv = min(sum(sum(sp_means(i,t:t+dt,8:9),2),3), sum(sum(sp_means(j,t:t+dt,8:9),2),3));
                    dv = dv / (5*varv);
                    
                    dist = max(dist, dsp*dv);
                end
            end
            d(i,j) = dist;
        end
    end
end


mask = d==-inf;
d(mask) = inf;


W = exp(-0.1*(d+d'));
D = diag(sum(W));
A = D^(-1/2)*(D-W)*D^(-1/2);
[V,L] = eigs(A,20,'SM');
%%

dims = 4;
data = zeros(size(V,1),dims);
for i=1:dims
    data(:,i) = V(:,i) / (abs(L(i,i)));%+0.0001);
end

a = kmeans(data, 20);

for t=1:T-1
    im = zeros(size(sp_labels,1), size(sp_labels,2));
    for k=1:K
        im(sp_indices{k,t}+1) = a(k);
    end
    figure(2);
    imagesc(im);
    title(['Time = ' num2str(t) '   v=' num2str(i)]);
    drawnow;
end

%%
for i=1:20
    for t=1:T-1
        im = zeros(size(sp_labels,1), size(sp_labels,2));
        for k=1:K
            im(sp_indices{k,t}+1) = V(k,i);
        end
        figure(2);
        imagesc(im);
        title(['Time = ' num2str(t) '   v=' num2str(i)]);
        drawnow;
    end
end