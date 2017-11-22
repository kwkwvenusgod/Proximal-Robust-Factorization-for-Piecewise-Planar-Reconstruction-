% cluster the means
% sp_means
% sp_indices
% sp_labels

z = randi(3, 20);
z = z(:);
means = rand(3,1,2)*10;
sp_means = randn([numel(z),1,2]) + means(z,:,:);
sp_labels = zeros(50);

load('seg_try.mat', 'sp_means', 'sp_indices', 'sp_labels');

F = size(sp_means,2);
X = size(sp_labels,1);
Y = size(sp_labels,2);

im = zeros(X,Y);

kappa = 0.1;
alpha = 0.0001;

d = [1 2 3 6 7];%1:5;
d = 1:7;
d = 6:7;
D = numel(d);

nu = D+2;

addpath('dpmm');
addpath('dpmm/utilities');
addpath('dpmm/distributions');
addpath('dpmm/dpmix');
addpath('dpmm/hdpmix');
addpath('dpmm/lda');
addpath('dpmm/mixtures');
addpath('dpmm/tests');
addpath('dpmm/datasets');

for f=5
    data = squeeze(sp_means(:,f,d));
    mask = any(data>=0,2) & all(data==data,2);
    data = data(mask,:);
    indices = sp_indices(mask,f);
    Ak = cellfun(@numel, indices);


    N = size(data,1);
    
    % set a random labeling
    z = randi(5, [N,1]);

    for its=1:100
        [z] = mydpmm(data', z);
        % convert the z's to an image
        for i=1:N
            im(indices{i}+1) = z(i);
        end
        imagesc(im);

        title(its);
        drawnow;
    end



    
%     theta = mean(data,1)';
%     Delta = eye(D)*10;
%     Delta = cov(data)/3;
%     
%     
%     N = size(data,1);
%     
%     % set a random labeling
%     z = randi(20, [N,1]);
%     
%     Nk = zeros(1,N);
%     pkappa = zeros(1,N);
%     pnu = zeros(1,N);
%     ptheta = zeros(D,N);
%     pDelta = zeros(D,D,N);
%     pcov = zeros(D,D,N);
%     
%     ncov = (kappa+1)*nu/(kappa*(nu-D-1))*Delta;
%     
%     for k=1:N
%         mask = z==k;
%         Nk(k) = sum(mask(:));
%         pkappa(k) = kappa + Nk(k);
%         pnu(k) = nu + Nk(k);
%         ptheta(:,k) = 1/pkappa(k) * (kappa*theta' + sum(data(mask,:),1));
%         pDelta(:,:,k) = 1/pnu(k) * (nu*Delta + kappa*(theta*theta') - pkappa(k)*(ptheta(:,k)*ptheta(:,k)') + data(mask,:).'*data(mask,:));
%         pDelta(:,:,k) = 1/2*(pDelta(:,:,k) + pDelta(:,:,k)');
%         
%         pcov(:,:,k) = (pkappa(k)+1)*pnu(k) / (pkappa(k)*(pnu(k)-D-1)) * pDelta(:,:,k);
%     end
%     
%     for its=1:100
%         for i=randperm(N)
%             k = z(i);
%             
%             pDelta(:,:,k) = (pDelta(:,:,k) * pnu(k) + pkappa(k)*(ptheta(:,k)*ptheta(:,k)') - data(i,:)'*data(i,:));
%             ptheta(:,k) = (ptheta(:,k) * pkappa(k) - data(i,:)') / (pkappa(k)-1);
%             pDelta(:,:,k) = (pDelta(:,:,k) - (pkappa(k)-1)*(ptheta(:,k)*ptheta(:,k)')) / (pnu(k)-1);
%             
%             pDelta(:,:,k) = 1/2*(pDelta(:,:,k) + pDelta(:,:,k)');
%             
%             Nk(k) = Nk(k)-1;
%             pkappa(k) = pkappa(k) - 1;
%             pnu(k) = pnu(k) - 1;
%             
%             if (Nk(k)==0)
%                 pDelta(:,:,k) = Delta;
%                 ptheta(:,k) = theta;
%             end
%             
%             pcov(:,:,k) = (pkappa(k)+1)*pnu(k) / (pkappa(k)*(pnu(k)-D-1)) * pDelta(:,:,k);
%             
%             Nz = nnz(Nk);
%             probs = zeros(1,Nz+1);
%             possible_ks = -ones(1,Nz+1);
%             possible_ks(1:end-1) = find(Nk>0);
%             count = 1;
%             for k=possible_ks(1:end-1)
%                 probs(count) = Nk(k)*mvnpdf(data(i,:)', ptheta(:,k), pcov(:,:,k));
%                 count = count + 1;
%             end
%             probs(end) = alpha*mvnpdf(data(i,:), theta', ncov);
%             
%             probs = cumsum(probs);
%             
%             k = find(rand*probs(end) <= probs, 1, 'first');
%             if (k==Nz+1)
%                 k = find(Nk==0,1,'first');
%             end
%             z(i) = k;
% 
%             pDelta(:,:,k) = (pDelta(:,:,k) * pnu(k) + pkappa(k)*(ptheta(:,k)*ptheta(:,k)') + data(i,:)'*data(i,:));
%             ptheta(:,k) = (ptheta(:,k) * pkappa(k) + data(i,:)') / (pkappa(k)+1);
%             pDelta(:,:,k) = (pDelta(:,:,k) - (pkappa(k)+1)*(ptheta(:,k)*ptheta(:,k)')) / (pnu(k)+1);
%             
%             pDelta(:,:,k) = 1/2*(pDelta(:,:,k) + pDelta(:,:,k)');
% 
%             Nk(k) = Nk(k)+1;
%             pkappa(k) = pkappa(k) + 1;
%             pnu(k) = pnu(k) + 1;
% 
%             pcov(:,:,k) = (pkappa(k)+1)*pnu(k) / (pkappa(k)*(pnu(k)-D-1)) * pDelta(:,:,k);
%             
%             
%             
% %             mask = z==k;
% %             Nk(k) = sum(mask(:));
% %             pkappa(k) = kappa + Nk(k);
% %             pnu(k) = nu + Nk(k);
% %             ptheta(:,k) = 1/pkappa(k) * (kappa*theta' + sum(data(mask,:),1));
% %             pDelta(:,:,k) = 1/pnu(k) * (nu*Delta + kappa*(theta*theta') - pkappa(k)*(ptheta(:,k)*ptheta(:,k)') + data(mask,:).'*data(mask,:));
% %             pDelta(:,:,k) = 1/2*(pDelta(:,:,k) + pDelta(:,:,k)');
% % 
% %             pcov(:,:,k) = (pkappa(k)+1)*pnu(k) / (pkappa(k)*(pnu(k)-D-1)) * pDelta(:,:,k);
%         end
%         
% %         tz = z - min(z);
% %         scatter(sp_means(:,1,1), sp_means(:,1,2), 20, tz / max(tz));
%         
%         % convert the z's to an image
%         for i=1:N
%             im(indices{i}+1) = z(i);
%         end
%         imagesc(im);
% 
%         title(its);
%         drawnow;
%         
%     end    
end