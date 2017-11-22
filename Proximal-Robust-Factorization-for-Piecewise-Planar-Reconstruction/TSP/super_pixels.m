function [label, mu] = super_pixels(oim, K)

IMG = IMG_init(oim,K);
[IMG.K, label, IMG.SP, IMG.SP_changed, IMG.max_UID] = local_move(IMG,1000);

label = prepareSegmentationLabels(label)-1;

if (nargout>1)
    mu_p = {IMG.SP.p_mu};
    mu_a = {IMG.SP.a_mu};
    
    empty_labels = find([IMG.SP.N]==0);
    if (~isempty(empty_labels))
        mu_p(empty_labels) = [];
        mu_a(empty_labels) = [];
    end
    K = numel(mu_a);
    C = numel(mu_a{1});
    mu_a = reshape(cell2mat(mu_a), [C,K]);
    mu_p = reshape(cell2mat(mu_p), [2,K]);
    mu = cat(1, mu_a, mu_p);
end