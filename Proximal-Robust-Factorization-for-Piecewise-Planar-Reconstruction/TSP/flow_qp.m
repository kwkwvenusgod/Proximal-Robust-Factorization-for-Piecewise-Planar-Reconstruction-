function f= flow_qp(IMG,IMG_old,prev_covariance)

IMG.SP.N
N_new = [IMG.SP(:).N];
ind_new = find(N_new);
ind_old = find([IMG_old.SP(:).N]);
Syy= prev_covariance(ind_new,ind_new)+1./diag(N_new(ind_new));
Sxy= prev_covariance(ind_old,ind_new);



obs_tmp = mean((reshape([IMG.SP(ind_new).p_mu],2,numel(ind_new))),2) -  mean((reshape([IMG_old.SP(ind_new).p_mu],2,numel(ind_new))),2);
obs_u = obs_tmp(1);
obs_v = obs_tmp(2);

f = Sxy*inv(Syy)*obs_u;

end
