path_suffix = '';
data_name = '';
data_path = [path_suffix,'\',data_name];
load(data_path);
[patch_similarity, patch_clique] = Construct_Patch_Similarity(uv, uvback, structure_mask, imageref);

 [ Wr_cvx,Er_cvx,N_cvx, Ur_cvx, Vr_cvx]...
     = Plane_Factorization( patch_flow_rank4, mask, r, N0,W_ini, E_ini,focal_length,...
     structure_mask ,masklabel,imageref,patchconf,uv,uvback);