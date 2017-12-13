data_name = 'chessboardforward';
data_path = ['.\data\', data_name,'\'];
focal_length = 356.2500;
%% process superpixel segments
n_seg = 300;
seg_path = [data_path, 'labels_',data_name ,'_K', num2str(n_seg), '.mat'];
segs = load(seg_path);
segs = segs.sp_labels;
consistent_label_segs = Find_Correspondences(segs);
label_map = segs(:,:,1);
label_ref = consistent_label_segs(:,:,1);
%% load flow
flowl2r_path = [data_path, 'flow_', data_name, 'l2r.mat'];
flowl2r = load(flowl2r_path);
flowl2r = flowl2r.flowl2r;
flowr2l_path = [data_path, 'flow_', data_name, 'r2l.mat'];
flowr2l = load(flowr2l_path);
flowr2l = flowr2l.flowr2l;
%% flow consistency 
frame_select = 9;
uv = flowl2r{frame_select};
uvback = flowr2l{frame_select};
%%
imagetype = 'ppm';
data_files = dir([data_path '*.',imagetype]);
image_ref = imread(fullfile(data_path, data_files(1).name));
image_set = cell(1,size(data_files,1));
for i = 1:size(data_files,1)
    image_set{i} = imread(fullfile(data_path, data_files(i).name));
end
%% scene constraint
patch_similarity = Construct_Patch_Similarity(uv, uvback, label_ref,image_ref);
%% affine estimation
[ patch_flow, patch_flow_rank4,label_set,label_map_consistent,residual] = Affine_Estimation(flowl2r,label_ref,patch_similarity);
%% flow confidence
conang=tan(18/180*pi);
fovx=ceil(conang*focal_length)*2;
fovy=ceil(conang*focal_length)*2;
patch_confidence_map = Calc_Flow_Confidence_in_Patch( flowl2r,flowr2l,image_ref,image_set,label_map_consistent,fovx,fovy );
%% good initialization
n_frame = size(patch_flow_rank4,1);
n_patch = size(patch_flow_rank4,2)/6;
fprintf('Getting a good initial estimate by using convex nuclear norm function in place of rank\n');
E0=zeros(size(patch_flow_rank4));
W0=patch_flow_rank4-E0;
mask = true(size(patch_flow_rank4));
top=1;bot=0;%define range we want to tune gamma_W
r = 6;% rank
[W_nuc,E_nuc,residue_nuc,iter,gamma_W]=optimal_gamma_W_initialization(patch_flow_rank4,mask,W0,E0,r,top,bot);
W_ini = W_nuc;
E_ini=E_nuc;
N0 = ceil(size(W_ini(:),1)*0.6);
%% optimization
[ Wr_cvx,Er_cvx,N_cvx, Ur_cvx, Vr_cvx]...
    = Plane_Factorization( patch_flow_rank4, mask, r, N0,W_ini, E_ini ,focal_length, label_map_consistent, label_set , label_ref ,image_ref ,patch_confidence_map,uv,uvback );
%% visualization
[ Depth_map ] = Depth_Completion( label_map, label_map_consistent, Vr_cvx, focal_length  );


%% load ground truth 
gt_path = ['.\groundtruth\', data_name];
translation_gt_path = [gt_path, '\Translation_Para.mat'];
rotation_gt_path = [gt_path, '\Rotation_Para.mat'];
translations_gt=load(translation_gt_path);
rotations_gt=load(rotation_gt_path);
translations_gt=struct2cell(translations_gt);
rotations_gt=struct2cell(rotations_gt);
[motion_gt, normal_gt] = Calc_GT(translations_gt, rotations_gt);
%% Evaluation motion
[ translation_error, translation_error_average, rotation_angular_error,rotation_angular_error_average, rotation_epe, rotation_epe_average] = Evaluation_Motion( motion_gt, Ur_cvx );
%% Evaluation normal
mask_path = ['.\data\', data_name,'\mask.jpg'];
mask_checkerboard = imread(mask_path);
normal_map = Calc_Normal(Vr_cvx, focal_length);
checkerboard_normal =  Calc_Checkerboard_Normal( normal_map,mask_checkerboard, label_map_consistent, label_set );
checkerboard_normal_erro = Evaluate_Normal( normal_gt, checkerboard_normal);
