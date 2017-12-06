function [Cons_A, Value_a, Cons_B, Value_b] = Construct_Equivalent_Constraint( n_patches,focal_len )
%% To Construct equivelent constraints  
% 
%  Cons_A * x = 0; Cons_B * x = Cons_d;
% 
%%
n_entry = 36;
%%
n_eq = 5;
Cons_A = zeros(n_patches*n_eq,n_entry*n_patches);
base_A = zeros(n_eq, n_entry);
base_A(1,1) = 1;base_A(1,8) = -1;
base_A(2,13) = 1; base_A(2, 20) = -1;
base_A(3,25) = 1; base_A(3, 32) = -1;
base_A(4, 15) = 1; base_A(4, 33) = -1;
base_A(5, 1) = 1; base_A(5, 15) = focal_len;

for i = 1:n_patches
    Cons_A((i-1)*n_eq+1:i*n_eq, (i-1)*n_entry+1: i*n_entry) = base_A;
end
Value_a = zeros(n_patches*n_eq,1);
%%
Cons_B=sparse(eye(n_entry*n_patches));
diag_ind=1:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;
diag_ind=8:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;
diag_ind=13:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;
diag_ind=15:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;
diag_ind=20:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;
diag_ind=25:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;
diag_ind=32:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;
diag_ind=33:n_entry:n_entry*n_patches;
ind=sub2ind(size(Cons_B),diag_ind,diag_ind);
Cons_B(ind)=0;

b=[zeros(4,1);-focal_len;zeros(4,1);focal_len;zeros(13,1);-1;zeros(5,1);1;zeros(6,1);];
Value_b=repmat(b,n_patches,1);
end

