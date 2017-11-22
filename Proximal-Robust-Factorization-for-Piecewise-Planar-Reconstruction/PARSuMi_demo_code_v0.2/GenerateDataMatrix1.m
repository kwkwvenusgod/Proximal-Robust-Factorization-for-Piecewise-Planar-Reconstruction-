function [Wcap,W,E,N,omega_E] = GenerateDataMatrix1(m,n,r,f,sigma, uni)
%GENERATEDATAMATRIX Summary of this function goes here
% (m,n): size of data matrix
% r: true rank of data matrix
% f: in (0,1], fraction of the sparse error
% 
% Wcap: data matrix corrupted with sparse error
% W: clean data matrix
% E: row-wise error

%% Generate data matrix W of rank r
W_left = rand(m,r);%-0.5;
W_right = rand(r,n);%-0.5;
W = W_left * W_right; %W should be rank r

%% Generate sparse error E, with fraction f of
%  the number of entries corrupted by uniform noise
M=m*n;
d = floor(f*M);
nf = floor(M);
rand_index = randperm(nf);
omega_E = sort(rand_index(1:d));
E = zeros(m,n);
n_error = size(omega_E,2);
for i=1:n_error
    E(omega_E(i)) = rand(1)*(uni(2)-uni(1))+uni(1);
end
%% Generate dense gaussian noise term N
N=randn(m,n)*sigma;
%% Generate observed data matrix Wcap
Wcap = W + E + N;


