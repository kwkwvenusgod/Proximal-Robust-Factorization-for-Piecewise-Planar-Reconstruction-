function [Wproj] = RankProjection(G,r)
%RANKPROJECTION Summary of this function goes here
% W: given matrix
% r: desired rank
%
% Wproj: projection of W onto the rank r space

[U,S,V] = svds(G,r);
Wproj = U*S*V';

end

