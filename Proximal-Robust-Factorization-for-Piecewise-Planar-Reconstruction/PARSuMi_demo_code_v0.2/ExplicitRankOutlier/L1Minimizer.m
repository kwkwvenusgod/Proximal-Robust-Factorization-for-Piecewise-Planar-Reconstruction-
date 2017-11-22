function E = L1Minimizer(G,beta)
%L1Minimizer Summary of this function goes here
%This function solves  min_E 1/2||E - G||^2 + beta*||E||_1
%by the soft thresholding operator

E = max(G-beta,0);
E = E + min(G+beta,0);

end

