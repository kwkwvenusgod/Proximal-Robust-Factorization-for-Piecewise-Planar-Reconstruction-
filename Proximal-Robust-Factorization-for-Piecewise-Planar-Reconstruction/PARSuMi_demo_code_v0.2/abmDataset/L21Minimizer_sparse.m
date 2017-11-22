function E = L21Minimizer_sparse(G,beta,m,n)
%L21MINIMIZER Summary of this function goes here
%This function solves the Moreau-Yoshida regularization
%of the L_{2,1} row-wise mix norm:
%min_E 1/2||E - G||^2 + mu ||E||_{2,1}

E = zeros(m,n);
for i=1:m
    for j=1:n
        g = norm(G(i,j),1);
        if(beta < g)
            E(i,j) = (1-beta/g)*G(i,j);
        end
    end
end

end

