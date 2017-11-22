function E = L0Minimizer(G,en)
%L0Minimizer Summary of this function goes here
%
%
        sorted=sort(abs(G(:)),'descend');              
        E_idx=abs(G)>sorted(en);
        E=zeros(size(G));
        E(E_idx)=G(E_idx);
end
