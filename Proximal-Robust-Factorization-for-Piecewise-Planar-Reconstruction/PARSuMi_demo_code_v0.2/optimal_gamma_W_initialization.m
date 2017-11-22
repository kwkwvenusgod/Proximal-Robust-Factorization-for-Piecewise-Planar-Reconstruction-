function [W_nuc,E_nuc,residue_nuc,iter,gamma_W]=optimal_gamma_W_initialization(Wcap,sample_mask,W_0,E_0,r,top,bot)
% This function runs nuclear norm intialization several times with different initialization gamma_W to obtain the best results
% r is desired rank
% top and bottom gives the range of tuning
% W_0 and E_0 areinitializations
% the returns values are the best gamma_W's return value
%top=1;
%bot=0;
[m n]=size(Wcap);
n=max(m,n);
gamma_E = 1/sqrt(n); %regularization weight for sparse error
%gamma_W = 5e-1; %regularization weight for nuclear norm approximation


for i=1:5
	gamma_W = (top+bot)/2;
	[W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(Wcap,sample_mask,W_0,E_0,gamma_W,gamma_E);
	sigma=svds(W_nuc,3*r);
	sig_large=sum(sigma>gamma_W);

	if sig_large > r 
		bot=gamma_W;
    elseif sig_large<r
		top=gamma_W;
    else
        %if i>=3
            break;
        %else
         %   top=gamma_W+1/32;
        %end
	end

end


if i<3 %we might want to further expand the range
    % refine upper bound
    upp_bound=gamma_W;
    low_bound=gamma_W;
    for j=1:2
        gamma_W= (top+upp_bound)/2;
        [W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(Wcap,sample_mask,W_0,E_0,gamma_W,gamma_E);
        sigma=svds(W_nuc,3*r);
        sig_large=sum(sigma>gamma_W);
        if sig_large==r
            upp_bound=gamma_W;
        else
            top=gamma_W;
        end
    end
    
    for j=1:2
        gamma_W= (bot+low_bound)/2;
        [W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(Wcap,sample_mask,W_0,E_0,gamma_W,gamma_E);
        sigma=svds(W_nuc,3*r);
        sig_large=sum(sigma>gamma_W);
        if sig_large==r
            low_bound=gamma_W;
        else
            bot=gamma_W;
        end
    end
    
    gamma_W=low_bound+(-low_bound+upp_bound)*0.35;
    [W_nuc,E_nuc,residue_nuc,iter] = NuclearNormAPG(Wcap,sample_mask,W_0,E_0,gamma_W,gamma_E);
end



end