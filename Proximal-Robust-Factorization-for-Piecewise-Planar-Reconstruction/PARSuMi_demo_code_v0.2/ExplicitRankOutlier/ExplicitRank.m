function [Wr,Er,residue_X,n_iter] = ExplicitRank(Wcap,sample_mask,r,W_0,E_0,gamma)
%ROSPER Summary of this function goes here
%   Wcap: observed data matrix
%   r: rank constraint i.e. rank(W) <= r
%   W_0: initialization for W
%   E_0: initialization for E
%
%   Wr: recovered data matrix

%% Initialization
tau_default = 1.0;
E_k = E_0;
W_k = W_0;
tol = 1e-6;
max_iter = 10^6;
has_converged = false;
verbose_freq = 1*10^2;
residue_X = zeros(max_iter,1);

fprintf('Tolerance: %.2e\n',tol);
fprintf('max_iter: %.2e\n',max_iter);
fprintf('\n');

%% Main loop
for i=1:max_iter            
%     if (i < 100); 
%        r = 3; %% important
%        tau = 1.1*tau_default; 
%     else
%        r = 4;
%        tau = tau_default;         
%     end    
    tau = tau_default;    
    
    % Optimize over W
    G_W = W_k - (1/tau)*ATAop((W_k + E_k - Wcap),sample_mask);    
    W_kplus1 = RankProjection(G_W,r);
    
    % Optimize over E
    G_E = E_k - (1/tau)*ATAop((W_kplus1 + E_k - Wcap),sample_mask);     
    E_kplus1 = L1Minimizer(G_E,gamma); 
     
    [has_converged,residue_X(i)] = CheckConvergence(E_kplus1,E_k,W_kplus1,W_k,tol);
    
    if(mod(i,verbose_freq) == 0)
        fprintf('iteration %d\n', i);
        fprintf('Residue of X = %.2e\n',residue_X(i));
        fprintf('Fitting residue = %.2e\n',FittingResidue(Wcap,sample_mask,W_k,E_k,tau));
        display_sv = false;
        if(display_sv)
            sv = svd(W_kplus1,'econ');
            fprintf('Singular values of W_kplus1:\n');
            disp(sv);
        end
        fprintf('\n');
    end
    
    if(has_converged) 
        residue_X = residue_X(1:i);
        break;
    else
        E_k = E_kplus1;
        W_k = W_kplus1;
    end
end

%% Exited main loop
if(~has_converged)
    fprintf('Not coverged!\n');  
end        

Wr = W_kplus1;
Er = E_kplus1;
n_iter = i;

end

function res = FittingResidue(Wcap,sample_mask,W,E,tau)

res = norm((1/tau)*ATAop((W+E-Wcap),sample_mask),'fro');

end

function [flag,residue_X] = CheckConvergence(E_cur,E_prev,W_cur,W_prev,tol)

flag = false;

residue_W = W_cur - W_prev;
residue_E = E_cur - E_prev;
residue_X = norm([residue_W; residue_E],'fro') / max(1, norm([W_cur;E_cur],'fro'));

if(residue_X < tol)
    flag = true;
end

end



