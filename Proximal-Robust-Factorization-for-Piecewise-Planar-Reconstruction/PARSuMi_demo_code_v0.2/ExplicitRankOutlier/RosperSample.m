function [Wr,Er,residue_X,i,gk] = RosperSample(Wcap,sample_mask,r,W_0,E_0,gamma)
%ROSPER Summary of this function goes here
%   Wcap: observed data matrix
%   r: rank constraint i.e. rank(W) <= r
%   W_0: initialization for W
%   E_0: initialization for E
%
%   Wr: recovered data matrix

%% Initialization
tau = 1;
E_k = E_0;
W_k = W_0;
one_over_tau = 1/tau;
tol = 1e-6;
max_iter = 10^5;
has_converged = false;
[m,n] = size(Wcap);
verbose_freq = 1*10^2;
residue_X = zeros(max_iter,1);

fprintf('Tolerance: %.2e\n',tol);
fprintf('tau: %.2f\n',tau);
fprintf('max_iter: %.2e\n',max_iter);
fprintf('\n');

%% Main loop
for i=1:max_iter            
    % Optimize over W
    G_W = W_k - one_over_tau*ATAop((W_k + E_k - Wcap),sample_mask);    
    W_kplus1 = RankProjection(G_W,r); % prox(G)==rankproj(G_W,r)
    gk=svd(G_W);
    % Optimize over E
    G_E = E_k - one_over_tau*ATAop((W_kplus1 + E_k - Wcap),sample_mask);     
    beta = gamma/tau;    
    E_kplus1 = L1Minimizer(G_E,beta); %prox(E)==soft thresholding
     
    [has_converged,residue_X(i)] = CheckConvergence(E_kplus1,E_k,W_kplus1,W_k,tol);
    
    if(mod(i,verbose_freq) == 0)
        fprintf('iteration %d\n', i);
        fprintf('Residue of X = %.2e\n',residue_X(i));
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

end

function [flag,residue_X] = CheckConvergence(E_cur,E_prev,W_cur,W_prev,tol)

flag = false;

residue_W = W_cur - W_prev;
residue_E = E_cur - E_prev;
%residue_X = norm([residue_W; residue_E],'fro') / max(1, norm([W_cur;E_cur],'fro'));
residue_X = (norm(residue_W,'fro')+sum(sum(abs(residue_E)))) / max(1, (norm(W_cur,'fro')+sum(sum(abs(E_cur)))));

if(residue_X < tol)
    flag = true;
end

end



