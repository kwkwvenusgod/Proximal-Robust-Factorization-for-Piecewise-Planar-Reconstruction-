function [Wr,Er,residue_X,residue_fit] = Rosper(Wcap,r,W_0,E_0)
%ROSPER Summary of this function goes here
%   Wcap: observed data matrix
%   r: rank constraint i.e. rank(W) <= r
%   W_0: initialization for W
%   E_0: initialization for E
%
%   Wr: recovered data matrix

%% Initialization
mu_0 = norm(Wcap,'fro');
mu_bar = 10^-3 * mu_0;
mu = mu_0;
nu = 0.7;
tau = 1;
E_k = E_0;
E_k_minus_1 = E_0;
W_k = W_0;
W_k_minus_1 = W_0;
t_k = 1;
t_k_minus_1 = 1;
one_over_tau = 1/tau;
tol = 10^-4;
max_iter = 10^5;
has_converged = false;
[m,n] = size(Wcap);
verbose_freq = 1*10^2;

%% Main loop
for i=1:max_iter
    if(mod(i,verbose_freq) == 0)
        fprintf('iteration %d\n', i);
        fprintf('Residue of X = %f\n',residue_X);
        fprintf('Fitting residue = %f\n',residue_fit);
        fprintf('\n');
    end
    
    t = (t_k_minus_1 - 1)/t_k;
    
    % Optimize over W    
%     YW_k = W_k + t*(W_k - W_k_minus_1);
%     YE_k = E_k + t*(E_k - E_k_minus_1);    
%     G_W = YW_k - one_over_tau*(YW_k + YE_k - Wcap);   
    G_W = W_k - one_over_tau*(W_k + E_k - Wcap);
    W_min = RankProjection(G_W,r);
    
    % Optimize over E
    beta = mu/tau;    
%     YW_k = W_min + t*(W_min - W_k_minus_1);    
%     G_E = YE_k - one_over_tau*(YW_k + YE_k - Wcap);       
    G_E = E_k - one_over_tau*(W_min + E_k - Wcap);
    E_min = L21Minimizer(G_E,beta,m,n);
     
    [has_converged,residue_X,residue_fit] = CheckConvergence(Wcap,E_min,E_k,W_min,W_k,tol);
    if(has_converged)
        break;
    else
        E_k_minus_1 = E_k;
        E_k = E_min;
        W_k_minus_1 = W_k;
        W_k = W_min;
        mu = max(nu*mu, mu_bar);
        t_k_minus_1 = t_k;
        t_k = 0.5*(1 + sqrt(1+4*t_k^2));
    end
end

%% Exited main loop
if(~has_converged)
    fprintf('Not coverged!\n');    
    residue = norm(Wcap - (W_min + E_min),'fro');   
    fprintf('Constraint residue = %.4f\n',residue);        
    fprintf('No. of iterations: %d\n',i);    
end

Wr = W_min;
Er = E_min;

end

function [flag,residue_X,residue_fit] = CheckConvergence(Wcap,E_cur,E_prev,W_cur,W_prev,tol)

flag = false;

residue_W = W_cur - W_prev;
residue_E = E_cur - E_prev;
residue_X = norm([residue_W; residue_E],'fro') / max(1, norm([W_prev;E_prev],'fro'));
numerator = abs(norm( W_cur + E_cur - Wcap,'fro') - norm( W_prev + E_prev - Wcap,'fro'));
denominator = norm(Wcap,'fro');
residue_fit = numerator  / denominator;

if((residue_X < tol) && (residue_fit < (5*tol)))
    flag = true;
end

end



