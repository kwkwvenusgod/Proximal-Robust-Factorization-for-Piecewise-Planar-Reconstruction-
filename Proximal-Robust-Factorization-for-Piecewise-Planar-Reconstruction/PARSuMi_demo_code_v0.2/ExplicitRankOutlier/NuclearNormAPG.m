function [W_bar,E_bar,residue_X,n_iter] = NuclearNormAPG(Wcap,sample_mask,W_0,E_0,gamma_W,gamma_E)
%ROSPER Summary of this function goes here
%   Wcap: observed data matrix
%   W_0: initialization for W
%   E_0: initialization for E
%   gamma_W: regularization weight for W
%   gamma_E: regularization weight for E
%
%   W_bar: recovered data matrix with missing entries
%   filled in and sparse outliers removed
%   E_bar: recovered sparse error matrix

%% Initialization
E_k = E_0;
E_k_minus_1 = E_0;
W_k = W_0;
W_k_minus_1 = W_0;
t_k = 1;
t_k_minus_1 = 1;

tol = 1e-6;
max_iter = 10^5;
has_converged = false;
verbose_freq = 1*10^2;
residue_X = zeros(max_iter,1);

fprintf('Tolerance: %.2e\n',tol);
fprintf('max_iter: %.2e\n',max_iter);
fprintf('\n');


%% Main loop
for i=1:max_iter      
    %Nesterov auxillary point
    t = (t_k_minus_1 - 1)/t_k;    
    Y_W_k = W_k + t*(W_k - W_k_minus_1);
    Y_E_k = E_k + t*(E_k - E_k_minus_1);        
    
    %Gradient evaluated at (Y_W_k,E_k) 
    grad_f_Y_W_k = ATAop((Y_W_k + E_k - Wcap),sample_mask);
    
    % Optimize over W
    G_W = Y_W_k - grad_f_Y_W_k;    
    W_kplus1 = SvThresholding(G_W,gamma_W);    
    
    %Gradient evaluated at (W_k,Y_E_k)
    grad_f_Y_E_k = ATAop((W_k + Y_E_k - Wcap),sample_mask);

    % Optimize over E
    G_E = Y_E_k - grad_f_Y_E_k;     
    E_kplus1 = L1Minimizer(G_E,gamma_E); 
    
    [has_converged,residue_X(i)] = CheckConvergence(Wcap,sample_mask,W_kplus1,E_kplus1,Y_W_k,Y_E_k,tol);    
%     [has_converged,residue_X(i)] = CheckConvergence(E_kplus1,E_k,W_kplus1,W_k,tol);
    
    if(mod(i,verbose_freq) == 0)
        fprintf('iteration %d\n', i);
        fprintf('Residue of X = %.2e\n',residue_X(i));
        fprintf('\n');
    end
    
    if(has_converged) 
        residue_X = residue_X(1:i);    
        
        debug_sv = false;
        if(debug_sv)
            Gsv = svd(G_W);
            figure(5);plot(Gsv,'*');
            title('Singular values of W proximal gradient at convergence');        
        end
        break;
    else
        E_k_minus_1 = E_k;        
        E_k = E_kplus1;
        W_k_minus_1 = W_k;        
        W_k = W_kplus1;
        t_k_minus_1 = t_k;
        t_k = 0.5*(1 + sqrt(1+4*t_k^2));        
    end
end

%% Exited main loop
if(~has_converged)
    fprintf('Not coverged!\n');  
end        

W_bar = W_kplus1;
E_bar = E_kplus1;    
n_iter = i;

end

%Stopping condition based on the Euclidean distance from the subdifferential, as in the APG paper
function [flag,residue] = CheckConvergence(Wcap,sample_mask,W_kplus1,E_kplus1,Y_W_k,Y_E_k,tol)

flag = false;
grad_f_kplus1 = ATAop((W_kplus1 + E_kplus1 - Wcap),sample_mask);
grad_f_Y_k = ATAop((Y_W_k + Y_E_k - Wcap),sample_mask);
S_W = (Y_W_k - W_kplus1) + grad_f_kplus1 - grad_f_Y_k;
S_E = (Y_E_k - E_kplus1) + grad_f_kplus1 - grad_f_Y_k;
S_kplus1 = S_W + S_E;
residue = norm(S_kplus1,'fro') / max(norm([W_kplus1;E_kplus1],'fro'),1);

if(residue < tol)
    flag = true;
end

end



