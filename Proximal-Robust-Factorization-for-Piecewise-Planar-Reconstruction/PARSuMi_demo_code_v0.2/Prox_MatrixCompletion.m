function [Wk_new, Nk, C, normWk, lambda]=Prox_MatrixCompletion(M, Wk, mask,N0,beta,lambda)

[m,n]=size(mask);
r=size(N0,2);
converged=0;
max_iter=500;
lambda=1e-6;
k=10;
times=0;

M=(M+beta*Wk)/(1+beta);
M_list=construct_M_list(M,mask);
N_old=N0;
[C]=compute_C(N_old,M_list,mask);
e_old=compute_err(N_old,M,C,mask);

while ~converged
    %[Ai_list,C]=compute_Ai_C(N,M_list,mask);
    [JTJ Jr]=compute_JTJ_Jr(N_old,M_list, mask);   
    
    
    while(1)

        delta=pinv((JTJ+lambda*eye(m*r)))*Jr;
%       delta= (JTJ+lambda*eye(m*r))\Jr;
%        delta=(JTJ+lambda*(eye(m*r).*JTJ))\Jr;

        K=reshape(delta,[m,r]); 
        N_new = orth(N_old+K);
        [C]=compute_C(N_new,M_list,mask);
        e_new=compute_err(N_new,M,C,mask);    
       
        if(e_new<e_old ||abs(e_new-e_old)<1e-8)
            lambda=lambda/k;
            break;
        else
            lambda=lambda*k;
        end
        
    end
    
    
    
    %check convergence over N
    N_new=orth(N_new);
    d=svd(N_old'*N_new);
    %if(e_old-e_new>=0&e_old-e_new<1e-10)

    
%    e(times+2)=e_new;
    N_old=N_new;
    e_old=e_new;
    times=times+1;   
       
    if(1-min(d)<1e-8)     
        converged=1;
        %fprintf('Converged.\n')
    elseif times>max_iter
        converged=1;
        %fprintf('Max Iteration reached for MC step. Quit.\n');
    else
        fprintf('[N][%d] Convergence residual = %.6f\n',times,min(d)); 
    end
end


Nk=N_new;
[C]=compute_C(Nk,M_list,mask);
Wk_new=Nk*C;
Wk_new=Wk_new(mask);
normWk=norm(C(:));

end