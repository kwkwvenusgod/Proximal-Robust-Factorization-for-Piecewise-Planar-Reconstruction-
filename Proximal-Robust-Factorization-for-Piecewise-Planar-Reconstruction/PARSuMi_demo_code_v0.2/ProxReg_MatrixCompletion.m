function [Wk_new, Nk, C, normWk, lambda]=ProxReg_MatrixCompletion(M, Wk, mask,N0,beta,lambda,lmbd)
%Note the input Wk is a matrix, while M is a vector.
%Output Wk_New is the vector supported on mask only.

[m,n]=size(mask);
H=zeros(m,n);%the \bar{H} in the paper
H(mask)=sqrt(1+beta);
H(~mask)=sqrt(lmbd*(1+beta));

r=size(N0,2);
converged=0;
max_iter=500;
%lambda=1e-6;
k=10;
times=0;

B=zeros(m,n);
B(mask)=(M+beta*Wk(mask))/sqrt(1+beta);
B(~mask)=lmbd*beta/sqrt(lmbd*(1+beta)+eps)*Wk(~mask);

% M=(M+beta*Wk)/(1+beta);
% M_list=construct_M_list(M,mask);
N_old=N0;
%[Ai_list,C]=compute_Ai_C(N_old,M_list,mask);
e_old=0.5*norm(H.*Wk-B,'fro')^2;
%compute_err(N_old,M,C,mask);

while ~converged
    %[Ai_list,C]=compute_Ai_C(N,M_list,mask);
    [JTJ, Jr]=compute_JTJ_Jr_H(N_old,B, H);   
    
    while(1)

        %imagesc((JTJ+lambda*eye(m*r)));
        delta=pinv((JTJ+lambda*eye(m*r)))*Jr;
%       delta= (JTJ+lambda*eye(m*r))\Jr;
%        delta=(JTJ+lambda*(eye(m*r).*JTJ))\Jr;

        K=reshape(delta,[m,r]); 
        N_new = N_old+K;        
        C=compute_C_H(N_new,B,H);
        e_new=0.5*norm(H.*(N_new*C)-B,'fro')^2; 
       
        if(e_new<e_old ||abs(e_new-e_old)<1e-10)
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
       
    if(1-min(d)<1e-10)     
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
C=compute_C_H(N_new,B,H);
Wk_new=Nk*C;
Wk_new=Wk_new(mask);
normWk=norm(C(:));

end