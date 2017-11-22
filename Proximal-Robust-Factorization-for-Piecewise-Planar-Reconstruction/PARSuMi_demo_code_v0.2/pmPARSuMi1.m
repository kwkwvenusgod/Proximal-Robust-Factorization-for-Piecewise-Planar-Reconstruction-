function [Wr,Er,N]=pmPARSuMi1(M,mask,r,N0,W_ini,E_ini,lmbd,sigma)
% Jan 4, 2014 update
% Partial Majorized PARSuMi with nonzero lambda
%INPUT:
%incomplete matrix: M
%sample mask: mask
%target rank: r
%max corruption: N0
%------------------
%INPUT:(optional)
% W_ini, E_ini: initial W and E (ensure they are feasible!)
% lmbd: small penalty on the unobserved entries as in eqn(1)
% sigma: estimated noise level (for a heuristic)
%------------------
%OUTPUT:
%completed matrix:Wr
%corrected corruption: Er
%subspace basis: N
% ||mask.*(Wr+Er-M)|| should be minimized
%-----------------------------------


% the Tikhonov Regularization term is assumed to be zero here


[m,n]=size(M);
if nargin <8
    sigma=0.001;
end
if nargin<7
    lmbd=1e-3*sigma/m/n*sum(mask(:));
    %this is a guideline to ensure the penelty on the magnitude of missing
    %data is smaller than the main penalty
end
if nargin<6
    E_ini=zeros(m,n);
end
if nargin<5
    [U S V]=svds(M,r);
    W_ini=U*S*V';
end
    
    
    

p=sum(mask(:));
iter=1;
maxiter=500;
beta1=1e-3/sqrt(max(m,n));
beta2=1e-3/sqrt(max(m,n));
tol=1e-6;

gamma=1/sqrt(max(m,n));
tau=1;
thresh=1;

vecM=M(mask);

Ek=E_ini(mask);
[N_ini, S, V]=svds(W_ini,r);
%Wk=M(mask)-Ek;

H=zeros(m,n); H(mask)=1;H(~mask)=sqrt(lmbd);

converged=0;
Nk=N_ini;
C=S*V';
Wk=Nk*C;Wk=Wk(mask);

lmax=3;
lambda=1e-6;


NOMISSING=~sum(mask==0);

while ~converged
   
    matWk=Nk*C;
   [Wk_new1, Nk1,C1,normWk1]=PM_step_H(vecM-Ek, Wk,r,matWk, mask,beta1,lmbd);
   if NOMISSING %no missing data
       flag=0;
   else

       [Wk_new, Nk, C, normWk,lambda]=ProxReg_MatrixCompletion(vecM-Ek, Nk*C, mask,Nk,beta1,lambda, lmbd);
       aa=0.5*norm(vecM-Wk_new-Ek)^2 +0.5*lmbd*norm(~mask.*(Nk*C),'fro')^2+...
           beta1/2*(norm(Wk_new-Wk)^2+lmbd*norm(~mask.*(Nk*C-matWk),'fro')^2);
       bb=0.5*norm(vecM-Wk-Ek)^2+0.5*lmbd*norm(~mask.*(matWk),'fro')^2;
       cc= 0.5*norm(vecM-Wk_new1-Ek)^2+0.5*lmbd*norm(~mask.*(Nk1*C1),'fro')^2+...
           beta1/2*(norm(Wk_new1-Wk)^2+lmbd*norm(~mask.*(Nk1*C1-matWk),'fro')^2);
       fprintf('aa=%f, bb=%f, cc=%f\n',aa,bb,cc);
       flag=(aa<=min(bb,cc));
   end
   if flag 
        if thresh>sigma
            thresh=gamma/tau;
            tau=tau/0.8;
           [Wk_new, C]=huber_fitting_iterative(vecM,mask,Nk,sigma,lmax);
        else
           %[Wk_new, C]=huber_fitting(vecM-Ek,mask,Nk,sigma);
           thresh=0;
        end
   else%use the majorization result
       if ~NOMISSING
        fprintf('[Caution]Majorization results are used.\n');
       end
       Wk_new=Wk_new1;
       Nk=Nk1;
       C=C1;
       normWk=normWk1;
       thresh=0;
   end

   
   [Ek_new, maskE, normEk]=Prox_ErrorCorrection(vecM-Wk_new,Ek,N0,beta2,thresh);
   %note that Wk_new and Ek_new are vectorized
   obj=0.5*norm(vecM-Wk_new-Ek)^2+0.5*lmbd*norm(~mask.*(matWk),'fro')^2;
   e1=norm(Wk_new-Wk)/normWk;
   %e2=norm(Ek_new-Ek)/normEk;
   normWr=norm(C(:));
   
   fprintf('[%d]e1=%.8f, obj=%.5f,||Er||=%.8f, ||Wr||_F=%f.\n',...
       iter,e1,obj,normEk,normWr);
   if e1<tol %&& e2<tol
       converged=1;
       fprintf('Overall converged.\n');
   elseif iter>=maxiter
       converged=1;
       fprintf('Max iterations reached. Quit.\n');
%    elseif sum(maskEk.*maskEk_old)==0
%        converged=1;
%        fprintf('Error support not changed, enter the final iteration.');
%     
   else
       Wk=Wk_new;
       Ek=Ek_new;
       iter=iter+1;
   end
   
   
end

%reconstruct Wr using mapping W(N)
N=Nk;
Wr=Nk*C;
Er=zeros(size(mask));
Er(mask)=Ek_new;

end