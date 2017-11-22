function [Wr,Er,N,e1F]=pmPARSuMi(M,mask,r,N0,W_ini,E_ini,sigma)
% Jan 4, 2014 update
% Partial Majorized PARSuMi
%INPUT:
%incomplete matrix: M
%sample mask: mask
%target rank: r
%max corruption: N0
%------------------
%OUTPUT:
%completed matrix:Wr
%corrected corruption: Er
%subspace basis: N
% ||mask.*(Wr+Er-M)|| should be minimized
%-----------------------------------
[m,n]=size(M);
p=sum(mask(:));
iter=1;
maxiter=1000;
beta1=1e-3/sqrt(max(m,n));
beta2=1e-3/sqrt(max(m,n));
tol=1e-6;

gamma=1/sqrt(max(m,n));
tau=1;
thresh=1;

vecM=M(mask);


%zero initialization
Ek=E_ini(mask);
[N_ini, S, V]=svds(W_ini,r);
%Wk=M(mask)-Ek;


converged=0;
%[Nk S V]=svds(W0,r);
Nk=N_ini;
%Nk=orth(Nk);
C=S*V';
Wk=Nk*C;Wk=Wk(mask);
% M_list=construct_M_list(Wk,mask);
% [Ai_list,C]=compute_Ai_C(Nk,M_list,mask);
lmax=3;
lambda=1e-6;

NOMISSING=~sum(mask==0);

while ~converged
   
    
   [Wk_new1, Nk1,C1,normWk1]=PM_step(vecM-Ek, Wk,r,Nk*C, mask,beta1);
   if NOMISSING %no missing data
       flag=0;
   else
       [Wk_new, Nk, C, normWk,lambda]=Prox_MatrixCompletion(vecM-Ek, Wk, mask,Nk,beta1,lambda);
       aa=(0.5*norm(vecM-Wk_new-Ek)^2 + beta1/2*norm(Wk_new-Wk)^2);
       bb=0.5*norm(vecM-Wk-Ek)^2;
       cc= 0.5*norm(vecM-Wk_new1-Ek)^2+beta1/2*norm(Wk_new1-Wk)^2;
       fprintf('aa=%f, bb=%f\n',aa,min(bb,cc));
       flag=(aa<=min(bb,cc));
   end
% flag=0;
% maxiter=1500;
   if flag 
        if thresh>sigma
            thresh=gamma/tau;
            tau=tau/0.7;
%           [Wk_new, C]=huber_fitting_iterative(vecM,mask,Nk,sigma,lmax);
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
   obj=0.5*norm(vecM-Wk_new-Ek)^2;
%    e=norm(reshape(Wk_new,m,n)-M_gt);
%    eF(iter)=e;
   e1=norm(Wk_new-Wk)/normWk;
   e1F(iter)=e1;
   %e2=norm(Ek_new-Ek)/normEk;
   normWr=norm(C(:));
   
   fprintf('[%d]e1=%.8f, obj=%.5f,||Er||=%.8f, ||Wr||_F=%f.\n',...
       iter,e1,obj,normEk,normWr);
   if e1<tol %&& e2<tol
       converged=0;
       fprintf('Overall converged.\n');
       break;
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