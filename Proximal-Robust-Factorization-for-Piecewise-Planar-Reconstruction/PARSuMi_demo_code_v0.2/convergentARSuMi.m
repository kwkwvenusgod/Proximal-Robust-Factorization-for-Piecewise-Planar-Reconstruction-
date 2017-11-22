function [Wr,Er,N]=convergentARSuMi(M,mask,r,N0,N_ini,E_ini,sigma)
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
maxiter=500;
beta1=1/sqrt(max(m,n));
beta2=1/sqrt(max(m,n));
tol=1e-6;

gamma=1/sqrt(max(m,n));
tau=1;
thresh=gamma;

vecM=M(mask);


%zero initialization
Ek=E_ini(mask);
Wk=vecM-Ek;
converged=0;
%[Nk S V]=svds(W0,r);
Nk=N_ini;
Nk=orth(Nk);
% M_list=construct_M_list(Wk,mask);
% [Ai_list,C]=compute_Ai_C(Nk,M_list,mask);
lmax=2;
lambda=1e-8;

while ~converged
   
   [Wk_new, Nk, C, normWk,lambda]=Prox_MatrixCompletion(vecM-Ek, Wk, mask,Nk,beta1,lambda);
   
       if thresh>sigma
            thresh=gamma/tau;
            tau=tau/0.7;
           [Wk_new, C]=huber_fitting_iterative(vecM,mask,Nk,sigma,lmax);
        else
           %[Wk_new, C]=huber_fitting(vecM-Ek,mask,Nk,sigma);
           thresh=0;
       end
   
   [Ek_new, maskE, normEk]=Prox_ErrorCorrection(vecM-Wk_new,Ek,N0,beta2,thresh);
   %note that Wk_new and Ek_new are vectorized
   obj=0.5*norm(vecM-Wk_new-Ek)^2;
   
   
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