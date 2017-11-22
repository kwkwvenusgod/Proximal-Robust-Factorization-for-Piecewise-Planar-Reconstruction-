function [x]=Gauss_Seidel(A,b,tol,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x,xx]=gauss_seidel(A,b,tol,[x0])
% This function should work in both Octave and Matlab
% Octave is free software available for download at 
% http://www.gnu.org/software/octave/download.html
%
% For a given A and b performs Gauss-Seidel iteration up to a tolerance tol 
% to find x from the system.
%
% Ax=b
% 
% INPUTS:
% A  [NxN]: Matrix on which to perform Gauss-Seidel iteration.
% b  [Nx1]: Right hand side vector.
% tol[1x1]: Relative tolerance in the L_infinity norm after which iteration will stop.
% x0 [Nx1]: (optional) Initial guess vector. If none is supplied x0=0 is used.
%
% WARNING: 
% There is no error checking. All input is assumed to be
% correct. You will get undefined behaviour (i.e., cryptic error
% messages or incorrect results) if your input is not within the
% parameters specified above. In particular, you will run into trouble
% if your A,b and x0 have inconsistent sizes.
%
% OUTPUT:
% xx [NxM]: Sequence of solutions between initial guess and final solution inclusive
%           its final width M will vary depending on the number of iterations.
%
% DIAGNOSTIC INFO:
% On the (j-1)th iteration, the following will be output to the terminal:
% [j] [max(abs(xx(:,j)-xx(:,j-1)))/max(abs(xx(:,j)))]
% 
% EXAMPLE:
% A=[10 -1 2 0;-1 11 -1 3;2  -1 10 -1;0 3 -1 8]
% b=[6 25 -11 15]
% tol=1e-3
% [xx]=gauss_seidel(A,b,tol);
%
% Kevin Mitchell kevmitch@math.sfu.ca 2009.06.09
% Please let me know if you have trouble running this code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% force b to be a column vector and find its length
b=b(:);
N=length(b);
% nargin tells you how many arguments were used to call the current function
if nargin<4
    %initialise x to a zero vector of the same length as b by default
    x=zeros(N,1);
else
    % force x to be a column vector if it was given as an argument
    x=x(:);
end

% initialise output array
xx=x;
% ensure that the initial delta is greater than the tolerance so that we
% perform at least one iteration
delta=tol*2;

% %The clear way
% while tol< delta
%     for ii=1:N
%         % update x one element at a time so that we can use updated values 
%         % note we skip the diagonal of A and put it in the denominator
%                   %new values              %old values 
%         x(ii)=( - A(ii,1:ii-1)*x(1:ii-1) - A(ii,ii+1:end)*x(ii+1:end) + b(ii))/A(ii,ii);
%     end
%     % compute the relative L_infinity change from last
%     delta=max(abs(x-xx(:,end)))/max(abs(x));
%     % extend the output array by adding the current vector onto it
%     xx(:,end+1)=x;
%     % print diagnostic info
%     fprintf('%20d %20g\n',size(xx,2),delta);
% end

%The clever way (this should be faster for larger matrices)
% extract diagonals from A
d=diag(A);
%D=diag(diag(A));
% sine flip A, and remove the diagonals
% then divide by the diagonals
T=bsxfun(@rdivide,(-A+diag(d)),d);
% divide b by the diagonals
c=b./d;
  while tol< delta
%     for ii=1:N
%         %update x one element at a time so we can make use of old values
%         x(ii)=T(ii,:)*x+c(ii);
%     end
    x=T*x+c;
    % compute the relative L_infinity change from last
    delta=max(abs(x-xx(:,end)))/max(abs(x));
    % extend the output array by adding the current vector onto it
    xx(:,end+1)=x;
    % print diagnostic info
    %fprintf('%20d %20g\n',size(xx,2),delta);
end