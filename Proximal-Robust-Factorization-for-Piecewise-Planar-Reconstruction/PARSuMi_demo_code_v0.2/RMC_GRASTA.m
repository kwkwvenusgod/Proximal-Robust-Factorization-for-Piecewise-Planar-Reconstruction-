function [W, E, t, conv]=RMC_GRASTA(M,mask,r, varargin)
tic;
cd grasta.1.2.0

if length(varargin)==1
    OPTIONS.Uinit=varargin{1};
else
    OPTIONS.Uinit=nan;
end

[m,n]=size(mask);
% We set the number of cycles and put the necessary parameters into OPTIONS

maxCycles                   = 100;    % the max cycles of robust mc
OPTIONS.QUIET               = 1;     % suppress the debug information
OPTIONS.MAX_LEVEL           = 20;    % For multi-level step-size,
OPTIONS.MAX_MU              = 15;    % For multi-level step-size
OPTIONS.MIN_MU              = 1;     % For multi-level step-size
OPTIONS.DIM_M               = m;  % your data's ambient dimension
OPTIONS.RANK                = r; % give your estimated rank
OPTIONS.ITER_MIN            = 20;    % the min iteration allowed for ADMM at the beginning
OPTIONS.ITER_MAX            = 20;    % the max iteration allowed for ADMM
OPTIONS.rho                 = 2;   % ADMM penalty parameter for acclerated convergence
OPTIONS.TOL                 = 1e-8;   % ADMM convergence tolerance
OPTIONS.USE_MEX             = 0;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.                                     
CONVERGE_LEVLE              = 20;    % If status.level >= CONVERGE_LEVLE, robust mc converges




index=find(mask>0);
[I,J] = ind2sub([m,n],index);
[J, inxs]=sort(J'); I=I(inxs)';
S=M(index);
% try
    [Usg, Vsg, Osg, conv] = grasta_mc(I,J,S,m,n,maxCycles,CONVERGE_LEVLE,OPTIONS);
% catch
%     fprintf('Run into numerical error\n');
%     Usg=rand(m,r);Vsg=rand(n,r);
%     Osg=zeros(m,n);
%     conv=0;
% end

t=toc;
W=Usg*Vsg';
E=Osg';
cd ..
end