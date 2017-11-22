function [ mask ] = GenerateSampleMask1( m,n,f,skewed )
%GENERATESAMPLEMASK1 Summary of this function goes here
%m and n is dimension of the matric
%f is the rate of missing data, skewed indicates whether the sampling is
%skewed towards diagonal, if not uniform sampling is assumed
%   Detailed explanation goes here

mask=ones(m,n);
M=m*n;
d = floor(f*M);
nf = floor(M);
if(~skewed)

    rand_index = randperm(nf);
    omega_E = sort(rand_index(1:d));

    mask(omega_E)=0;
    mask=logical(mask);

else %if diagonal band shape
    %segment the matrix into N*N block shape
    alpha=m/n;
    map=zeros(m,n);
    %assume this banded structure features uniform decaying
    %use bisection to find a good step size st
    lo=0;
    hi=1;
    
    while abs(hi-lo)>1e-4
        st=(lo+hi)/2;
        sm=0;
        for i=1:m
            for j=1:n
                x=abs(i-j*alpha);
                %sigma=m/2;
                %map(i,j)=1/sqrt(2*pi)/st*exp(-(3*x/m)^2/2/st^2);
                map(i,j)=max(1-st*x^2,0);
                sm=sm+map(i,j);         
            end
        end
        rate=sm/M;
        if rate<1-f %need to decrease st
            hi=st;
        else
            lo=st;
        end
    end
    %map denotes the probability each node got sampled
    mask=map>rand(m,n);  
end


end
