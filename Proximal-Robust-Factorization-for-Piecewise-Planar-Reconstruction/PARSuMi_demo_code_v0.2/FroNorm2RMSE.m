function [ RMSE ] = FroNorm2RMSE( Z ,m,n)

RMSE=sqrt(Z.*Z/m/n);


end
