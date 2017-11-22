function  plot_traj( Wr )
%PLOT_TRAJ Summary of this function goes here
%   Detailed explanation goes here
[m n]=size(Wr);
hold off;
plot(Wr(1:2:m,1:n),Wr(2:2:m,1:n)) % plots all tracks
axis ij  % y in -ve direction


end
