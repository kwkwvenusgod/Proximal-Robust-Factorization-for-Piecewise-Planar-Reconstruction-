function [distances, indicesx, indicesy] = generate_distances()

W = 500;

[yg xg] = meshgrid(-W:W, -W:W);

distances = sqrt(xg.^2+yg.^2);
indicesx = int32(xg(:));
indicesy = int32(yg(:));

[distances, i] = sort(distances(:));
indicesx = indicesx(i);
indicesy = indicesy(i);