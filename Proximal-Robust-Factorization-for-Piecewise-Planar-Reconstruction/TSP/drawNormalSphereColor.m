clear;clc;
radius = 100;
normalMap = zeros(2*radius+1,2*radius+1,3);
[x y] = meshgrid(1:2*radius+1,1:2*radius+1);
x = x(:);y=y(:);
center = [radius+1;radius+1];
ptnum = length(x);
z = zeros(ptnum,1);
for i=1:length(x)
    sq = radius^2 - (x(i) - center(1))^2 - (y(i) - center(2))^2;
    if sq > 0
        z(i) = sqrt(sq);
        nomraltmp = [(y(i) - center(2)) -(x(i) - center(1)) z(i)];
        normalMap(x(i),y(i),:) = nomraltmp/norm(nomraltmp);
    else
        normalMap(x(i),y(i),:) = [0 0 0];
    end
end
normalMapColor = (normalMap + 1)/2;
figure,imshow(normalMapColor);



% n = 1000;
% [X,Y,Z] = sphere(n);
% surf(X(n/2+1:n+1,:),Y(n/2+1:n+1,:),Z(n/2+1:n+1,:));
% hemiSpX = X(n/2+1:n+1,:);hemiSpY = Y(n/2+1:n+1,:);hemiSpZ = Z(n/2+1:n+1,:);
% hemiSpX = hemiSpX(:);hemiSpY = hemiSpY(:);hemiSpZ = hemiSpZ(:);
% normalMap = zeros(n+1,n+1,3);
% for i = 1:length(hemiSpX)
%     idx = int32((1-hemiSpY(i))*n/2+1);
%     idy = int32((hemiSpX(i)+1)*n/2+1);
%     normalMap(idx,idy,:) = [hemiSpX(i) hemiSpY(i) hemiSpZ(i)];
% end
% normalMap = (normalMap + 1)/2;
% figure,imshow(normalMap);

%surf(X,Y,Z);