% function to warp images with different dimensions
function foo = warpImage(im,vx,vy,xx,yy)

[height,width]=size(vx);

XX=xx+vx;
YY=yy+vy;
mask=XX<1 | XX>width | YY<1 | YY>height;
XX=min(max(XX,1),width);
YY=min(max(YY,1),height);

foo=interp2(xx,yy,im,XX,YY);
foo(mask)=-1;

