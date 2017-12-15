function Visualization_3D( image_ref,depth,focal_length)
    %VISUALIZATION_3D Summary of this function goes here
    %   Detailed explanation goes here

    H=size(image_ref,1);
    W=size(image_ref,2);
    [x,y]=meshgrid(1:W,1:H);
    x=x(:);
    y=y(:);
    imsize=H*W;
    ratio_depth=depth(:)/focal_length;
    % comVal=median(depth(:));
    % depth(find(depth>3*comVal))=3*comVal;
    % depth(find(depth<comVal/3))=comVal/3;
    image_ref=double(image_ref)/255;
    im1=image_ref(:,:,1);
    im1= fliplr(im1);
    im2=image_ref(:,:,2);
    im2= fliplr(im2);
    im3=image_ref(:,:,3);
    im3= fliplr(im3);

    im1=im1(:);
    im2=im2(:);
    im3=im3(:);
    X=zeros(imsize,1);
    Y=zeros(imsize,1);
    C=zeros(imsize,3);

    for i=1:imsize
        ratio=ratio_depth(i);
        Xtmp=ratio*(W/2-x(i));
        Ytmp=ratio*(H/2-y(i));

        c1=im1(i);
        c2=im2(i);
        c3=im3(i);
        colortmp=[c1,c2,c3];


        X(i)=Xtmp;
        Y(i)=Ytmp;

        C(i,:)=colortmp;


    end
    campValX=median(abs(X));
    campValY=median(abs(Y));

    figure,scatter3(X,Y,depth(:),[],C,'.')

    axis([-3*campValX,3*campValX,-3*campValY,3*campValY]);
    view([0,-0.0,1]);

end

