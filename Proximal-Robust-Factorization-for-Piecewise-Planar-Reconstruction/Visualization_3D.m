function Visualization_3D( label_map,label_set,depth,image_ref,focal_length )
    %VISUALIZATION_3D Summary of this function goes here
    %   Detailed explanation goes here
    count=1;
    H=size(label_map,1);
    W=size(label_map,2);
    im=double(image_ref)/255;
    im1=im(:,:,1);
    im2=im(:,:,2);
    im3=im(:,:,3);
    X=0;
    Y=0;
    Z=0;
    C=[0,0,0];

    for i=1:length(label_set)
        ll=find(label_map==label_set(i));
        [ytmp,xtmp]=find(label_map==label_set(i));
        ytmp=H/2*ones(length(ll),1)-ytmp;
        xtmp=W/2*ones(length(ll),1)-xtmp;
        Ztmp=depth(ll);

        c1=im1(ll);
        c2=im2(ll);
        c3=im3(ll);
        colortmp=[c1,c2,c3];
        Xtmp=xtmp.*Ztmp/focal_length;
        Ytmp=ytmp.*Ztmp/focal_length;

        X=[X;-Xtmp];
        Y=[Y;Ytmp];
        Z=[Z;Ztmp];
        C=[C;colortmp];
        %   scatter3(-Xtmp,Ytmp,Ztmp,[],color,'*') ;hold on
        count=count+length(ll);
    end
    figure,scatter3(X,Y,Z,[],C,'*')
    view([0,0,1]);
    hold off
end

