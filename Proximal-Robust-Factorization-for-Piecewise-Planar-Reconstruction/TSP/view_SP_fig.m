function [] = view_SP(sp_labels, root, files, frames)

oim = imread([root files(frames(1)).name]);

figure(1);
im = double(oim)/255;

xdim = size(im,1);
ydim = size(im,2);

borders = is_border_valsIMPORT(double(sp_labels(:,:,1)));
im = setPixelColors(im, find(borders), [0 1 0]);
hold off;
image(im);
theTitle = '(s)ingle, (l)ine, (r)ectangle, (p)olygon, (f)reehand';
title(theTitle);
button = 1;
count = 1;
SPs = [];
button = 115;
% p==112
% l==108
% r==114
% s==115
% f==102
while (button~=3 && button~=27)
    if (button==115)
        % single
        button = 1;
        while (button==1)
            [y, x, button] = ginput(1);
            if (button==1 && all(SPs(:)~=sp_labels(round(x), round(y),1)))
                SPs(count) = sp_labels(round(x), round(y), 1);
                im = setPixelColors(im, find(sp_labels(:,:,1)==SPs(count)), [0 1 0]);
                image(im);
                title(theTitle);
                count = count + 1;
            end
        end
    elseif (button==108)
        % line
        image(im);
        points = 1;
        while ~isempty(points)
            [points, button] = imlinesegment();
            
            count_start = count;
            
            N = size(points,1);
            for i=1:N-1
                xdelta = abs(points(i,2) - points(i+1,2));
                ydelta = abs(points(i,1) - points(i+1,1));
                
                nmid = round(max(ydelta,xdelta));
                x = linspace(points(i,2),points(i+1,2),nmid);
                y = linspace(points(i,1),points(i+1,1),nmid);
                
                values = unique(sp_labels(sub2ind([xdim,ydim], round(x),round(y))));
                SPs(count:count+numel(values)-1) = values;
                
                count = count + numel(values);
            end
            
            mask = false(xdim,ydim);
            for j=count_start:count-1;
                mask = mask | sp_labels(:,:,1)==SPs(j);
            end
            im = setPixelColors(im, find(mask), [0 1 0]);
            image(im);
        end
    elseif (button==102 || button==112 || button==114)
        count_start = count;
        if (button==114)
            % rectangle
            h = imrect();
        elseif (button==112)
            % polygon
            h = impoly();
        elseif (button==102)
            % freehand
            h = imfreehand();
        end
        mask = createMask(h);
        values = unique(sp_labels(mask));
        SPs(count:count+numel(values)-1) = values;
        count = count + numel(values);
        delete(h);
        button = 115;
        if (count_start~=count)
            mask = false(xdim,ydim);
            for j=count_start:count-1;
                mask = mask | sp_labels(:,:,1)==SPs(j);
            end
            im = setPixelColors(im, find(mask), [0 1 0]);
            image(im);
        end
    else
        error('Wrong Button!');
    end
    title(theTitle);
end
drawnow;

SPs = unique(SPs);
maxSP = max(max(sp_labels(:,:,1)));

N = numel(SPs);
for f=1:numel(files)
    %borders = false(xdim, ydim);
    borders = is_border_valsIMPORT(double(sp_labels(:,:,f)).*ismember(sp_labels(:,:,f), SPs));
%     for n=1:N
%         borders = borders | is_border_valsIMPORT(double(sp_labels(:,:,f)==SPs(n)));
%     end
    im = im2double(imread([root files(frames(f)).name]));
    im = setPixelColors(im, find(borders), [0 1 0]);
    
    newBorders = is_border_valsIMPORT(double(sp_labels(:,:,f)).*(sp_labels(:,:,f)>maxSP));
    im = setPixelColors(im, find(newBorders & ~borders), [0 0 1]);
    
    image(im);
    title(['Frame ' num2str(f) ' - ' num2str(numel(unique(sp_labels(:,:,f))))]);
    waitforbuttonpress;
end
