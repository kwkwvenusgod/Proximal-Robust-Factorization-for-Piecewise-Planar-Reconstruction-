% root = '/afs/csail.mit.edu/u/j/jchang7/Desktop/Latex/jchang7_2013_cvpr_sp/figures/images/';
% gbh = '/csail/fisher5/projects/Vseg_sp/saves/girl_gbh/';
% swa = '/csail/fisher5/projects/Vseg_sp/saves/girl_swa/';
% tsp = '/csail/fisher5/projects/Vseg_sp/saves/girl_tsp/';
% 
% count = 1;
% for i=0:2:8
%     copyfile([gbh num2str(i,'%02d') '.png'], [root 'girl_gbh_' num2str(count) '.png']);
%     copyfile([swa num2str(i,'%02d') '.png'], [root 'girl_swa_' num2str(count) '.png']);
%     copyfile([tsp num2str(i,'%02d') '.png'], [root 'girl_tsp_' num2str(count) '.png']);
%     count = count + 1;
% end

%%
root = '/afs/csail.mit.edu/u/j/jchang7/Desktop/Latex/jchang7_2013_cvpr_sp/figures/images/';
gbh = '/csail/fisher5/projects/Vseg_sp/saves/phone_gbh/';
swa = '/csail/fisher5/projects/Vseg_sp/saves/phone_swa/';
tsp = '/csail/fisher5/projects/Vseg_sp/saves/phone_tsp/';

W = 1000;

count = 1;
files = dir([tsp '*.png']);
for i=0:10:30
    im = imread([gbh files(i+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root 'phone_gbh_' num2str(count) '.png']);
    
    im = imread([swa files(i+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root 'phone_swa_' num2str(count) '.png']);
    
    im = imread([tsp files(i+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root 'phone_tsp_' num2str(count) '.png']);
    
    
    count = count + 1;
end


%%
name = 'dog'
root = '/afs/csail.mit.edu/u/j/jchang7/Desktop/Latex/jchang7_2013_cvpr_sp/figures/images/';
tsp = ['/csail/fisher5/projects/Vseg_sp/saves/' name '_tsp/'];
files = dir([tsp '*.png']);
for count=1:3
    im = imread([tsp files(20*(count-1)+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root name '_tsp_' num2str(count) '.png']);
end

name = 'car1_stabilized'
root = '/afs/csail.mit.edu/u/j/jchang7/Desktop/Latex/jchang7_2013_cvpr_sp/figures/images/';
tsp = ['/csail/fisher5/projects/Vseg_sp/saves/' name '_tsp/'];
files = dir([tsp '*.png']);
for count=1:3
    im = imread([tsp files(10*(count-1)+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root name '_tsp_' num2str(count) '.png']);
end

name = 'girl'
root = '/afs/csail.mit.edu/u/j/jchang7/Desktop/Latex/jchang7_2013_cvpr_sp/figures/images/';
tsp = ['/csail/fisher5/projects/Vseg_sp/saves/' name '_tsp/'];
files = dir([tsp '*.png']);
for count=1:3
    im = imread([tsp files(7*(count-1)+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root name '_tsp_' num2str(count) '.png']);
end

name = 'sample1'
root = '/afs/csail.mit.edu/u/j/jchang7/Desktop/Latex/jchang7_2013_cvpr_sp/figures/images/';
tsp = ['/csail/fisher5/projects/Vseg_sp/saves/' name '_tsp/'];
files = dir([tsp '*.png']);
for count=1:3
    im = imread([tsp files(4*(count-1)+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root name '_tsp_' num2str(count) '.png']);
end

%%
name = 'penguin'
root = '/afs/csail.mit.edu/u/j/jchang7/Desktop/Latex/jchang7_2013_cvpr_sp/figures/images/';
tsp = ['/csail/fisher5/projects/Vseg_sp/saves/' name '_tsp/'];
files = dir([tsp '*.png']);
for count=1:3
    im = imread([tsp files(13*(count-1)+1).name]);
    if (size(im,2)>W)
        im = imresize(im, [nan,W]);
    end
    imwrite(im, [root name '_tsp_' num2str(count) '.png']);
end