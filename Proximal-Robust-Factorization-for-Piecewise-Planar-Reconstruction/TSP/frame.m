obj=VideoReader('C:\Users\IBM\Videos\Logitech Webcam\Video 3.mp4');

nFrames=obj.NumberOfFrames;
path=['.\sequences\Room\Frame'];

for i= 1 : nFrames
    imPath=char(strcat(path,int2str(i),'.bmp'));
    im=read(obj,i);
    imwrite(im,imPath,'bmp');
end