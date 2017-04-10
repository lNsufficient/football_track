%%kmeans segmentation test
clear all
close all
load('testImg.mat');
load('mask.mat');

fieldmask = uint8(zeros(size(img)));
fieldmask(:,:,1) = mask;
fieldmask(:,:,2) = mask;
fieldmask(:,:,3) = mask;

mImg = fieldmask.*img;
% imagesc(img)
%% 
close all

figure
imagesc(mImg)
%img = double(img);
%figure
rmask = uint8(zeros(size(img))).*mImg;
rmask(:,:,1) = 1;
rmask = rmask.*fieldmask;

gmask = uint8(zeros(size(img))).*mImg;
gmask(:,:,2) = 1;
gmask = gmask.*fieldmask;

bmask = uint8(zeros(size(img)));
bmask(:,:,3) = 1;
bmask = bmask.*fieldmask;

r = img.*rmask;
g = img.*gmask;
b = img.*bmask;
figure
imagesc(r);
figure
imagesc(g);
figure
imagesc(b);
figure
imagesc(r+g+b)

r = r(:,:,1);
g = g(:,:,1);
b = b(:,:,1);
%% kmeans segmentation
nsegs = 20;
rgb = double([r(:) g(:) b(:)]);
[IDX, C] = kmeans(rgb, nsegs);
%% 
segments = reshape(IDX,[1028,5300]);
figure
imagesc(segments)

%% this did not work!












