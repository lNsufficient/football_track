%% tester script
close all
v = VideoReader('stitched.m4v');
v.currentTime = 13;

corners = [240 817;
           1910 78;
           3525 90;
           5130 818];

img1 = readFrame(v);
img2 = readFrame(v);
img3 = readFrame(v);

diff1 = imcomplement(abs(img1-img2));
diff2 = imcomplement(abs(img2-img3));
diff3 = imcomplement(abs(img1-img3));
%%
imagesc(img1)
hold on
plot(corners(:,1),corners(:,2))