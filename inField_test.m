%% inField test

corners = [240 817;
           1910 78;
           3525 90;
           5130 818];

v = VideoReader('stitched.mkv');
v.currentTime = 13;
img1 = rgb2gray(readFrame(v));
%%
D = size(img1);
coords = zeros;
mask = uint8(zeros(D));
for y = 1:D(1)
    for x = 1:D(2)
        X = [x,y];
        mask(y,x) = inField(corners,X);
    end
end

%%
spy(mask)
figure
maskedImg = uint8(mask).*(img1);
imagesc(maskedImg)
save('mask','mask')
hold on
plot(corners(:,1),corners(:,2),'o');