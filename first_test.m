new_run = 0;

if new_run
    clear;
    v = VideoReader('stitched.m4v');
end

v.currentTime = 13; 

do_plot = 0;
if do_plot
    while hasFrame(v)
        figure(1)
        im = readFrame(v);
        %imagesc(im);
        %figure(2)
        imbw = rgb2gray(im);
        imagesc(imbw);
        colormap('gray')
        pause;
        title(v.CurrentTime)
    end
end
%% Find waldo approach
ball_xy_time = [3128 142 13.04; 3041 152 13.84]
makebw = 1

v.CurrentTime = ball_xy_time(1, 3);
im = readFrame(v);
if makebw
    im = rgb2gray(im);
end
im = double(im);
ball = [ball_xy_time(i,1) - ballsize, ball_xy_time(i,2)-ballsize, ballsize*2, ballsize*2];
small_ball = imcrop(im, ball);

figure(2)
imagesc(small_ball);

patch = fliplr(flipud(small_ball));

figure(5)
imagesc(patch)
colormap('gray')

%%

for i = 2:size(ball_xy_time,1)
    v.CurrentTime = ball_xy_time(i, 3)
    im = readFrame(v);
    if makebw
        im = rgb2gray(im);
    end
    im = double(im)
%    for ballsize = 3:7

        

        
        figure(3)
        imagesc(im);
        v.CurrentTime
        colormap('gray')
        
        AA = sum(sum(patch.*patch));
        if makebw == 0
            AA = sum(AA)
        end
        BBt = im.*im;
        if makebw
            BB = conv2(BBt,ones(size(patch)),'valid');
            AB = conv2(im,patch,'valid');
        else
            BB = convn(BBt,ones(size(patch)),'valid');
            AB = convn(im,patch,'valid');
        end
        res = AA - 2*AB + BB;
        [vals, m] = min(res);
        [~, n] = min(vals);
        m = m(n);
        [m_im, n_im] = size(im);
        [m_conv, n_conv] = size(AB);
        edge_m = (m_im - m_conv)/2;
        edge_n = (n_im - n_conv)/2;
        found_m_n = [m+edge_m n+edge_n]
        
        
        
        figure(2)
        imagesc(small_ball);
 
        figure(4)
        str = sprintf('m: %d n: %d', m, n);
        imagesc(res)
        title(str)
        hold on;
        
        antiRes = (max(max(res)) - res).^7;
        
        figure(6)
        clf;
        imagesc(antiRes)
        title(str)
        colormap('gray')
        
        
        plot(n, m, '*r')
        if makebw
            colormap('gray')
        end
        pause;

        
%    end
    
end


%%
corners_im = [1909 78; 3526 90; 5131 819; 238 819];
swedbank_stadion_width = 105; %meter
swedbank_stadion_height = 68; %meter
corners_malmo = [0 0; 
    swedbank_stadion_width 0;
    swedbank_stadion_width swedbank_stadion_height;
    0 swedbank_stadion_height];

corners_malmo = corners_malmo*25; %meter to pixel - never set to anything larger than 40!!!

x1 = [corners_im'; [1 1 1 1]];
x2 = [corners_malmo'; [1 1 1 1]];

H = getH(x1, x2, getN(x1), getN(x2))
figure(2)
plot(corners_malmo(:,1), corners_malmo(:,2), '*r')
hold on;
hx1 = pflat(H*x1); 
plot(hx1(1,:), hx1(2,:), 'ob')

%tform = maketform('projective',H');
tform = projective2d(H');
topdown = imwarp(im,tform);
figure(3)
imagesc(topdown)
%%
