%clc;
close all;
clear;

path = 'DATASET1/TESTING/CAMERA1_JPEGS/';
frameIdComp = 4;
str = ['%s%.' num2str(frameIdComp) 'd.%s'];

nFrame = 3065;
step = 5;

vid3D = zeros([576 768 nFrame/step]);
d = zeros([576 768 nFrame/step]);
for k = 1 : 1 : nFrame/step
    str1 = sprintf(str,path,k,'jpg');
    img = imread(str1);
    vid3D(:,:,k) = rgb2gray(img);
end

bkg = median(vid3D,3);
%figure; imagesc(uint8(bkg)); colormap gray;

ths = 29;

bw = (abs(vid3D(:, :, 1) - bkg) > ths);
bw_final = bwareaopen(bw, 70);
se = strel('disk', 5);
bw_final = imclose(bw_final,se);
se = strel('disk', 8);
bw_final = imdilate(bw_final,se);
se = strel('disk', 3);
bw_final = imopen(bw_final,se);
bw_final = bwareaopen(bw_final, 150);
d(:, :, 1) = bw_final;
previous_bw = bw_final;

for k = 2 : (nFrame/step - 1)
    bw = (abs(vid3D(:, :, k) - bkg) > ths);
    bw_final = bwareaopen(bw, 70);
    se = strel('disk', 5);
    bw_final = imclose(bw_final,se);
    se = strel('disk', 8);
    bw_final = imdilate(bw_final,se);
    se = strel('disk', 3);
    bw_final = imopen(bw_final,se);
    bw_final = bwareaopen(bw_final, 150);
    d(:, :, k) = bw_final + previous_bw;
    previous_bw = bw_final;
end

for k = 1 : size(d,3)
    % Find and label the different regions
    [lb, num]= bwlabel(d(:, :, k));

    % Get the stats of each label
    stats = regionprops(lb);

    % Compute area for each region
    objects = [stats.Area];

    imagesc(uint8(vid3D(:, :, k))); colormap gray; hold on;
    if num > 0
        for i = 1 : num
            boundingBox = stats(i).BoundingBox;
            if (abs(boundingBox(3)/boundingBox(4) - 1) < 0.09)
                continue;
            end
            if (boundingBox(3)/boundingBox(4) > 1) %boundingBox(3) = width; boundingBox(4) = height. When width > height, it is a car
                rectangle('Position', boundingBox, 'EdgeColor','r', 'LineWidth', 2);
            else
                rectangle('Position', boundingBox, 'EdgeColor','b', 'LineWidth', 2);
            end
        end
    end
    drawnow;
    hold off;
end