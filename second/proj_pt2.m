close all;
clear all;

path = 'DATASET1/TRAINING/CAMERA1_JPEGS/';
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
figure; imagesc(uint8(bkg)); colormap gray;

%ths = graythresh(bkg)*255;
ths = 20;
figure;
for k = 1 : (nFrame/step - 1)
    bw = (abs(vid3D(:, :, k) - bkg) > 29);
    bw_final = bwareaopen(bw, 70);
    bw_final = bwmorph(bw_final,'hbreak'); % Removes H-connected pixels
    bw_final = bwmorph(bw_final,'spur'); % Removes spur pixels
    se = strel('disk', 5);
    bw_final = imclose(bw_final,se);
    bw_final = imdilate(bw_final,se);
    se = strel('disk', 3);
    bw_final = imopen(bw_final,se);
    bw_final = bwareaopen(bw_final, 150);
    bw_final = bwmorph(bw_final,'remove'); % Removes interior pixels
    bw_final = imfill(bw_final,'holes'); % Fills the regions
    d(:, :, k) = bw_final;
    
    % Find and label the different regions
    [lb, num]= bwlabel(bw_final);
    %figure, imshow(mat2gray(lb));

    % Get the stats of each label
    stats = regionprops(lb);

    % Compute area for each region
    objects = [stats.Area];

    if num > 0
        for i = 1 : num % For each sorted index of the objects
                    boundingBox = stats(i).BoundingBox; % Get the bounding box
                    %imshow(bw_final);
                    imagesc(uint8(vid3D(:, :, k))); colormap gray;
                    rectangle('Position', boundingBox, 'EdgeColor','r', 'LineWidth', 3);
                    drawnow, hold on;
        end
    end
    hold off;
end