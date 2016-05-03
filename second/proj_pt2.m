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
    bw_final = bwareaopen(bw, 75); % Remove all object containing fewer than 20 pixels
    bw_final = bwmorph(bw_final,'hbreak'); % Removes H-connected pixels
    bw_final = bwmorph(bw_final,'spur'); % Removes spur pixels
    se = strel('disk', 5);
    bw_final = imclose(bw_final,se);
    bw_final = imdilate(bw_final,se);
    bw_final = bwmorph(bw_final,'remove'); % Removes interior pixels
    bw_final = imfill(bw_final,'holes'); % Fills the regions
    d(:, :, k) = bw_final;
    imshow(bw_final); drawnow
    hold off;
end