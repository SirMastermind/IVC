clear all;
close all;

% Load image
image = imread('Moedas1.jpg');
%imshow(image); title('Original');

%Convert image to gray
image_gray = rgb2gray(image);

% Use Otsu to get the threshold
bw = im2bw(image_gray, graythresh(image_gray));
%figure,imshow(bw);

% Use closure to get the shapes well defined
se = strel('disk', 15);
bw_final = imclose(bw,se);
bw_final = imopen(bw_final,se);
%figure,imshow(bw_final);

% Find and label the different regions
[lb, num]= bwlabel(bw_final);
%figure, imshow(mat2gray(lb));

% Get the stats of each label
stats = regionprops(lb);
    
% Compute area for each region
objects = [stats.Area];

string = sprintf('The number of objects in the image is %d.', length(objects));
disp(string);