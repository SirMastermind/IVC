clear all;
close all;

image = imread('images/MoedasA.jpg');

image_gray = rgb2gray(image);

bw = im2bw(image_gray, graythresh(image_gray));

% Use closure to get the shapes well defined
se = strel('disk', 5);
bw_final = imclose(bw,se);
bw_final = bwmorph(bw_final,'hbreak'); % Removes H-connected pixels
bw_final = bwmorph(bw_final,'spur'); % Removes spur pixels
bw_final = bwareaopen(bw_final, 20); % Remove all object containing fewer than 20 pixels
bw_final = bwmorph(bw_final,'bridge'); % Removes bridge pixels
bw_final = bwmorph(bw_final,'remove'); % Removes interior pixels
bw_final = imfill(bw_final,'holes'); % Fills the regions

[centers,radii] = imfindcircles(bw_final, [50 150], 'Sensitivity', 0.9);
figure('Name','Circles found','NumberTitle','off'), imshow(image_gray);
hold on;
viscircles(centers, radii,'EdgeColor','b');

[X,Y] = meshgrid(1:size(bw_final,2), 1:size(bw_final,1));
IDs = zeros(size(bw_final));
for idx = 1 : numel(radii)
    r = radii(idx);
    cen = centers(idx,:);

    loc = (X - cen(1)).^2 + (Y - cen(2)).^2 <= r^2;
    IDs(loc) = idx;
end

figure, imshow(IDs, []);

figure, imshow(IDs > 0);

% Find and label the different regions
[lb, num]= bwlabel(bw_final);

% Get the stats of each label
stats = regionprops(lb);
    
% Compute area for each region
objects = [stats.Area];