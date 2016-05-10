%clc;
close all;
clear;

path = 'DATASET1/TESTING/CAMERA1_JPEGS/';
frameIdComp = 4;
str = ['%s%.' num2str(frameIdComp) 'd.%s'];

nFrame = 3065;
step = 5;
pathing = zeros(2, nFrame, 3);
nF = 0;
nOD = 0;

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
    nOD = max(nOD, length(objects));
    % Compute centroid for each region
    centroids = zeros(length(objects), 2); % To save (sum of lines, sum of columns) for each label
    for i = 1 : size(lb,1) % For each lines
        for j = 1 : size(lb,2) % For each column
            if lb(i,j) ~= 0 % If it's not background
                centroids(lb(i,j),1) = centroids(lb(i,j),1) + i; % Sum the lines
                centroids(lb(i,j),2) = centroids(lb(i,j),2) + j; % Sum the columns
            end
        end
    end

    for l = 1 : length(objects) % For each object
        centroids(l,1) = centroids(l,1)/objects(l); % lines' = sum(lines)/area
        centroids(l,2) = centroids(l,2)/objects(l); % columns' = sum(columns)/area
        pathing(1, nF, l) = centroids(l,1);
        pathing(2, nF, l) = centroids(l,2);
    end
    

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
    
    for i = 1 : num
        x_plot = [];
        y_plot = [];
        for j = 1 : nF
            x_plot = [ x_plot pathing(1, j, i) ];
            y_plot = [ y_plot pathing(2, j, i) ];
        end
        if(stats(i).BoundingBox(3) / stats(i).BoundingBox(4) > 1)
            plot(y_plot, x_plot, 'r.', 'MarkerSize', 5);
        else
            plot(y_plot, x_plot, 'b.', 'MarkerSize', 5);
        end
        drawnow;
    end
    hold off;
    nF = nF+1;
end