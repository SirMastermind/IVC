%clc;
close all;
clear;
beep on;

mode = 'picture'; % picture, movie
show = 'boxes'; % boxes, path, plot, speed, areas, positions

nFrame = 3065;
step = 5;

pathing = zeros(2, nFrame, 5);
speeds = zeros(1, nFrame, 5);
nF = 1;
nOD = 0;
merge = false;
split = false;
split_warning = false;
split_count = 0;

maxObjs = 5;

if(strcmp(mode,'picture'))
    path = 'DATASET1/TESTING/CAMERA1_JPEGS/';
    frameIdComp = 4;
    str = ['%s%.' num2str(frameIdComp) 'd.%s'];

    vid3D = zeros([576 768 nFrame/step]);
    d = zeros([576 768 nFrame/step]);
    for k = 1 : 1 : nFrame/step
        str1 = sprintf(str,path,k,'jpg');
        img = imread(str1);
        vid3D(:,:,k) = rgb2gray(img);
    end
else
    path = 'DATASET1/TRAINING/TRAINING/';
    file = 'camera1.mov';
    str1 = strcat(path, file);
    vid = VideoReader(str1);
    nFrame = 120*25;
    step = 5;

    vid3D = zeros([vid.Height vid.Width nFrame/step]);
    d = zeros([576 768 nFrame/step]);
    k=1;
    for i = 1:step:nFrame
        img = read(vid,i);
        vid3D(:,:,k)=rgb2gray(img);
        k = k+1;
    end
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
    [lb, num]= bwlabel(d(:, :, k));
    if(num > maxObjs)
        maxObjs = num;
    end
end

beep;

switch show
    case 'boxes'
        prev_num = 0;
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
                    split = false;
                    boundingBox = stats(i).BoundingBox;
                    if (abs(boundingBox(3)/boundingBox(4) - 1) < 0.09)
                        continue;
                    end
                    if (objects(1) > 5000 && num < prev_num)
                        merge = true;
                    elseif(merge && num > prev_num)
                        merge = false;
                        split = true;
                        split_count = 50;
                    end
                    if (boundingBox(3)/boundingBox(4) > 1) %boundingBox(3) = width; boundingBox(4) = height. When width > height, it is a car
                        color = 'r';
                    else
                        color = 'b';
                    end
                    rectangle('Position', boundingBox, 'EdgeColor',color, 'LineWidth', 2);
                    if(merge)
                        str = 'MERGE';
                        t = text(boundingBox(3), boundingBox(1),str);
                        s = t.Color;
                        t.Color = [1.0 0.0 0.0];
                        s = t.FontSize;
                        t.FontSize = 25;
                    end
                    if(not(split_warning) && (split || split_count > 0))
                        str = 'SPLIT';
                        t = text(boundingBox(2), boundingBox(1),str);
                        s = t.Color;
                        t.Color = [0.0 0.0 1.0];
                        s = t.FontSize;
                        t.FontSize = 25;
                        split_count = split_count - 1;
                        split = false;
                        split_warning = true;
                    end
                end
                split_warning = false;
            end
            drawnow;
            hold off;
            prev_num =  num;
        end
    case 'path'
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
    case 'speed'
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
                if(nF <= 5) %Nas primeiras 5 frames considera-se que nao ha speed
                    speeds(1, nF, l) = 0;
                else
                    speeds(1, nF, l) = 100 * (sqrt((pathing(1, nF, l) - pathing(1, nF-5, l))^2 + (pathing(2, nF, l) - pathing(2, nF-5, l))^2) / 5);
                end
%                 if (speeds(1, nF, l) > 10000)
%                     disp(strcat('Speed = ', num2str(speeds(1, nF, l))));
%                 end
            end

            imagesc(uint8(vid3D(:, :, k))); colormap gray; hold on;
            
            if num > 0
                for i = 1 : num
                    boundingBox = stats(i).BoundingBox;
                    speed = speeds(1, nF, i);
                    if (abs(boundingBox(3)/boundingBox(4) - 1) < 0.09)
                        continue;
                    end
                    if (speed > 10000) %It is a bike
                        rectangle('Position', boundingBox, 'EdgeColor', 'g', 'LineWidth', 2);
                    elseif (boundingBox(3)/boundingBox(4) > 1) %boundingBox(3) = width; boundingBox(4) = height. When width > height, it is a car
                        rectangle('Position', boundingBox, 'EdgeColor','r', 'LineWidth', 2);
                    else
                        rectangle('Position', boundingBox, 'EdgeColor','b', 'LineWidth', 2);
                    end
                end
            end
            drawnow;

            hold off;
            nF = nF+1;
        end
    case 'areas'
        numbers = zeros(size(d,3), maxObjs);
        for k = 1 : size(d,3)       
            % Find and label the different regions
            [lb, num]= bwlabel(d(:, :, k));
            % Get the stats of each label
            stats = regionprops(lb);

            % Compute area for each region
            objects = [stats.Area];
            for i = 1 : length(objects)
                numbers(k,i) = objects(i);
            end
        end
        figure; hold on;
        %plot(numbers);
        plot(numbers(:,1),'y-');
        plot(numbers(:,2),'r--');
        plot(numbers(:,3),'g:');
        plot(numbers(:,4),'b--o');
        plot(numbers(:,5),'k-*');
        legend('1','2','3','4','5');
        hold off;
    case 'positions'
        centroids = zeros(maxObjs, 2, size(d,3));
        for k = 1 : size(d,3)
            [lb, num]= bwlabel(d(:, :, k));
            stats = regionprops(lb);
            % Compute area for each region
            objects = [stats.Area];
            for i = 1 : size(lb,1) % For each lines
                for j = 1 : size(lb,2) % For each column
                    if lb(i,j) ~= 0 % If it's not background
                        centroids(lb(i,j),1,k) = centroids(lb(i,j),1,k) + i; % Sum the lines
                        centroids(lb(i,j),2,k) = centroids(lb(i,j),2,k) + j; % Sum the columns
                    end
                end
            end
            for i = 1 : length(objects) % For each object
                centroids(i,1,k) = centroids(i,1,k)/objects(i); % lines' = sum(lines)/area
                centroids(i,2,k) = centroids(i,2,k)/objects(i); % columns' = sum(columns)/area
            end
        end
        beep;
        subplot(2, 1, 1), hold on;
        for k = 1 : size(d,3)
            plot(k, centroids(1,1,k), 'b.',  'LineWidth', 2);
            plot(k, centroids(2,1,k), 'ro', 'LineWidth', 2);
            plot(k, centroids(3,1,k), 'g*', 'LineWidth', 2);
            plot(k, centroids(4,1,k), 'y.', 'LineWidth', 2);
            plot(k, centroids(5,1,k), 'ko', 'LineWidth', 2);
        end
        legend('1','2','3','4','5');
        hold off;
        beep;
        subplot(2, 1, 2), hold on;
        for k = 1 : size(d,3)
            plot(k, centroids(1,2,k), 'b.',  'LineWidth', 2);
            plot(k, centroids(2,2,k), 'ro', 'LineWidth', 2);
            plot(k, centroids(3,2,k), 'g*', 'LineWidth', 2);
            plot(k, centroids(4,2,k), 'y.', 'LineWidth', 2);
            plot(k, centroids(5,2,k), 'ko', 'LineWidth', 2);
        end
        legend('1','2','3','4','5');
        hold off;
        beep;
        figure, hold on;
        for k = 1 : size(d,3)
            plot(centroids(1,2,k), centroids(1,1,k), 'b.',  'LineWidth', 2);
            plot(centroids(2,2,k), centroids(2,1,k), 'ro', 'LineWidth', 2);
            plot(centroids(3,2,k), centroids(3,1,k), 'g*', 'LineWidth', 2);
            plot(centroids(4,2,k), centroids(4,1,k), 'y.', 'LineWidth', 2);
            plot(centroids(5,2,k), centroids(5,1,k), 'ko', 'LineWidth', 2);
        end
        legend('1','2','3','4','5');
        hold off;
end
beep;