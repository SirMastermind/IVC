% Clear workspace
clc;
close all;
clear;

% Enable sound
beep on;

% Show bar
prompt = {'Enter source type (video or picture):','Enter mode (box, path or graph):'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'video','box'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

type = answer(1);
mode = answer(2);

% Variables
k = 1;
alfa = 0.05;
ths = 29;
n_train = 300;

nFrame = 3065;
step = 5;

nF = 1;
nOD = 0;
merge = false;
split = false;
split_warning = false;
merge_warning = false;
split_count = 0;

nFrame = 3000;
step = 1;
prev_num = 0;

pathing = zeros(2, nFrame, 5);
d = zeros([576 768 nFrame/step]);
previous_bw = zeros([576 768]);

if(strcmp(type,'picture')) % If the user wants picture type
    path = 'DATASET1/TESTING/CAMERA1_JPEGS/';
    frameIdComp = 4;
    str2 = ['%s%.' num2str(frameIdComp) 'd.%s'];
    vid3D = zeros([576 768 nFrame/step]);
    bkg = zeros(576, 768);
elseif(strcmp(type,'video')) % If the user wants movie type
    path = 'DATASET1/TRAINING/TRAINING/';
    file = 'camera1.mp4';
    str1 = strcat(path, file);
    vid = VideoReader(str1);
    vid3D = zeros([vid.Height vid.Width nFrame/step]);
    bkg = zeros(vid.Height, vid.Width);
else
    errordlg('Invalid type','Type Error');
end


h = waitbar(0, 'Getting the background, please wait...');

if(strcmp(type,'picture')) % If the user wants picture type
   for i = 1 : step : n_train
        str1 = sprintf(str2,path,i,'jpg');
        img = imread(str1);
        vid3D(:,:,k) = rgb2gray(img);
        bkg = alfa * double(vid3D(:,:,k)) + (1-alfa) * double(bkg);
        waitbar(i/n_train, h); 
        k = k + 1;
    end
else
    for i = 1 : step : n_train
        img = read(vid,i);
        vid3D(:,:,k) = rgb2gray(img);
        bkg = alfa * double(vid3D(:,:,k)) + (1-alfa) * double(bkg);
        waitbar(i/n_train, h); 
        k = k + 1;
    end
end

close(h);
beep;
if(strcmp(mode,'box'))
    figure('Name','Applying algorithm with box mode','NumberTitle','off'), hold on;
    for i = n_train : step : nFrame
        if(strcmp(type,'picture'))
            str1 = sprintf(str2,path,i,'jpg');
            img = imread(str1);
        else
            img = read(vid,i);
        end

        vid3D(:,:,k) = rgb2gray(img);

        bw = (abs(vid3D(:,:,k) - bkg) > ths);
        bw_final = bwareaopen(bw, 85);
        se = strel('disk', 5);
        bw_final = imclose(bw_final,se);
        se = strel('disk', 5);
        bw_final = imdilate(bw_final,se);
        se = strel('disk', 5);
        bw_final = imclose(bw_final,se);
        bw_final = bwareaopen(bw_final, 150);
        bw_image = (bw_final + previous_bw) > 0;
        previous_bw = bw_final;

        [lb, num]= bwlabel(bw_image);
        stats = regionprops(lb);
        objects = [stats.Area];

        imshow(img); hold on;

        if num > 0
            for j = 1 : num
                split = false;
                boundingBox = stats(j).BoundingBox;
                if (find(objects(:) > 4500) & num < prev_num)
                    merge = true;
                elseif(merge && num > prev_num)
                    merge = false;
                    split = true;
                    split_count = 15;
                end
%                 if (abs(boundingBox(3)/boundingBox(4) - 1) < 0.08)
%                     continue;
%                 end
                if (boundingBox(3)/boundingBox(4) > 1) %boundingBox(3) = width; boundingBox(4) = height. When width > height, it is a car
                    t = text(boundingBox(1), boundingBox(2) - 12, 'Car');
                    t.Color = [1.0 0.0 0.0];
                    t.FontSize = 16;
                    color = 'r';
                elseif (abs(boundingBox(3)/boundingBox(4) - 1) < 0.2)
                    t = text(boundingBox(1), boundingBox(2) - 12, 'Other');
                    t.Color = [0.0 1.0 0.0];
                    t.FontSize = 16;
                    color = 'g';
                else
                    t = text(boundingBox(1), boundingBox(2) - 12, 'Person');
                    t.Color = [0.0 0.0 1.0];
                    t.FontSize = 16;
                    color = 'b';
                end
                rectangle('Position', boundingBox, 'EdgeColor',color, 'LineWidth', 2);
                if(not(merge_warning) && merge)
                    str = 'MERGE';
                    t = text(0, 500,str);
                    s = t.Color;
                    t.Color = [1.0 1.0 0.0];
                    s = t.FontSize;
                    t.FontSize = 25;
                    merge_warning = true;
                end
                if(not(split_warning) && (split || split_count > 0))
                    str = 'SPLIT';
                    t = text(0, 420,str);
                    t.Color = [0.0 1.0 1.0];
                    t.FontSize = 25;
                    split_count = split_count - 1;
                    split = false;
                    split_warning = true;
                end
            end
            split_warning = false;
            merge_warning = false;
        end
        drawnow;
        hold off;
        prev_num =  num;
        k = k + 1;
    end
elseif(strcmp(mode,'path'))
    figure('Name','Applying algorithm with path mode','NumberTitle','off');
    %figure;
    for i = n_train : step : nFrame
        if(strcmp(type,'picture'))
            str1 = sprintf(str,path,i,'jpg');
            img = imread(str1);
        else
            img = read(vid,i);
        end

        vid3D(:,:,k) = rgb2gray(img);

        bw = (abs(vid3D(:,:,k) - bkg) > ths);
        bw_final = bwareaopen(bw, 85);
        se = strel('disk', 5);
        bw_final = imclose(bw_final,se);
        se = strel('disk', 5);
        bw_final = imdilate(bw_final,se);
        se = strel('disk', 5);
        bw_final = imclose(bw_final,se);
        bw_final = bwareaopen(bw_final, 150);
        bw_image = (bw_final + previous_bw) > 0;
        previous_bw = bw_final; 
        
        [lb, num]= bwlabel(bw_image);
        stats = regionprops(lb);
        objects = [stats.Area];
        
        nOD = max(nOD, length(objects));
        
        centroids = zeros(length(objects), 2); % To save (sum of lines, sum of columns) for each label
        for a = 1 : size(lb,1) % For each lines
            for j = 1 : size(lb,2) % For each column
                if lb(a,j) ~= 0 % If it's not background
                    centroids(lb(a,j),1) = centroids(lb(a,j),1) + a; % Sum the lines
                    centroids(lb(a,j),2) = centroids(lb(a,j),2) + j; % Sum the columns
                end
            end
        end

        for l = 1 : length(objects) % For each object
            centroids(l,1) = centroids(l,1)/objects(l); % lines' = sum(lines)/area
            centroids(l,2) = centroids(l,2)/objects(l); % columns' = sum(columns)/area
            pathing(1, nF, l) = centroids(l,1);
            pathing(2, nF, l) = centroids(l,2);
        end

        %imshow(img); hold on;
        imagesc(uint8(vid3D(:, :, k))); colormap gray; hold on;
        axis off;

        for a = 1 : num
            x_plot = [];
            y_plot = [];
            for j = 1 : nF
                x_plot = [ x_plot pathing(1, j, a) ];
                y_plot = [ y_plot pathing(2, j, a) ];
            end
            if(stats(a).BoundingBox(3) / stats(a).BoundingBox(4) > 1)
                plot(y_plot, x_plot, 'r.', 'MarkerSize', 5);
            else
                plot(y_plot, x_plot, 'b.', 'MarkerSize', 5);
            end
            drawnow;
        end
        hold off;
        nF = nF + 1;
        k = k + 1;
    end
end
beep;