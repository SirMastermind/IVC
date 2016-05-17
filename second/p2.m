%clc;
close all;
clear;
beep on;

mode = 'movie'; % picture, movie

nFrame = 3065;
step = 5;

nF = 1;
merge = false;
split = false;
split_warning = false;
merge_warning = false;
split_count = 0;

path = 'DATASET1/TRAINING/TRAINING/';
file = 'camera1.mp4';
str1 = strcat(path, file);
vid = VideoReader(str1);
nFrame = 120*25;
step = 5;
prev_num = 0;

vid3D = zeros([vid.Height vid.Width nFrame/step]);
d = zeros([576 768 nFrame/step]);
previous_bw = zeros([576 768]);
bkg = zeros(vid.Height, vid.Width);

k = 1;
alfa = 0.05;
ths = 29;

for i = 1 : step : nFrame
    img = read(vid,i);
    vid3D(:,:,k) = rgb2gray(img);
    if i < 300
        bkg = alfa * double(vid3D(:,:,k)) + (1-alfa) * double(bkg);
    else
        bw = (abs(vid3D(:,:,k) - bkg) > ths);
        bw_final = bwareaopen(bw, 70);
        %bw_final = bwmorph(bw_final, 'close');
        se = strel('disk', 5);
        bw_final = imclose(bw_final,se);
        se = strel('disk', 7);
        bw_final = imdilate(bw_final,se);
        se = strel('disk', 5);
        bw_final = imopen(bw_final,se);
        bw_final = bwareaopen(bw_final, 150);
        bw_image = (bw_final + previous_bw) > 0;
        previous_bw = bw_final;
       
        [lb, num]= bwlabel(bw_image);
        stats = regionprops(lb);
        objects = [stats.Area];

        imagesc(uint8(vid3D(:, :, k))); colormap gray; hold on;
        axis off;

        if num > 0
            for j = 1 : num
                split = false;
                boundingBox = stats(j).BoundingBox;
                if (find(objects(:) > 4500) & num < prev_num)
                    if merge == false
                        disp('There was a merge.');
                    end
                    merge = true;
                elseif(merge && num > prev_num)
                    merge = false;
                    if split == false
                        disp('There was a split.');
                    end
                    split = true;
                    split_count = 15;
                end
                if (abs(boundingBox(3)/boundingBox(4) - 1) < 0.08)
                    continue;
                end
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
    end
    k = k + 1;
end

beep;