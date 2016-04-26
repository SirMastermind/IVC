clear all;
close all;
warning off;

divider = sprintf('<----------------------------------------------->');

disp('Welcome to the first project!');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while(true)
    disp(' ');
    disp('Please, choose a picture:');
    disp('     1 - Moedas1.jpg;');
    disp('     2 - Moedas2.jpg;');
    disp('     3 - Moedas3.jpg;');
    disp('     4 - Moedas4.jpg;');
    disp('     5 - Moedas5.jpg;');
    disp('     6 - Moedas6.jpg;');
    disp('     7 - Moedas7.jpg;');
    disp('     8 - Moedas8.jpg;');
    disp('     9 - Moedas9.jpg;');
    disp('     10 - cube.png;');
    disp('     11 - ring.png;');
    disp('     0 - Exit Matlab.');
    disp(' ');
    option = input('Your choice: ');
    disp(' ');
    switch option
        case 1
            image = imread('images/Moedas1.jpg');
            break;
        case 2
            image = imread('images/Moedas2.jpg');
            break;
        case 3
            image = imread('images/Moedas3.jpg');
            break;
        case 4
            image = imread('images/Moedas4.jpg');
            break;
        case 5
            image = imread('images/Moedas5.jpg');
            break;
        case 6
            image = imread('images/Moedas6.jpg');
            break;
        case 7
            image = imread('images/Moedas7.jpg');
            break;
        case 8
            image = imread('images/Moedas8.jpg');
            break;
        case 9
            image = imread('images/Moedas9.jpg');
            break;
        case 10
            image = imread('images/cube.png');
            break;
        case 11
            image = imread('images/ring.png');
            break;
        case 0
            exit;
        otherwise
            disp('Please, insert a valid image.');
    disp(divider);
    end
end

disp(' ');
disp('The program started processing your image...');
disp(' ');

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert image to gray
image_gray = rgb2gray(image);

% Use Otsu to get the threshold
bw = im2bw(image_gray, graythresh(image_gray));
%figure,imshow(bw);

% Use closure to get the shapes well defined
se = strel('disk', 3);
bw_final = imclose(bw,se);
bw_final = bwmorph(bw_final,'hbreak'); % Removes H-connected pixels.
bw_final = bwmorph(bw_final,'spur'); % Removes spur pixels.
bw_final = bwmorph(bw_final,'clean'); % Removes isolated pixels (individual 1s that are surrounded by 0s).
bw_final = bwmorph(bw_final,'remove'); % Removes interior pixels.
bw_final = imfill(bw_final,'holes'); % Fills the regions.
%figure,imshow(bw_final);

% Find and label the different regions
[lb, num]= bwlabel(bw_final);
%figure, imshow(mat2gray(lb));

% Get the stats of each label
stats = regionprops(lb);
    
% Compute area for each region
objects = [stats.Area];
   
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

for k = 1 : length(objects) % For each object
    centroids(k,1) = centroids(k,1)/objects(k); % lines' = sum(lines)/area
    centroids(k,2) = centroids(k,2)/objects(k); % columns' = sum(columns)/area
end

bw_centroids = bw_final;
for i = 1 : size(centroids,1)
    bw_centroids(uint16(centroids(i,1)), uint16(centroids(i,2))) = 0; % Mark the centroid with a black pixel
end

% Compute perimeter for each region
perimeters8 = bwperim(bw_final, 8); % bwperim(binary image, connectivity between neighbors)

% Compute distance for each object to other
distance = zeros(size(centroids,1), size(centroids,1));
for i = 1 : size(distance,1) % For each lines
    for j = 1 : size(distance,2) % For each column
        distance(i,j) = sqrt((centroids(j,1) - centroids(i,1))^2 + (centroids(j,2) - centroids(i,2))^2); % Vector's length
        distance(j,i) = sqrt((centroids(j,1) - centroids(i,1))^2 + (centroids(j,2) - centroids(i,2))^2); % Mirror
    end
end
    
% Compute differences between objects to find similarity between regions
perimeters4 = bwperim(bw_final, 4); % Find perimeters for each region with 4 neighbors
labeled_perimeters8 = bwlabel(perimeters8 - perimeters4); % Find perimeters for each region with only oblique neighbors
labeled_perimeters4 = bwlabel(perimeters4); % Label the regions of the perimters
individual_circularities2 = zeros(1, length(objects)); % To save the computations
individual_perimeters = zeros(1,length(objects)); % Vector to save each object's perimeter
for i = 1 : length(individual_perimeters)
    count8s = 0; % Counter perimeter8 - perimeter4
    count4s = 0; % Counter perimeter4
    radius_mean = 0; %Counter for radius_mean
    radius_var = 0; %Counter for radius_var
    for j = 1 : size(labeled_perimeters8,1)
        for k = 1 : size(labeled_perimeters8,2)
            radius_mean = radius_mean + norm([j-centroids(i,1) k-centroids(i,2)]); %Sum part of the radius mean calculus
            if labeled_perimeters8(j,k) == i
                count8s = count8s + 1;
            end
            if labeled_perimeters4(j,k) == i
                count4s = count4s + 1;
            end
        end
    end
    
    radius_mean = (1/(j*k)) * radius_mean; % Final calculus of radius mean
    
    % New iteration to calculate radius variation
    for j = 1 : size(labeled_perimeters8,1)
        for k = 1 : size(labeled_perimeters8,2)
            radius_var = radius_var + (norm([j-centroids(i,1) k-centroids(i,2)]) - radius_mean)^2; % Sum part of the radius variation calculus
        end
    end
    
    radius_var = sqrt((1/(j*k)) * radius_var); % Final calculus of radius variation
    individual_circularities2(i) = radius_mean / radius_var; % Compute individual circularity for each region
    individual_perimeters(i) = count4s + sqrt(2) * count8s; % Calculate object-i's perimeter
end

individual_circularities = zeros(1,length(objects));
for i = 1 : length(individual_circularities)
    individual_circularities(i) = (individual_perimeters(i)^2) / objects(i); % Calculate object-i's circularity
end

% Compute quadtree
output_size = [power(2,nextpow2(size(image_gray, 1))), power(2,nextpow2(size(image_gray, 2)))]; % Output size must be power of 2
S = qtdecomp(imresize(image_gray, output_size) , 0.27); % Resize it to the output size and decomposition in quadtree
blocks = repmat(uint8(0),size(S)); % Get the blocks

for dim = [512 256 128 64 32 16 8 4 2 1]; % For each dimension of power 2
  numblocks = length(find(S==dim)); % Find blocks the size of the current dimension 
  if (numblocks > 0)        
    values = repmat(uint8(1),[dim dim numblocks]); % New block with the dimensions
    values(2:dim,2:dim,:) = 0; % Put it to zeros in the middle
    blocks = qtsetblk(blocks,S,dim,values); % Set the blocks with the new values
  end
end

% Put the blocks' perimeters to 1
blocks(end,1:end) = 1;
blocks(1:end,end) = 1;

blocks = imresize(blocks, [size(image_gray, 1), size(image_gray, 2)]); % Resize for the previous size

% Compute derivative in y
Hy = fspecial('sobel');
imgdy = filter2(Hy,image(:,:,1));

% Compute derivative in x
Hx = fliplr(Hy');
imgdx = filter2(Hx,image(:,:,1));

% Compute module
modG = sqrt(imgdy.^2 + imgdx.^2);

% Hough Transform
hough_image = edge(image_gray, 'canny');
[H,T,R] = hough(hough_image,'RhoResolution',0.5,'ThetaResolution',0.5);

% Limited theta range Hough transform
[H_h,T_h,R_h] = hough(image_gray, 'Theta', 44:0.5:46);

% Gradient magnitude and gradient direction
[Gmag, Gdir] = imgradient(image_gray,'sobel');

time = toc; % To measure the time it took to process the image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Program's flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Image processed!');
disp(sprintf('It took %.2f seconds.', time));
disp(' ');
disp('This is the main menu.');

while(true)
    % User interface
    disp(' ');
    disp('Please, input your option:');
    disp('     1 - Show objects;');
    disp('     2 - Count objects;');
    disp('     3 - Show centroids;');
    disp('     4 - Show perimeters (processed binary image);');
    disp('     5 - Show perimeters (original image);');
    disp('     6 - Show areas by grayscale objects;');
    disp('     7 - Show areas by labels;');
    disp('     8 - Show distance between objects;');
    disp('     9 - Show similar objects;');
    disp('     10 - Show edges;');
    disp('     11 - Show quadtree;');
    disp('     12 - Show derivative in x;');
    disp('     13 - Show derivative in y;');
    disp('     14 - Show module;');
    disp('     15 - Show boundaries (unprocessed binary image);');
    disp('     16 - Show boundaries (processed binary image);');
    disp('     17 - Show histogram RGB;');
    disp('     18 - Show histogram grayscale;');
    disp('     19 - See specific object from mask;');
    disp('     20 - Hough transform;');
    disp('     21 - Limited theta range Hough transform;');
    disp('     22 - Show gradient magnitude and direction;');
    disp('     23 - Show skeleton of the unprocessed binary image;');
    disp('     24 - Show skeleton of the processed binary image;');
    disp('     0 - Exit.');
    disp(' ');
    option = input('Your option: ');
    disp(' ');
    
    % Program's flow
    switch option
        case 1
            close all;
            figure('Name','Objects detected','NumberTitle','off'), hold on;
            imshow(image_gray.*uint8(bw_final));
            disp(divider);
        case 2
            close all;
            string = sprintf('The number of objects in the image is %d.', length(objects));
            disp(string);
            disp(divider);
        case 3
            close all;
            for i = 1 : size(centroids,1)
                string = sprintf('Object %d has centroid in (%0.2f, %0.2f).', i, centroids(i,1), centroids(i,2));
                disp(string);
            end
            figure('Name','Centroids of the objects','NumberTitle','off'), imshow(image_gray);
            hold on;
            x_plot = [];
            y_plot = [];
            for i = 1 : length(centroids)
                x_plot = [ x_plot centroids(i, 1) ];
                y_plot = [ y_plot centroids(i, 2) ];
            end
            plot(y_plot, x_plot, 'r.', 'MarkerSize', 30);
            legend('Centroid');
            disp(divider);
        case 4
            close all;
            figure('Name','Perimeters of the objects in binary image','NumberTitle','off'), imshow(perimeters8); 
            disp(divider);
        case 5
            close all;
            figure('Name','Perimeters of the objects in original image','NumberTitle','off'), imshow(image);
            hold on
            x_perimeters = [];
            y_perimeters = [];
            for i = 1:size(perimeters8,1)
                for j = 1:size(perimeters8,2)
                    if perimeters8(i,j) == 1
                        x_perimeters = [x_perimeters i];
                        y_perimeters = [y_perimeters j];
                    end
                end
            end
            plot(y_perimeters, x_perimeters, 'r.', 'LineWidth', 5);
            set(gca,'Color','None');
            legend('Perimeter');
            hold off;
            disp(divider);
        case 6
            close all;
            figure('Name','Area in grayscale of the objects','NumberTitle','off'), imshow(image_gray.*uint8(bw_final));
            for i=1 : size(centroids, 1)
                x1 = centroids(i,2);
                y1 = centroids(i,1);
                str = num2str(objects(i));
                t = text(x1,y1,str);
                s = t.Color;
                t.Color = [1.0 0.0 0.0];
                s = t.FontSize;
                t.FontSize = 30;
            end
            figure('Name','Increasing area from left to right','NumberTitle','off'), hold on;
            [sorted, indexes] = sort(objects); % Sort areas
            for i = 1 : length(indexes) % For each sorted index of the objects
                boundingBox = stats(indexes(i)).BoundingBox; % Get the bounding box
                part = imcrop(image, boundingBox); % Crop the original image to keep only the original bounding box
                subplot(1,length(indexes),i), imshow(part); % Put the small picture in one full figure
                title_part = sprintf('%d', i); % Name it
                title(title_part);
            end
            disp(divider);
        case 7
            close all;
            figure('Name','Areas with different colors by labels','NumberTitle','off'), imshow(mat2gray(lb));
            for i=1 : size(centroids, 1)
                x1 = centroids(i,2);
                y1 = centroids(i,1);
                str = num2str(objects(i));
                t = text(x1,y1,str);
                s = t.Color;
                t.Color = [1.0 0.0 0.0];
                s = t.FontSize;
                t.FontSize = 30;
            end
            disp(divider);
        case 8
            close all;
            disp('Click on two objects.');
            figure('Name','Click on two objects','NumberTitle','off'), imshow(image), hold on;
            N = 0;
            but = 1;
            while (but == 1)
                [ci,li,but] = ginput(1);
                if but == 1 % Left click
                    N     = N+1;
                    cp(N) =  ci;
                    lp(N) =  li;
                    plot(ci,li,'r.','MarkerSize',8); drawnow;
                end
                if size(cp,2) == 2
                    break;
                end
            end
            close;
            disp('Please, wait a second.');
            figure('Name','Distance between centroids of two objects','NumberTitle','off'), imshow(image), hold on;
            object_1 = lb(uint16(lp(1)),uint16(cp(1)));
            object_2 = lb(uint16(lp(2)),uint16(cp(2)));
            centroid_1 = centroids(object_1, :);
            centroid_2 = centroids(object_2, :);
            plot(centroid_1(2), centroid_1(1), 'r.-', 'MarkerSize', 8); drawnow;
            plot(centroid_2(2), centroid_2(1), 'r.-', 'MarkerSize', 8); drawnow;
            plot([centroid_1(2) centroid_2(2)], [centroid_1(1) centroid_2(1)], 'r.-', 'MarkerSize', 5); drawnow;
            distance_centroids = sqrt((centroid_1(1) - centroid_2(1))^2 + (centroid_1(2) - centroid_2(2))^2);
            hold on;
            x1 = centroid_1(1);
            y1 = centroid_1(2);
            str = ['Distance is ', num2str(uint16(distance_centroids))];
            t = text(x1,y1,str);
            s = t.Color;
            t.Color = [1.0 0.0 0.0];
            s = t.FontSize;
            t.FontSize = 25;

            legend('Distance between selected dots');
            
            disp('Right click to close the image.');
            but = 32;
            while (but == 32)
                [ci,li,but] = ginput(1);
                if but == 32 % Right click
                    break;
                end
            end
            close; % Closes picture
            clear cp lp;
            disp(divider);
        case 9
            close all;
            while(true)
                disp(' ');
                disp('Please, choose a method:');
                disp('     1 - Circularity 1;');
                disp('     2 - Circularity 2;');
                disp('     0 - Go to previous menu.');
                disp(' ');
                option = input('Your method: ');
                switch option
                    case 1
                        close all;
                        disp('Click on an object to see which are the most similiar to it');
                        figure('Name','Similarity using circularity 1','NumberTitle','off'), imshow(image), hold on;
                        N = 1;
                        but = 1;
                        [ci,li,but] = ginput(1);
                        if but == 1 % Left click
                            cp(N) =  ci;
                            lp(N) =  li;
                        end

                        object_1 = lb(uint16(lp(1)),uint16(cp(1)));
                        circularity = individual_circularities(object_1);

                        %Create a vector with the module diferences between the
                        %circularity of the chosen object and all the other objects
                        difs = zeros(1, length(individual_circularities)); 
                        for i=1 : length(difs)
                            difs(i) = abs(circularity - individual_circularities(i));
                        end

                        %Sort the labels (indexes), from the most similiar to the less,
                        %and prints in the center of each object
                        [sorted, indexes] = sort(difs);
                        for i=1 : length(indexes)
                            x1 = centroids(indexes(i),2);
                            y1 = centroids(indexes(i),1);
                            str = [num2str(i)];
                            t = text(x1,y1,str);
                            s = t.Color;
                            t.Color = [1.0 0.0 0.0];
                            s = t.FontSize;
                            t.FontSize = 35;
                        end
                        figure('Name','Similarity decreasing from left to right','NumberTitle','off'), hold on;
                        for i = 1 : length(indexes) % For each sorted index of the objects
                            boundingBox = stats(indexes(i)).BoundingBox; % Get the bounding box
                            part = imcrop(image, boundingBox); % Crop the original image to keep only the original bounding box
                            subplot(1,length(indexes),i), imshow(part); % Put the small picture in one full figure
                            title_part = sprintf('%d', i); % Name it
                            title(title_part);
                        end
                    case 2
                        close all;
                        disp('Click on an object to see which are the most similiar to it');
                        figure('Name','Similarity using circularity 2','NumberTitle','off'), imshow(image), hold on;
                        N = 1;
                        but = 1;
                        [ci,li,but] = ginput(1);
                        if but == 1 % Left click
                            cp(N) =  ci;
                            lp(N) =  li;
                        end

                        object_1 = lb(uint16(lp(1)),uint16(cp(1)));
                        circularity = individual_circularities2(object_1);

                        %Create a vector with the module diferences between the
                        %circularity of the chosen object and all the other objects
                        difs = zeros(1, length(individual_circularities2)); 
                        for i=1 : length(difs)
                            difs(i) = abs(circularity - individual_circularities2(i));
                        end

                        %Sort the labels (indexes), from the most similiar to the less,
                        %and prints in the center of each object
                        [sorted, indexes] = sort(difs);
                        for i=1 : length(indexes)
                            x1 = centroids(indexes(i),2);
                            y1 = centroids(indexes(i),1);
                            str = [num2str(i)];
                            t = text(x1,y1,str);
                            s = t.Color;
                            t.Color = [1.0 0.0 0.0];
                            s = t.FontSize;
                            t.FontSize = 35;
                        end
                        figure('Name','Similarity decreasing from left to right','NumberTitle','off'), hold on;
                        for i = 1 : length(indexes) % For each sorted index of the objects
                            boundingBox = stats(indexes(i)).BoundingBox; % Get the bounding box
                            part = imcrop(image, boundingBox); % Crop the original image to keep only the original bounding box
                            subplot(1,length(indexes),i), imshow(part); % Put the small picture in one full figure
                            title_part = sprintf('%d', i); % Name it
                            title(title_part);
                        end
                    case 0
                        close all;
                        break;
                    otherwise
                        close all;
                        disp('Please, insert a valid method.');
                disp(divider);
                end
            end   
        case 10
            close all;
            while(true)
                disp(' ');
                disp('Please, choose a method:');
                disp('     1 - Canny;');
                disp('     2 - Log;');
                disp('     3 - Prewitt;');
                disp('     4 - Roberts;');
                disp('     5 - Sobel;');
                disp('     6 - Zerocross.');
                disp('     0 - Go to previous menu.');
                option = input('Your method: ');
                switch option
                    case 1
                        close all;
                        figure('Name','Edges using canny','NumberTitle','off'), imshow(edge(image_gray,'Canny'));
                    case 2
                        close all;
                        figure('Name','Edges using log','NumberTitle','off'), imshow(edge(image_gray,'log'));
                    case 3
                        close all;
                        figure('Name','Edges using prewitt','NumberTitle','off'), imshow(edge(image_gray,'Prewitt'));
                    case 4
                        close all;
                        figure('Name','Edges using roberts','NumberTitle','off'), imshow(edge(image_gray,'Roberts'));
                    case 5
                        close all;
                        figure('Name','Edges using sobel','NumberTitle','off'), imshow(edge(image_gray,'Sobel'));
                    case 6
                        close all;
                        figure('Name','Edges using zerocross','NumberTitle','off'), imshow(edge(image_gray,'zerocross'));
                    case 0
                        close all;
                        break;
                    otherwise
                        close all;
                        disp('Please, insert a valid method.');
                disp(divider);
                end
            end
            disp(divider);
        case 11
            close all;
            figure('Name','Quadtree of processed binary image','NumberTitle','off'), imshow(blocks,[]);
            disp(divider);
        case 12
            close all;
            figure('Name','Derivative in x','NumberTitle','off'), imshow(mat2gray(abs(imgdx)));
            disp(divider);
        case 13
            close all;
            figure('Name','Derivative in y','NumberTitle','off'), imshow(mat2gray(abs(imgdy)));
            disp(divider);
        case 14
            close all;
            figure('Name','Module','NumberTitle','off'), imshow(mat2gray(modG));
            disp(divider);
        case 15
            close all;
            % Compute boundaries
            [B,L,N,A] = bwboundaries(bw);
            figure('Name','Boundaries of unprocessed binary image','NumberTitle','off'); imshow(bw);
            hold on;
            % Loop through object boundaries
            for k = 1:N
                % Boundary k is the parent of a hole if the k-th column
                % of the adjacency matrix A contains a non-zero element
                if (nnz(A(:,k)) > 0)
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
                    % Loop through the children of boundary k
                    for l = find(A(:,k))'
                        boundary = B{l};
                        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
                    end
                end
            end
            legend('External boundary','Inner boundary');
            disp(divider);
        case 16
            close all;
            % Compute boundaries
            [B,L,N,A] = bwboundaries(bw_final);
            figure('Name','Boundaries of processed binary image','NumberTitle','off'); imshow(bw_final);
            hold on;
            % Loop through object boundaries
            for k = 1 : N
                % Boundary k is the parent of a hole if the k-th column
                % of the adjacency matrix A contains a non-zero element
                if (nnz(A(:,k)) > 0)
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
                    % Loop through the children of boundary k
                    for l = find(A(:,k))'
                        boundary = B{l};
                        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
                    end
                end
            end
            legend('External boundary','Inner boundary');
            disp(divider);
        case 17
            close all;            
            R = imhist(image(:,:,1));
            G = imhist(image(:,:,2));
            B = imhist(image(:,:,3));
            figure('Name','RGB histogram','NumberTitle','off');
            plot(R,'r');
            hold on;
            plot(G,'g');
            plot(B,'b');
            legend('Red channel','Green channel','Blue channel');
            hold off;
            disp(divider);
        case 18
            close all;
            figure('Name','Grayscale histogram','NumberTitle','off'), imhist(image_gray);
            disp(divider);
        case 19
            close all;
            disp('Click on an area to see the object behind it.');
            figure('Name','See behind mask','NumberTitle','off'), imshow(bw_final), hold on;
            N = 1;
            but = 1;
            [ci,li,but] = ginput(1);
            if but == 1 % Left click
                cp(N) =  ci;
                lp(N) =  li;
            end
            object_to_find = lb(uint16(lp(1)),uint16(cp(1))); % Get the label
            image_object = (lb == object_to_find); % Get only the mask for the object of interest
            figure, imshow(image_gray.*uint8(image_object)); % And operation
            disp(divider);
        case 20
            close all;
            figure('Name','Hough transform','NumberTitle','off');
            subplot(2,1,1);
            imshow(image);
            subplot(2,1,2);
            imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,'InitialMagnification','fit');
            xlabel('\theta'), ylabel('\rho');
            axis on, axis normal, hold on;
            colormap(hot);
            disp(divider);
        case 21
            close all;
            figure('Name','Limited Theta Range Hough Transform of Gantrycrane Image','NumberTitle','off');
            imshow(imadjust(mat2gray(H_h)),'XData',T_h,'YData',R_h,'InitialMagnification','fit');
            xlabel('\theta'), ylabel('\rho');
            axis on, axis normal;
            colormap(hot);
            disp(divider);
        case 22
            close all;
            figure('Name','Gradient magnitude and gradient direction using Sobel method','NumberTitle','off');
            imshowpair(Gmag, Gdir, 'montage');
            axis off;
            disp(divider);
        case 23
            close all;
            figure('Name','Skeleton of the unprocessed binary image','NumberTitle','off'), imshow(bwmorph(bw,'skel',Inf));
        case 24
            close all;
            figure('Name','Skeleton of the processed binary image','NumberTitle','off'), imshow(bwmorph(bw_final,'skel',Inf));
        case 0
            close all;
            disp('Goodbye.');
            break;
        otherwise
            disp('Please, insert a valid option.');
            disp(divider);
    end
end