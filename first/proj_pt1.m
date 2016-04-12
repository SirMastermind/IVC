clear all;
close all;
warning off;

divider = sprintf('<----------------------------------------------->');

disp('Welcome to the first project!');
disp(' ');
disp('The program started processing your image...');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
for i = 1 : length(centroids) - 1
    bw_centroids(uint16(centroids(i,1)), uint16(centroids(i,2))) = 0; % Mark the centroid with a black pixel
end

% Compute perimeter for each region
perimeters = bwperim(bw_final, 8); % bwperim(binary image, connectivity between neighbors)

% Compute distance for each object to other
distance = zeros(size(centroids,1), size(centroids,1));
for i = 1 : size(distance,1) % For each lines
    for j = 1 : size(distance,2) % For each column
        distance(i,j) = sqrt((centroids(j,1) - centroids(i,1))^2 + (centroids(j,2) - centroids(i,2))^2); % Vector's length
        distance(j,i) = sqrt((centroids(j,1) - centroids(i,1))^2 + (centroids(j,2) - centroids(i,2))^2); % Mirror
    end
end
    
% Compute differences between objects to find similarity

% TODO

% Compute quadtree
output_size = [power(2,nextpow2(size(image_gray, 1))), power(2,nextpow2(size(image_gray, 2)))]; % Output size must be power of 2
S = qtdecomp(imresize(image_gray, output_size) , 0.27); % Resize it to the output size and decomposition in quadtree
blocks = repmat(uint8(0),size(S)); % Get the blocks

for dim = [512 256 128 64 32 16 8 4 2 1];    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)        
    values = repmat(uint8(1),[dim dim numblocks]);
    values(2:dim,2:dim,:) = 0;
    blocks = qtsetblk(blocks,S,dim,values);
  end
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Input and program's flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Image processed!');
disp(' ');
disp('This is the main menu.');

while(true)
    % User interface
    disp(' ');
    disp('Please, input your option:');
    disp('     1 - Show original;');
    disp('     2 - Count objects;');
    disp('     3 - Show centroids;');
    disp('     4 - Show perimeters;');
    disp('     5 - Show areas;');
    disp('     6 - Show distance between objects;');
    disp('     7 - Show similar objects;');
    disp('     8 - Show edges;');
    disp('     9 - Show quadtree;');
    disp('     10 - Show derivative in x;');
    disp('     11 - Show derivative in y;');
    disp('     12 - Show module;');
    disp('     13 - Show boundaries (unprocessed binary image);');
    disp('     14 - Show boundaries (processed binary image);');
    disp('     0 - Exit.');
    disp(' ');
    option = input('Your option: ');
    disp(' ');
    
    % Program's flow
    switch option
        case 1
            close all;
            figure, imshow(image);
            disp(divider);
        case 2
            string = sprintf('The number of objects in the image is %d.', length(objects));
            disp(string);
            disp(divider);
        case 3
            close all;
            for i = 1 : length(centroids) - 1
                string = sprintf('Object %d has centroid in (%f, %f).', i, centroids(i,1), centroids(i,2));
                disp(string);
            end
            figure, imshow(bw_centroids);
            disp(divider);
        case 4
            close all;
            figure, imshow(perimeters); 
            disp(divider);
        case 5
            close all;
            figure, imshow(mat2gray(lb));
            disp(divider);
        case 6
            close all;
            figure, imshow(bw_centroids), hold on;
            N = 0;
            but = 1;
            while (but == 1)
                [ci,li,but] = ginput(1);
                but;
                if but == 1 % Add point
                    N     = N+1;
                    cp(N) =  ci;
                    lp(N) =  li;
                    plot(ci,li,'r.','MarkerSize',8); drawnow;
                    if N > 1
                        plot(cp(:),lp(:),'r.-','MarkerSize',8); drawnow;
                    end
                end
                if size(cp,2) == 2
                    break;
                end
            end
%             object_1 = lb(uint16(lp(1)),uint16(cp(1)));
%             object_2 = [ lp(2) cp(2)];
            
            
            disp(divider);
        case 7
            close all;
            disp('Missing implementation');
            disp(divider);
        case 8
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
                        figure, imshow(edge(image_gray,'Canny'));
                    case 2
                        close all;
                        figure, imshow(edge(image_gray,'log'));
                    case 3
                        close all;
                        figure, imshow(edge(image_gray,'Prewitt'));
                    case 4
                        close all;
                        figure, imshow(edge(image_gray,'Roberts'));
                    case 5
                        close all;
                        figure, imshow(edge(image_gray,'Sobel'));
                    case 6
                        close all;
                        figure, imshow(edge(image_gray,'zerocross'));
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
        case 9
            close all;
            figure, imshow(blocks,[]);
            disp(divider);
        case 10
            close all;
            figure, imshow(mat2gray(abs(imgdx)));
            disp(divider);
        case 11
            close all;
            figure, imshow(mat2gray(abs(imgdy)));
            disp(divider);
        case 12
            close all;
            figure, imshow(mat2gray(modG));
            disp(divider);
        case 13
            close all;
            % Compute boundaries
            [B,L,N,A] = bwboundaries(bw);
            figure; imshow(bw); hold on;

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
            disp(divider);
        case 14
            close all;
            % Compute boundaries
            [B,L,N,A] = bwboundaries(bw_final);
            figure; imshow(bw_final); hold on;

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
            disp(divider);
        case 0
            close all;
            disp('Googbye.');
            disp(divider);
            break;
        otherwise
            disp('Please, insert a valid option.');
            disp(divider);
    end
end