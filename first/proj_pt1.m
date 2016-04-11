clear all;
close all;

divider = sprintf('---------------');

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
for i = 1 : length(centroids)
    bw_centroids(uint16(centroids(i,1)), uint16(centroids(i,2))) = 0; % Mark the centroid with a black pixel
end

% Compute perimeter for each region
perimeters = bwperim(bw_final, 8); % bwperim(binary image, connectivity between neighbors)

% Compute distance for each object to other

% TODO
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Input and program's flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to the first project!');
disp(' ');
out = false;
while(true)
    % User interface
    disp(' ');
    disp('Please, input your option:');
    disp('     0 - Show original;');
    disp('     1 - Count objects;');
    disp('     2 - Show centroids;');
    disp('     3 - Show perimeters;');
    disp('     4 - Show areas;');
    disp('     5 - Show distance between objects;');
    disp('     6 - Show similar objects;');
    disp('     7 - Show edges;');
    disp('     8 - Show quadtree;');
    disp('     9 - Exit.');
    disp(' ');
    option = input('Your option: ');
    disp(' ');
    
    % Program's flow
    switch option
        case 0
            close all;
            figure, imshow(image);
            disp(divider);
        case 1
            string = sprintf('The number of objects in the image is %d.', length(objects));
            disp(string);
            disp(divider);
        case 2
            close all;
            for i = 1 : length(centroids)
                string = sprintf('Object %d has centroid in (%f, %f).', i, centroids(i,1), centroids(i,2));
                disp(string);
            end
            figure, imshow(bw_centroids);
            disp(divider);
        case 3
            close all;
            figure, imshow(perimeters); 
            disp('Missing implementation');
            disp(divider);
        case 4
            close all;
            figure, imshow(mat2gray(lb));
            disp(divider);
        case 5
            disp('Missing implementation');
            disp(divider);
        case 6
            close all;
            disp('Missing implementation');
            disp(divider);
        case 7
            close all;
            figure, imshow(edge(image_gray,'Prewitt'));
            disp(divider);
        case 8
            close all;
            figure, imshow(blocks,[]);
            disp(divider);
        case 9
            close all;
            disp('Googbye.');
            disp(divider);
            break;
    end
end