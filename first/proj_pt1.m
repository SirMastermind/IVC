clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load image
image = imread('Moedas1.jpg');
% imshow(image); title('Original');

%Convert image to gray
image_gray = rgb2gray(image);

% Use Otsu to get the threshold
thr = graythresh(image_gray)*255;
bw = image_gray>thr;
% figure,imshow(bw);

% Use closure to get the shapes well defined
se = strel('disk', 5);
bw_final = imclose(bw,se);
% figure,imshow(bw_final);

% Find and label the different regions
[lb, num]= bwlabel(bw_final);
%figure, imshow(mat2gray(lb));

% Get the area of each object
stats = regionprops(lb);
objects = [stats.Area];
objects = objects(objects > 1); % Take out single trash pixels
    
% Compute area for each region 
   
% Compute centroid for each region
    
% Compute perimeter for each region
    
% Compute distance for each object to other
    
% Compute differences between objects to find similarity
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Input and program's flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to the first project!');
disp(' ');
out = false;
while(true)
    % User interface
    disp('Please, input your option:');
    disp('     1 - Count objects;');
    disp('     2 - Show centroids;');
    disp('     3 - Show perimeters;');
    disp('     4 - Show areas;');
    disp('     5 - Show distance between objects;');
    disp('     6 - Show similar objects;');
    disp('     7 - Exit.');
    disp(' ');
    option = input('Your option: ');

    % Program's flow
    switch option
        case 1
            disp('The number of objects in the image is:');
            disp(length(objects));
        case 2
            close all;
            disp('Missing implementation')
        case 3
            close all;
            disp('Missing implementation')
        case 4
            close all;
            figure, imshow(mat2gray(lb));
        case 5
            disp('Missing implementation')
        case 6
            close all;
            disp('Missing implementation')
        case 7
            close all;
            disp('Googbye.');
            break;
    end
end