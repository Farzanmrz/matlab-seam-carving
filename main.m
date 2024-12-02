%% Q2 - Image Resizing
clc; clearvars

% Read Images
im = imread('img.jpeg');
im2 = imread('img2.jpeg');

% Variable to store image to process based on user selection
% Variable to store original image
im_used = uint8.empty;
im_og = uint8.empty;

% Store image selection
im_selected = input("Enter 1 for image1, 2 for image2: ");

% Set image being processed as per input
if im_selected == 1
    im_og = im;
    im_used = im;
elseif im_selected == 2
    im_og = im2;
    im_used = im2;
else
    disp("Incorrect input");
    return;
end

% Get original height and width
[h, w, dim] = size(im_used);

% Ask user the new height
h_new = input("Enter new height: ");

% Ask user the new width
w_new = input("Enter new width: ");

% Declare variable to store resultant nn/bilinear result image
im_nn = uint8.empty;
im_bilinear = uint8.empty;

% Check if image is being downsampled
if h*w > h_new*w_new

    % Calculate how many times image being halved
    halvings = min(floor(log2(w/w_new)), floor(log2(h/h_new)));

    % Get kernel size adjusted for image being halved else default 3
    k = max(3, 2 * ceil(2 * 2^halvings) + 1);

    % Make sure the kernel is odd
    if mod(k, 2) == 0
        k = k + 1;
    end

    % Apply gaussian filtering 
    im_used = imgaussfilt(im_used, 2^halvings, 'FilterSize',k);
    
end

% Apply nearest neighbor
for y_new = 1:h_new
    for x_new = 1:w_new
        % Original Image ideal location
        x = round(x_new * (w/w_new));
        y = round(y_new * (h/h_new));
        
        
        % Assign the pixel value from nearest neighbor in the original image
        im_nn(y_new, x_new, :) = im_used(y, x, :);
    end
end

% Apply bilinear interpolation
for y_new = 1:h_new
    for x_new = 1:w_new

        % Original Image ideal location shift to avoid integer locations
        x = ((x_new-0.5) * (w/w_new))+0.5;
        y = ((y_new-0.5) * (h/h_new))+0.5;

        % Correct indices if they go out of bounds
        x = max(min(x, w), 1);
        y = max(min(y, h), 1);

        % Get x1,x2,y1,y2 for the four surrounding pixels
        x1 = floor(x);
        x2 = x1 + 1;
        y1 = floor(y);
        y2 = y1 + 1;

        % Correct for edge cases
        if x2 > w
            x2 = w;
            x1 = w - 1;
        end
        if y2 > h
            y2 = h;
            y1 = h - 1;
        end

        % Declare f(A),f(B),f(C),f(D) variables
        f_a = double(im_used(y1, x1, :));
        f_b = double(im_used(y1, x2, :));
        f_c = double(im_used(y2, x1, :));
        f_d = double(im_used(y2, x2, :));

        
        % Interpolate in x-direction
        f_xy1 = ((x2 - x) * f_a) + ((x - x1) * f_b);
        f_xy2 = ((x2 - x) * f_c) + ((x - x1) * f_d);

        % Interpolate in y-direction
        f_xy = ((y2 - y) * f_xy1) + ((y - y1) * f_xy2);

        im_bilinear(y_new, x_new, :) = uint8(f_xy);
    end
end

% Smooth afterward for upsampling
if h*w < h_new*w_new
    im_nn = imgaussfilt(im_nn, 2);
    im_bilinear = imgaussfilt(im_bilinear, 2);
end

% Show images
figure;
imshow(im_og)
figure;
imshow(im_nn);
figure;
imshow(im_bilinear);

%% Q3 - Energy Function
clc; clearvars

% Read images
im = imread('img.jpeg');
im2 = imread('img2.jpeg');

% Use 3a function to get magnitude/energy of images
im_energy = get_energy(im);
im2_energy = get_energy(im2);

% Show images in uint8
figure;imshow(uint8(im_energy));
figure;imshow(uint8(im2_energy));


%% Q4 - Optimal Seam
clc; clearvars;

% Read images
im = imread('img.jpeg');
im2 = imread('img2.jpeg');

% Use 3a func to get energy matrices
im_energy = get_energy(im);
im2_energy = get_energy(im2);

% Use 4a func to get cost matrices
im_c = get_cm(im_energy);
im2_c = get_cm(im2_energy);

% Use 4b func to get optimal seam vectors
im_os = get_os(im_c);
im2_os = get_os(im2_c);

% Use 4c func to get RGB image overlayed with optimal seam
im_seam = get_overlay(im,im_os,size(im,2));
im2_seam = get_overlay(im2,im2_os,size(im2,2));

% Display both images with seam superimposed
figure;imshow(im_seam);
figure;imshow(im2_seam);

%% Q5 - Seam Carving
clc; clearvars; 

% Read images
im = imread('img.jpeg');
im2 = imread('img2.jpeg');

% Get original width of image
im_width = size(im, 2);
im2_width = size(im2, 2);

% Output video showing seam carving for both images
get_video(im,im_width);
get_video(im2,im2_width);

%% Function Definitions

% Q3a - Return the gradient matrix: Input RGB matrix
function energy = get_energy(im)

    % Convert to grayscale then smooth with gaussian
    im_sm = imgaussfilt(rgb2gray(im), 2);

    % Get directional gradients
    [im_gx,im_gy] = imgradientxy(im_sm, 'sobel');

    % Get magnitude/energy of images
    energy = imgradient(im_gx,im_gy);

end

% Q4a - Return cost matrix. Input: Energy/Gradient matrix
function cost_matrix = get_cm(im_energy)
    
    % Initialize size of the cost matrix
    cost_matrix = zeros(size(im_energy,1),size(im_energy,2));

    % Set first row same as energy matrix
    cost_matrix(1,:) = im_energy(1,:);

    % Loop through all rows start from second since first is obtained
    for i = 2:size(im_energy,1)

        % Loop through all columns
        for j = 1:size(im_energy,2)
            
            % If left-edge then 
            if j == 1
    
                % Use up-middle and up-right pixel
                cost_matrix(i,j) = im_energy(i,j) + min(cost_matrix(i-1,j),cost_matrix(i-1,j+1));
            
            % If right-edge then
            elseif j == size(cost_matrix,2)
    
                % Use up-left and up-middle pixel
                cost_matrix(i,j) = im_energy(i,j) + min(cost_matrix(i-1,j-1),cost_matrix(i-1,j));
    
            % If no edge cases default
            else
    
                % Use up-left, up-middle and up-right pixel
                cost_matrix(i,j) = im_energy(i,j) + min([cost_matrix(i-1,j-1),cost_matrix(i-1,j),cost_matrix(i-1,j+1)]);
            
            end
    
        end
    end


end

% Q4b - Return optimal seam vector. Input: Cost matrix
function seam = get_os(im_c)

    % Initial seam vector
    seam = zeros(size(im_c,1),1);

    % Get the minimum value,column index of the least element in last row
    [~, min_idx] = min(im_c(end,:));

    % Loop through all rows start from last decrement by 1 till second row
    for i = size(im_c,1):-1:2
    
        % Set 1 on minimum value index
        seam(i) = min_idx;
    
        % If left edge
        if min_idx == 1
    
            % Take middle and left
            indices = [min_idx, min_idx+1];
    
        % If right edge
        elseif min_idx == size(im_c,2)
    
            % Take left and middle
            indices = [min_idx-1, min_idx];
    
        % Not edge
        else
    
            % Take all left,middle and right
            indices = [min_idx-1, min_idx, min_idx+1];
        end
        
        % Find minimum value in previous row using the indices found
        [~, curr_idx] = min(im_c(i-1, indices));
    
        % Update minimum index
        min_idx = indices(curr_idx);
    
    end

    % Mark the first row pixel with last value of min_idx
    seam(1) = min_idx;
    
end

% Q4c - Return padded RGB image with optimal seam overlayed
% Input: RGB matrix, Optimal seam vector
function seam_image = get_overlay(im, im_os, im_width)

    % Initialize overlayed image matrix
    seam_image = im;

    % Set the red channel to max value rest all to 0
    for i = 1:length(im_os)
        seam_image(i, im_os(i), 1) = 255;  
        seam_image(i, im_os(i), 2) = 0;
        seam_image(i, im_os(i), 3) = 0;    
    end 

    % Get current width to pad 0s
    curr_width = size(seam_image,2);

    % Check if current width less than original width 
    % so everywhere except start
    if curr_width < im_width

        % Add column of 0s at the current column next to the produced image
        seam_image = [seam_image, uint8(zeros(size(seam_image,1),im_width - curr_width, 3))];
    end

end

% Q5a - Remove seam and return new image matrix 
% Input: RGB matrix, Optimal seam vector
function new_image = remove_seam(im, im_os)

    % Initialize new image size 1 pixel less wide
    new_image = uint8(zeros(size(im,1),size(im,2)-1,size(im,3)));

    % Go through all rows of the image
    for i = 1:size(im, 1)

        % Keep every column besides the one with seam pixel on all channels
        new_image(i, :, :) = im(i, [1:im_os(i)-1, im_os(i)+1:end], :);
    end

end

% Q5b - Write video to file. Input: RGB matrix, original width
function vid = get_video(im,im_width)
    
    % Set output filename to the variable name for smooth running
    name = [inputname(1), '.mp4'];
    
    % Create a VideoWriter object
    im_vid = VideoWriter(name, 'MPEG-4');
    im_vid.FrameRate = 30;
    
    % Open the VideoWriter for writing
    open(im_vid);
    
    % Loop till width - 1 since have to be 1 pixel wide
    for i = 1: im_width - 1
        
        % Use 3a func to get energy matrix
        im_energy = get_energy(im);

        % Use 4a func to get cost matrix
        im_c = get_cm(im_energy);

        % Use 4b func to get optimal seam matrix
        im_os = get_os(im_c);

        % Use 4c func to get RGB image overlayed with optimal seam
        im_seam = get_overlay(im,im_os,im_width);

        % Write the overlayed image matrix to video
        writeVideo(im_vid, im_seam);

        % Use 5a to do seam removal
        im_new = remove_seam(im, im_os);

        % Set the new image to the original image
        im = im_new;
    end

    % Close the output video
    close(im_vid);

end