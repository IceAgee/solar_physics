close all; clear all; clc

load('full_disk_sub_map.mat')
data = double(data);  % Convert the data to double precision
norm_data1 = (data - nanmin(data)) ./ (nanmax(data) - nanmin(data)) * 255;  % Normalize the data to 0-255
norm_data2 = norm_data1;  % Backup the normalized data
% Apply a nonlinear logarithmic transformation to enhance the contrast
norm_data1 = 30 * (log10(norm_data1 + 1)) + 0.3 .* (norm_data1 + 0.01).^0.1;
I = uint8(255 - ((norm_data1)));  % Invert the image (bright becomes dark, dark becomes bright)
w = fspecial('gaussian',[42 42], 10);  % Create a 42x42 Gaussian filter with a standard deviation of 10
I = imfilter(I, w);  % Apply the Gaussian filter to the image

gmag = imgradient(I);  % Compute the gradient magnitude of the image (used for edge detection)
se = strel('disk', 5);  % Create a disk-shaped structural element with radius 5
Io = imopen(I, se);  % Apply morphological opening (removes small bright regions)
Ie = imerode(I, se);  % Apply morphological erosion
Iobr = imreconstruct(Ie, I);  % Perform morphological reconstruction after erosion
Ioc = imclose(Io, se);  % Apply morphological closing (fills small dark regions)
Iobrd = imdilate(Iobr, se);  % Apply dilation to the image
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));  % Perform reconstruction using the complement of the dilated image
Iobrcbr = imcomplement(Iobrcbr);  % Complement the result to restore the original image direction

fgm = imregionalmax(Iobrcbr);  % Detect regional maxima (used for foreground marking before watershed)
I2 = labeloverlay(I, fgm);  % Overlay the foreground maxima on the image
se2 = strel(ones(5, 5));  % Create a 5x5 structural element (for further morphological operations)
fgm2 = imclose(fgm, se2);  % Apply morphological closing to the foreground maxima
fgm3 = imerode(fgm2, se2);  % Apply erosion to the foreground maxima
fgm4 = bwareaopen(fgm3, 20);  % Remove small foreground regions (less than 20 pixels)
I3 = labeloverlay(I, fgm4);  % Overlay the cleaned foreground maxima on the image

bw = imbinarize(Iobrcbr);  % Binarize the image (convert to black and white)
D = bwdist(bw);  % Compute the Euclidean distance transform (distance to the nearest background pixel)
DL = watershed(D);  % Apply the watershed algorithm to the distance transform
bgm = DL == 0;  % Mark the background as 0 (edges of the watershed)

gmag2 = imimposemin(gmag, bgm | fgm4);  % Impose minima on the gradient magnitude to guide the watershed
L = watershed(gmag2);  % Apply the final watershed segmentation
labels = imdilate(L == 0, ones(3,3)) + 2 * bgm + 3 * fgm4;  % Create a label map for segmentation
I4 = labeloverlay(I, labels);  % Overlay the segmentation labels on the original image

close all
figure

% Modify the watershed labels for boundary highlighting
Lnew = L;
Lnew(Lnew > 0) = 4;  % Set all non-boundary regions to 4
Lnew(Lnew == 0) = 255;  % Set boundary regions to 255 (white)
% Remove boundary pixels along the edges of the image (prevent false detection)
size_Lnew = size(Lnew);
Lnew(1:3,:) = 4;
Lnew(size_Lnew(1)-2:size_Lnew(1),:) = 4;
Lnew(:,1:3) = 4;
Lnew(:,size_Lnew(2)-2:size_Lnew(2)) = 4;

Lnew(Lnew == 4) = 0;  % Set non-boundary regions back to 0 (black)
boundary_index = find(Lnew == 255);  % Find the boundary pixel indices

boundary_width = 1;  % Set the boundary width to 1 pixel
% Expand the boundary region by 1 pixel
for ii = 1:boundary_width
    Lnew(boundary_index + ii) = 255;
    Lnew(boundary_index - ii) = 255;
    Lnew(boundary_index + ii * size_Lnew(1)) = 255;
    Lnew(boundary_index - ii * size_Lnew(1)) = 255;
end

Lnew = uint8(Lnew);  % Convert the label map to uint8 type (for overlaying on the image)

save('Lnew_full_sub_map.mat', "Lnew", '-mat')  % Save the boundary label map

% Overlay the boundary labels on the original image for visualization
himage = pcolor(uint8(norm_data2) + Lnew);
himage.LineStyle = 'none';  % Remove the grid lines from the plot
title('Chromospheric Network Boundaries')  % Set the title of the plot
