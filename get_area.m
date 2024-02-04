function[cut,area] = get_area(YBest,XBest,reference,image,scale)

Best = [YBest, XBest];
% Calculate the bounding box area
shape = size(reference) * scale;
% Compute the bounding box coordinates
area = [Best - shape/2; Best + shape/2];
area = max(area, 1);
% Ensure that the bounding box coordinates do not exceed the image size
area(:, 1) = min(area(:, 1), size(image, 1));
area(:, 2) = min(area(:, 2), size(image, 2));
% Create a binary mask for the tumor area within the bounding box
%area_tumor = zeros(size(cutted_tumor_image));
%area_tumor(area(1, 1):area(2, 1), area(1, 2):area(2, 2)) = 1;
cut = image(area(1, 1):area(2, 1), area(1, 2):area(2, 2));
end