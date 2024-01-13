% Input
% image:    One layer of the volume
% mask:     One layer of the mask

function [] = display_segmented_images(image,mask)

%
figure;
subplot(2,1,1)
imshow(image)
subplot(2,1,2)
image = imfuse(image,mask);
imshow(image);

end