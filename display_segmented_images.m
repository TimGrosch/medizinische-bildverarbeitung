function [] = display_segmented_images(image,mask)

image = imfuse(image,mask);
imshow(image);

end