%% import of data

% preparing the table
XL_table = readtable("patients_25.xlsx");
XL_table = XL_table(strcmp(XL_table.DatensatzVerwenden, "Y"),:);

% paths to save and load the data from
Path_Cases = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Cases";
Path_Matrices = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Matrices";

% interpolates, cuts and saves all masks and volumes in the target path
prepare_all_patients(XL_table, Path_Cases, Path_Matrices);


%% Loading of example data

% only valid for case 00003
example_coronal_layer = 305
example_image = squeeze(V(:,example_coronal_layer,:));
example_mask = squeeze(mask(:,example_coronal_layer,:));

%% Visualisation

% Visualisation of 1 layer of the image
% and overlay of the mask over the image
display_segmented_images(example_image, example_mask);


%% Edge Detection

img = rescale(example_image(:,1:256));
img = imgaussfilt(img,2);

% three different edge detection algorithms
sob = rescale(sobel_filter(img));
gpb = rescale(gPb(img));
can = rescale(ImprovedCanny(img,'rich'));

figure
subplot(3,1,1)
title('Sobel')
imshow(sob,[])
subplot(3,1,2)
title('Canny')
imshow(can,[])
subplot(3,1,3)
title('general Probability of Boundary')
imshow(gpb,[])

% for a solid edge detection, the three filter algorithms get combined
super_edge = compareEdges(sob, gpb, can);

% getting rid of artifacts and smaller edges
super_edge = bwareaopen(super_edge,50);

figure;
imshow(super_edge,[]);



