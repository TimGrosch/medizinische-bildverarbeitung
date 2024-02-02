%% import of data

% preparing the table
XL_table = readtable("patients_25.xlsx");
XL_table = XL_table(strcmp(XL_table.DatensatzVerwenden, "Y"),:);

% paths to save and load the data from
Path_Cases = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Cases";
Path_Matrices = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Matrices";
% for laptop
% Path_Cases = "C:\Users\Tim\Documents\MATLAB\Cases";
% Path_Matrices = "C:\Users\Tim\Documents\MATLAB\Matrices";

%% extracting data

% interpolates, cuts and saves all masks and volumes in the target path
prepare_all_patients(XL_table, Path_Cases, Path_Matrices);


%% Loading of example data

% enter the case id without the zeroes
Case_ID = 3;
added_zeros = 5 - length(num2str(Case_ID));
path = append(string(Case_ID),'.mat');
for i = 1:added_zeros
    path = append('0',path);
end

% loads the case
load(append(Path_Matrices, "\", path));

% extracts the image and mask of the target case
[example_coronal_layer, side_tumor] = get_data(Case_ID, XL_table);
example_image = squeeze(V(:,example_coronal_layer,:));
example_mask = squeeze(mask(:,example_coronal_layer,:));

right = 1:256;
left = 257:512;

if side_tumor(1) == 'r'
    side_cut = right;
else
    side_cut = left;
end



%% Visualisation

% Visualisation of 1 layer of the image
% and overlay of the mask over the image
display_segmented_images(example_image, example_mask);


%% Edge Detection

diff_image = imgaussfilt(example_image,2);
diff_image = diff_image(:,side_cut);

% three different edge detection algorithms
sob = rescale(sobel_filter(diff_image));
gpb = rescale(gPb(diff_image));
can = rescale(ImprovedCanny(diff_image,'rich'));

% figure
% subplot(3,1,1)
% title('Sobel')
% imshow(sob,[])
% subplot(3,1,2)
% title('Canny')
% imshow(can,[])
% subplot(3,1,3)
% title('general Probability of Boundary')
% imshow(gpb,[])

% for a solid edge detection, the three filter algorithms get combined
super_edge = compareEdges(sob, gpb, can);

% getting rid of artifacts and smaller edges
super_edge(super_edge > 0) = 1;
super_edge = bwmorph(super_edge,"thin",inf);

% Metric
metric = nnz(imfill(super_edge,"holes"))/nnz(super_edge);

if metric > 2.5
    super_edge = bwareaopen(super_edge,50);
else
    super_edge = clean_up(super_edge);
    super_edge = bwareaopen(super_edge,50);
    super_edge = bwmorph(super_edge,"thin",inf);
end

figure;
imshow(super_edge,[]);


%% Localization

% import of example reference
example_reference = rgb2gray(imread("KidneyCoronal.png"));

% cutting the img to only analyze the left kidney

%%%%%%%%%%% Choose which to use
%%%%%%%%%%% one fills up the edge of the kidney, the other does not
%%%%%%%%%%% depending on which changes have to be made in find_object

% img = imfill(rescale(super_edge(:,side_cut)),'holes');   % Filled
% img = rescale(super_edge(:,side_cut));                 % Not filled
img = imfill(rescale(super_edge),'holes');

tic
[target_marked,reference_marked,YBest,XBest,ang,scale,score] = find_object(img, example_reference);
toc

figure;
imshow(target_marked);




%% Chan-Vese segmentation of the Kidney

kidney_seg = activecontour(example_image(:,side_cut),reference_marked,'Chan-vese');
%%%%%
% delete every form except for largest
kidney_seg = bwareaopen(kidney_seg,1000);
%%%%%
kidney_seg = clean_up(kidney_seg,3);

figure
imshow(example_image(:,side_cut),[])
hold on
visboundaries(kidney_seg,'Color','r');

%% Sorensen-Dice similarity coefficient for image segmentation

cutted_mask = logical(example_mask(:,side_cut));
similarity_2D = dice(kidney_seg, cutted_mask);

figure;
subplot(3,1,1);
imshow(cutted_mask);
title('Mask');
subplot(3,1,2);
imshow(kidney_seg);
title('Kidney');
subplot(3,1,3);
imshowpair(kidney_seg, cutted_mask);
title(['Dice Index = ' num2str(similarity_2D)]);

%% Tumor detection (loop)

% import of example tumor reference



example_tumor_reference = rgb2gray(imread("Circle.png"));
tumor_score = 0;


tic
for n = 100:400

    tumor_image = squeeze(V(:,n,:));
    tumor_mask = squeeze(mask(:,n,:));
    cutted_tumor_image = tumor_image(:,side_cut);
    cutted_tumor_image(kidney_seg == 0) = 0;
    
    
    % tumor edge detection
    % same like the kidney

    cutted_tumor_image = imgaussfilt(cutted_tumor_image,2);
    sob = rescale(sobel_filter(cutted_tumor_image));
    gpb = rescale(gPb(cutted_tumor_image));
    can = rescale(ImprovedCanny(cutted_tumor_image,'rich'));
    super_tumor_edge = compareEdges(sob, gpb, can);
    
    super_tumor_edge(super_tumor_edge > 0) = 1;
    super_tumor_edge = bwmorph(super_tumor_edge,"thin",inf);
    metric = nnz(imfill(super_tumor_edge,"holes"))/nnz(super_tumor_edge);
    if metric > 2.5
        super_tumor_edge = bwareaopen(super_tumor_edge,50);
    else
        super_tumor_edge = clean_up(super_tumor_edge);
        super_tumor_edge = bwareaopen(super_tumor_edge,50);
        super_tumor_edge = bwmorph(super_tumor_edge,"thin",inf);
    end

    img = imfill(rescale(super_tumor_edge),'holes');

    % ght for each slice
    [target_marked,reference_marked,YBest,XBest,ang,scale,score] = find_object(img, example_tumor_reference);

    if score > tumor_score
        tumor_score = score
        target_slice = n
        tumor_marked = target_marked;
    end

end
toc

figure;
imshow(tumor_marked);

%% Tumor detection

% for case 3
chosen_slice = 270;
kmeans = 5;

% import of example tumor reference
example_tumor_reference = rgb2gray(imread("Circle.png"));

% cuts and filters one side of the image
tumor_image = squeeze(V(:,chosen_slice,:));
tumor_mask = squeeze(mask(:,chosen_slice,:));
cutted_tumor_image = tumor_image(:,side_cut);
cutted_tumor_image = imgaussfilt(cutted_tumor_image,2);

display_segmented_images(tumor_image, tumor_mask);


% cuts the image around the kidney
Best = [YBest, XBest];
% Calculate the bounding box area
shape = size(example_reference) * scale * 2;
% Compute the bounding box coordinates
area = [Best - shape/2; Best + shape/2];
area = max(area, 1);
% Ensure that the bounding box coordinates do not exceed the image size
area(:, 1) = min(area(:, 1), size(cutted_tumor_image, 1));
area(:, 2) = min(area(:, 2), size(cutted_tumor_image, 2));
% Create a binary mask for the tumor area within the bounding box
area_tumor = cutted_tumor_image(area(1, 1):area(2, 1), area(1, 2):area(2, 2));
cutted_tumor_image(area_tumor == 0) = 0;



% k-means-algorithm
cutted_tumor_image = single(area_tumor);
pre_segmented = imsegkmeans(cutted_tumor_image,kmeans);

figure;
subplot(2,1,1)
imshow(area_tumor,[])
subplot(2,1,2)
imshow(pre_segmented,[])

% Threshold the image
threshold = 1; % Adjust this threshold as needed
binary__kmeans = pre_segmented <= threshold;

figure;
imshow(binary__kmeans,[])

% tumor edge detection
% same like the kidney
sob = rescale(sobel_filter(binary__kmeans));
gpb = rescale(gPb(double(binary__kmeans)));
can = rescale(ImprovedCanny(binary__kmeans,'rich'));
super_tumor_edge = compareEdges(sob, gpb, can);
% getting rid of artifacts and smaller edges
super_tumor_edge(super_tumor_edge > 0) = 1;
super_tumor_edge = bwmorph(super_tumor_edge,"thin",inf);

% Metric
metric = nnz(imfill(super_tumor_edge,"holes"))/nnz(super_tumor_edge);
if metric > 2.5
    super_tumor_edge = bwareaopen(super_tumor_edge,50);
else
    % super_tumor_edge = clean_up(super_tumor_edge);
    super_edge = bwareaopen(super_edge,50);
    % super_edge = bwmorph(super_edge,"thin",inf);
end



figure;
subplot(4,1,1)
imshow(sob)
title('Sobel')
subplot(4,1,2)
imshow(gpb)
title('GPB')
subplot(4,1,3)
imshow(can)
title('Canny')
subplot(4,1,4)
imshow(super_tumor_edge)
title('super')

%%

img = imfill(rescale(super_tumor_edge),'holes');

[target_marked,reference_marked,Y_Best,X_Best,ang,scale,score] = find_object(img, example_tumor_reference);

figure;
imshow(target_marked);

