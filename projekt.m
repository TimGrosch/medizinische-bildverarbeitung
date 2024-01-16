%% import of data

% preparing the table
XL_table = readtable("patients_25.xlsx");
XL_table = XL_table(strcmp(XL_table.DatensatzVerwenden, "Y"),:);


%% extracting data

% paths to save and load the data from
Path_Cases = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Cases";
Path_Matrices = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Matrices";

% interpolates, cuts and saves all masks and volumes in the target path
prepare_all_patients(XL_table, Path_Cases, Path_Matrices);


%% Loading of example data

% enter the case id without the zeroes
Case_ID = 132;
added_zeros = 5 - length(Case_ID);
path = append(string(Case_ID),'.mat');
for i = 1:added_zeros
    path = append('0',path);
end

% extracts the image and mask of the target case
example_coronal_layer = get_data(Case_ID, XL_table);
example_image = squeeze(V(:,example_coronal_layer,:));
example_mask = squeeze(mask(:,example_coronal_layer,:));

%% Visualisation

% Visualisation of 1 layer of the image
% and overlay of the mask over the image
display_segmented_images(example_image, example_mask);


%% Edge Detection

diff_image = imgaussfilt(example_image,2);
diff_image = diff_image(:,257:end);

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
%%%%%%
% Metric
metric = nnz(imfill(super_edge,"holes"))/nnz(super_edge);

if metric > 2.5
    super_edge = bwareaopen(super_edge,50);
else
    super_edge = clean_up(super_edge);
    super_edge = bwareaopen(super_edge,50);
    super_edge = bwmorph(super_edge,"thin",inf);
end

% edges_cut = 0;
% counter = 1;
% 
% while ~edges_cut
%     edge = super_edge - bwmorph(super_edge,'endpoints');
%     if edge == super_edge
%         edges_cut = 1;
%     end
%     super_edge = edge;
%     counter = counter + 1;
% end

figure;
imshow(super_edge,[]);


%% Localization

% import of example reference
example_reference = rgb2gray(imread("KidneyCoronal.png"));

% cutting the img to only analyze the left kidney

%%%%%%%%%%% Choose which to use
%%%%%%%%%%% one fills up the edge of the kidney, the other does not
%%%%%%%%%%% depending on which changes have to be made in find_object

% img = imfill(rescale(super_edge(:,257:end)),'holes');   % Filled
% img = rescale(super_edge(:,257:end));                 % Not filled
img = imfill(rescale(super_edge),'holes');

tic
[target_marked,reference_marked,YBest,XBest,ang,scale,score] = find_object(img, example_reference);
toc

figure;
imshow(target_marked);

%% Chan-Vese segmentation of the Kidney

kidney_seg = activecontour(example_image(:,257:end),reference_marked,'Chan-vese');
%%%%%
% delete every form except for largest
kidney_seg = bwareaopen(kidney_seg,1000);
%%%%%
kidney_seg = clean_up(kidney_seg,3);
figure
imshow(example_image(:,257:end),[])
hold on
visboundaries(kidney_seg,'Color','r'); 
