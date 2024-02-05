%% import of data

% preparing the table
XL_table = readtable("patients_25.xlsx");
XL_table = XL_table(strcmp(XL_table.DatensatzVerwenden, "Y"),:);

% paths to save and load the data from
%Path_Cases = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Cases";
%Path_Matrices = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Matrices";
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
%load(append(Path_Matrices, "\", path));
load(path);

% extracts the image and mask of the target case
[example_coronal_layer, side_tumor] = get_data(Case_ID, XL_table);
side_tumor = 'links';
example_image = squeeze(V(:,example_coronal_layer,:));
example_mask = squeeze(mask(:,example_coronal_layer,:));

right = 1:256;
left = 257:512;

if side_tumor(1) == 'r'
    side_cut = right;
else
    side_cut = left;
end

clear i added_zeros path side_tumor


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

figure
subplot(2,2,1)
imshow(sob,[])
title('Sobel')
subplot(2,2,2)
imshow(can,[])
title('Canny')
subplot(2,2,3)
imshow(gpb,[])
title('general Probability of Boundary')


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

subplot(2,2,4)
imshow(super_edge,[]);
title('Super Edge')


clear diff_image gpb sob can metric

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
imshow(img,[]);
hold
visboundaries(reference_marked,'Color','b')
title('Localization of the kidney')

clear img ang score



%% Chan-Vese segmentation of the Kidney

kidney_seg = activecontour(example_image(:,side_cut),reference_marked,'Chan-vese','SmoothFactor',0.5);
%%%%%
% delete every form except for largest
%%%%NEWLY ADDED 03.02
max_size = 0;
list_objects = bwconncomp(kidney_seg);
list_objects = list_objects.PixelIdxList;
for j = 1:length(list_objects)
    arr = list_objects(j);
    curr = length(arr{1});
    if curr > max_size
        max_size = curr;
    end
end
kidney_seg = bwareaopen(kidney_seg, max_size-1);
%%%%%
%kidney_seg = clean_up(kidney_seg,3);
kidney_seg = imfill(kidney_seg,'holes');

figure
imshow(example_image(:,side_cut),[])
hold on
visboundaries(kidney_seg,'Color','r');
title('Segmentation of the kidney')

clear arr curr j i list_objects max_size

%% Sorensen-Dice similarity coefficient for image segmentation

cutted_mask = logical(example_mask(:,side_cut));
similarity_2D = dice(kidney_seg, cutted_mask);

figure;
subplot(3,1,1);
imshow(cutted_mask);
title('Provided Kidney Mask');
subplot(3,1,2);
imshow(kidney_seg);
title('Segmented Kidney');
subplot(3,1,3);
imshowpair(kidney_seg, cutted_mask);
title(['Dice Index = ' num2str(similarity_2D)]);

%%
similarity = dice(mask,mask_segmented)
%% SKIP (Tumor detection (loop))

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

%% Expansion of the kidney mask into 3D
%%%%%%%NEWLY ADDED 02.03.
generated_mask_right = expansion_volume(kidney_seg,example_coronal_layer,V,'left');
generated_mask_right(:,example_coronal_layer,:) = kidney_seg;

%% Tumor detection

% for case 3
chosen_slice = 270;
kmeans = 4;

% import of example tumor reference
example_tumor_reference = rgb2gray(imread("Circle.png"));

% cuts and filters one side of the image
tumor_image = squeeze(V(:,chosen_slice,:));
tumor_mask = squeeze(mask(:,chosen_slice,:));
cutted_tumor_image = tumor_image(:,side_cut);
cutted_tumor_image = imgaussfilt(cutted_tumor_image,2);

%display_segmented_images(tumor_image, tumor_mask);


% cuts the image around the kidney
[tumor_focused, area_tumor] = get_area(YBest,XBest,example_reference,cutted_tumor_image,scale*2.5);


% k-means-algorithm
cutted_tumor_image = single(tumor_focused);
%cutted_tumor_image = single(cutted_tumor_image);
pre_segmented = imsegkmeans(cutted_tumor_image,kmeans);

figure;
subplot(1,2,1)
imshow(tumor_focused,[])
title('Cutout of slice')
subplot(1,2,2)
imshow(pre_segmented,[])
title(append('k-means segmentated (Number of clusters = ',num2str(kmeans),')'))

%%  SKIP
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

super_tumor_edge = bwareaopen(super_tumor_edge,50);


%%%%NEWLY ADDED 03.02.
close all
super_tumor_edge = bwmorph(super_tumor_edge,'close');
super_tumor_edge = bwmorph(super_tumor_edge,'thin',inf);
super_tumor_edge = bwareaopen(super_tumor_edge,20);

figure;
subplot(2,2,1)
imshow(sob)
title('Sobel')
subplot(2,2,2)
imshow(gpb)
title('GPB')
subplot(2,2,3)
imshow(can)
title('Canny')
subplot(2,2,4)
imshow(super_tumor_edge)
title('super')

clear sob can gpb 

%% SKIP

img = imfill(rescale(super_tumor_edge),'holes');
%img = super_tumor_edge;
%img(img > 0) = 1;

[target_marked_t,rmt,Y_Best_t,X_Best_t,ang_t,scale_t,score_t] = find_object(img, example_tumor_reference);

figure;
imshow(target_marked_t);

%% SKIP
target_marked_t = zeros(size(tumor_image(:,side_cut)));
reference_marked_t( ~any(reference_marked_t,2), : ) = [];  %rows
reference_marked_t( :, ~any(reference_marked_t,1) ) = [];  %columns

[cut, area_t] = get_area(Y_Best_t+area_tumor(1,1),X_Best_t+area_tumor(1,2),reference_marked_t,tumor_image(:,side_cut),1);

target_marked_t(area_t(1, 1):area_t(2, 1)-1, area_t(1, 2):area_t(2, 2)-1) = reference_marked_t;

%% SKIP
tumor_seg = activecontour(tumor_image(:,side_cut),target_marked_t,'edge','SmoothFactor',0.8,'ContractionBias',-0.3);
%tumor_seg = activecontour(tumor_image(:,side_cut),target_marked_t,'edge');

figure
imshow(tumor_image(:,side_cut),[])
hold on
visboundaries(tumor_seg,'Color','r');
visboundaries(target_marked_t,'Color','b');


%% Idea of segmenting by kmeans and kicking the tissue of kidney

binary_masks_tumor = zeros([136,kmeans,111]);
cor_coef = zeros([4,2]);
kidney_mask = squeeze(generated_mask_right(:,chosen_slice,:));

for i = min(min(pre_segmented)):max(max(pre_segmented))
    inbetween = pre_segmented;
    inbetween(inbetween ~= i) = 0;
    binary_masks_tumor(:,i,:) = inbetween;
    cor_coef(i,1) = i;
    % cor_coef(i,2) = max(max(xcorr2(kidney_mask,inbetween)));
    cor_coef(i,2) = mean2(xcorr2(kidney_mask,inbetween));
end

% figure
% for i = 1:4
% subplot(2,2,i)
% imshow(squeeze(binary_masks_tumor(:,i,:)),[])
% end

% for i = 1:2
[to_delete, cor] = max(cor_coef);
cor_coef(cor_coef(:,1) == to_delete,:) = [];
binary_masks_tumor = binary_masks_tumor(:,cor_coef(:,1),:);
% end
% figure
% for i = 1:2
% subplot(1,3,i)
% imshow(squeeze(binary_masks_tumor(:,i,:)),[])
% end

%% Localization of the tumor

% in this case the score matches the desired output but needs to be adjust
% if more cases need higher robustness
close all
example_tumor_reference = edge(imresize(rgb2gray(imread("Circle.png")),0.5),'sobel');
bestscore = 0;
variables_localization_tumor = zeros(3);

for i = 1:3
    edge_for_bin = edge(squeeze(binary_masks_tumor(:,i,:)),"sobel");
    [ttm,trm,tY,tX,tang,tscale,tscore] = find_object(edge_for_bin, example_tumor_reference);
    if tscore > 2
        % reference_marked_t = trm;
        variables_localization_tumor(1,i) = tY;
        variables_localization_tumor(2,i) = tX;
        variables_localization_tumor(3,i) = tscore;
        % Y_Best_t = tY;
        % X_Best_t = tX;
        % bestscore = tscore;
    end
    figure
    imshow(ttm,[])
    title(tscore)
end

%%
d = length(binary_masks_tumor);
for i = 1:length(variables_localization_tumor)-1
    for j = i+1:length(variables_localization_tumor)
        a = [variables_localization_tumor(1,i) variables_localization_tumor(2,i)];
        b = [variables_localization_tumor(1,j) variables_localization_tumor(2,j)];
        
        if norm(b-a) < d
            d = norm(b-a);
            tissue = i;
        end
    end
end

[ttm,reference_marked_t,Y_Best_t,X_Best_t,tang,tscale,tscore] = find_object(edge(squeeze(binary_masks_tumor(:,tissue,:)),'sobel'), example_tumor_reference);
reference_marked_t = imfill(reference_marked_t,'holes');


%% Putting the mask in the same size as the original image again,
% (It has been cut to the area around the kidney)

target_marked_t = zeros(size(tumor_image(:,side_cut)));
reference_marked_t( ~any(reference_marked_t,2), : ) = [];  %rows
reference_marked_t( :, ~any(reference_marked_t,1) ) = [];  %columns

[cut, area_t] = get_area(Y_Best_t+area_tumor(1,1),X_Best_t+area_tumor(1,2),reference_marked_t,tumor_image(:,side_cut),1);

target_marked_t(area_t(1, 1):area_t(2, 1)-1, area_t(1, 2):area_t(2, 2)-1) = reference_marked_t;

%% Segmentation of the tumor 
tumor_seg = activecontour(tumor_image(:,side_cut),target_marked_t,'edge','SmoothFactor',0.8,'ContractionBias',-0.3);
%tumor_seg = activecontour(tumor_image(:,side_cut),target_marked_t,'edge');

figure
imshow(tumor_image(:,side_cut),[])
hold on
visboundaries(tumor_seg,'Color','r');
%visboundaries(target_marked_t,'Color','b');
title('Tumor segmented')

%% 3D expansion of tumor
% 257 - 289

% TODO: Cutoffs don't work, better one to find or arbitrary for-loop

generated_mask_right_tumor = zeros(size(generated_mask_right));
generated_mask_right_tumor(:,chosen_slice,:) = tumor_seg;
current_slice = chosen_slice;
new_seg = tumor_seg;
scale_ref = 1.75;
example_tumor_reference = imresize(rgb2gray(imread("Circle.png")),0.5);

while true
    current_slice = current_slice - 1;
    tic
    % [tm,rm,Y_tumor,X_tumor,ant,sct,scot] = find_object(new_seg, example_tumor_reference);
    % 
    % %
    % rm = imresize(rm,scale_ref);
    % rm( ~any(rm,2), : ) = [];  %rows
    % rm( :, ~any(rm,1) ) = [];  %columns
    % 
    % % rm = imresize(example_tumor_reference,sct);
    % % rm( ~any(rm,2), : ) = [];  %rows
    % % rm( :, ~any(rm,1) ) = [];  %columns
    % 
    % tm = zeros(size(tm));
    % [cut, area_t] = get_area(Y_tumor,X_tumor,rm,tm,1);
    % 
    % tm(area_t(1, 1):area_t(2, 1)-1, area_t(1, 2):area_t(2, 2)-1) = rm;
    % tm = bwmorph(tm,'bridge');
    % tm = imfill(tm,'holes');

    inbetween = activecontour(squeeze(V(:,current_slice,side_cut)),new_seg,'edge','SmoothFactor',0.8,'ContractionBias',0);
    toc
    % if scale_ref > 1 && nnz(inbetween) < nnz(new_seg)*0.85
    %     scale_ref = scale_ref * 0.95;
    % end
    if nnz(inbetween) > 2*nnz(new_seg)
        fprintf('Tumor stops at slice: %d\n',current_slice)
        break;
    elseif nnz(inbetween) == 0
        fprintf('Tumor stops at slice: %d\n',current_slice)
        break;
    else
        new_seg = inbetween;
        generated_mask_right_tumor(:,current_slice,:) = new_seg;
    end

    close all
    imshow(squeeze(V(:,current_slice,side_cut)),[])
    hold on
    visboundaries(new_seg,'Color','r');
    hold off
end

current_slice = chosen_slice;
new_seg = tumor_seg;
scale_ref = 1;

while true
    current_slice = current_slice + 1;
    tic
    % [tm,rm,Y_tumor,X_tumor,ant,sct,scot] = find_object(new_seg, example_tumor_reference);
    % 
    % %
    % rm = imresize(rm,scale_ref);
    % rm( ~any(rm,2), : ) = [];  %rows
    % rm( :, ~any(rm,1) ) = [];  %columns
    % 
    % % rm = imresize(example_tumor_reference,sct);
    % % rm( ~any(rm,2), : ) = [];  %rows
    % % rm( :, ~any(rm,1) ) = [];  %columns
    % 
    % tm = zeros(size(tm));
    % [cut, area_t] = get_area(Y_tumor,X_tumor,rm,tm,1);
    % 
    % tm(area_t(1, 1):area_t(2, 1)-1, area_t(1, 2):area_t(2, 2)-1) = rm;
    % tm = bwmorph(tm,'bridge');
    % tm = imfill(tm,'holes');

    inbetween = activecontour(squeeze(V(:,current_slice,side_cut)),new_seg,'edge','SmoothFactor',0.8,'ContractionBias',0);
    toc
    % tic
    % [tm,rm,Y_tumor,X_tumor,ant,sct,scot] = find_object(new_seg, example_tumor_reference);
    % 
    % rm = imresize(rm,scale_ref);
    % rm( ~any(rm,2), : ) = [];  %rows
    % rm( :, ~any(rm,1) ) = [];  %columns
    % 
    % tm = zeros(size(tm));
    % [cut, area_t] = get_area(Y_tumor,X_tumor,rm,tm,1);
    % 
    % tm(area_t(1, 1):area_t(2, 1)-1, area_t(1, 2):area_t(2, 2)-1) = rm;
    % 
    % inbetween = activecontour(squeeze(V(:,current_slice,side_cut)),tm,'edge','SmoothFactor',0.8,'ContractionBias',0.1);
    % toc
    % if scale_ref > 1 && nnz(inbetween) < nnz(new_seg)*0.85
    %     scale_ref = scale_ref * 0.95;
    % end
    if nnz(inbetween) > 2*nnz(new_seg)
        fprintf('Tumor stops at slice: %d\n',current_slice)
        break;
    elseif nnz(inbetween) == 0
        fprintf('Tumor stops at slice: %d\n',current_slice)
        break;
    else
        new_seg = inbetween;
        generated_mask_right_tumor(:,current_slice,:) = new_seg;
    end

    close all
    imshow(squeeze(V(:,current_slice,side_cut)),[])
    hold on
    visboundaries(new_seg,'Color','r');
    hold off
end

%%

