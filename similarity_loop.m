%% 2D-Loop

% preparing the table
XL_table = readtable("patients_25.xlsx");
XL_table = XL_table(strcmp(XL_table.DatensatzVerwenden, "Y"),:);

similarity_table = table(Case_ID,similarity_2D)

%%
IDs = XL_table.CaseID;
num_cases = 5;

% paths to save and load the data from
Path_Cases = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Cases";
Path_Matrices = "C:\Users\Tim\Documents\MATLAB\Medizinische Bildverarbeitung\Matrices";

for i = 1:25

    ID = XL_table.CaseID(i,1);


    % Loading of example data
    % enter the case id without the zeroes
    Case_ID = ID
    added_zeros = 5 - length(num2str(Case_ID));
    path = append(string(Case_ID),'.mat');
    for j = 1:added_zeros
        path = append('0',path);
    end
    
    % loads the case
    load(append(Path_Matrices, "\", path));
    
    % extracts the image and mask of the target case
    example_coronal_layer = get_data(Case_ID, XL_table);
    example_image = squeeze(V(:,example_coronal_layer,:));
    example_mask = squeeze(mask(:,example_coronal_layer,:));
    
    % Edge Detection
    
    diff_image = imgaussfilt(example_image,2);
    diff_image = diff_image(:,257:end);
    
    % three different edge detection algorithms
    sob = rescale(sobel_filter(diff_image));
    gpb = rescale(gPb(diff_image));
    can = rescale(ImprovedCanny(diff_image,'rich'));
    
    % for a solid edge detection, the three filter algorithms get combined
    super_edge = compareEdges(sob, gpb, can);
    
    % getting rid of artifacts and smaller edges
    super_edge(super_edge > 0) = 1;
    super_edge = bwmorph(super_edge,"thin",inf);

    metric = nnz(imfill(super_edge,"holes"))/nnz(super_edge);
    
    if metric > 2.5
        super_edge = bwareaopen(super_edge,50);
    else
        super_edge = clean_up(super_edge);
        super_edge = bwareaopen(super_edge,50);
        super_edge = bwmorph(super_edge,"thin",inf);
    end
    
    % Localization
    
    % import of example reference
    example_reference = rgb2gray(imread("KidneyCoronal.png"));
    
    % cutting the img to only analyze the left kidney
    
    img = imfill(rescale(super_edge),'holes');
    
    tic
    [target_marked,reference_marked,YBest,XBest,ang,scale,score] = find_object(img, example_reference);
    toc

    % Chan-Vese segmentation of the Kidney
    
    kidney_seg = activecontour(example_image(:,257:end),reference_marked,'Chan-vese');
    kidney_seg = bwareaopen(kidney_seg,1000);
    kidney_seg = clean_up(kidney_seg,3);
    
    figure
    imshow(example_image(:,257:end),[])
    hold on
    visboundaries(kidney_seg,'Color','r');
    
    % Sorensen-Dice similarity coefficient for image segmentation
    
    cutted_mask = logical(example_mask(:,257:end));
    similarity_2D = dice(kidney_seg, cutted_mask)
    
    figure;
    subplot(3,1,1);
    title(['Mask']);
    imshow(cutted_mask);
    subplot(3,1,2);
    title('Kidney');
    imshow(kidney_seg);
    subplot(3,1,3);
    imshowpair(kidney_seg, cutted_mask);
    title(['Dice Index = ' num2str(similarity_2D)]);

    cell = {ID, similarity_2D};
    similarity_table = [similarity_table;cell];



end

filename = 'similarty_table.xlsx';
writetable(similarity_table,filename,'Sheet','MyNewSheet','WriteVariableNames',false);

