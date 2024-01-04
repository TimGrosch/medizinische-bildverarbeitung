
%Step 1: Datenimport

%Constants
PATHXCLTABLE = "C:\Users\Tim\Documents\MATLAB\patients_25.xlsx";
PATHCASES = "C:\Users\Tim\Documents\MATLAB\Cases\case_";

%Read Data
data = readtable(PATHXCLTABLE);
patientData = data(strcmp(data.DatensatzVerwenden, "Y"),:);
patientIDs = patientData.CaseID;
layers_x = patientData.NiereMaximalBeiKoronalerSchicht_x_;
layers_y_min = patientData.gew_hlteAxialeSchichten_z_;
layers_y_max = patientData.Var11;
% Ich habe die Datensätze limitiert, damit das Skript schneller läuft
amountPatients = 5; %length(patientIDs);

cases = cell(1,amountPatients);
masks = cell(1,amountPatients);

for i = 1:amountPatients
    id = patientIDs(i,1);
    layer_x = layers_x(i,1);
    layer_y_min = layers_y_min(i,1);
    layer_y_max = layers_y_max(i,1);
    id_str = string(id);
    id_len = strlength(id_str);
    adding_zeros = 5 - id_len;

    for n = 1:adding_zeros
        id_str = append("0", id_str);
    end

    finalPathImages = append(PATHCASES, id_str, "\imaging.nii.gz");
    finalPathMasks = append(PATHCASES, id_str, "\segmentation.nii.gz"); 

    cases{i} = double(niftiread(finalPathImages));
    masks{i} = double(niftiread(finalPathMasks));
    
    %Function_show2DPlot(cases{i}, masks{i}, layer_x, layer_y_min, layer_y_max)

    data_row = patientData(i,:);
    interpolated = interpolate_nifti_z(cases{i}, data_row);
    img = Function_show2DPlot(interpolated, masks{i}, layer_x, layer_y_min, layer_y_max);


    canny = edge(img,"canny");
    %Canny Kantendetektion
    figure();
    subplot(2,2,1)
    imshow(canny)
    subplot(2,2,2)
    sobel = edge(img,"sobel");
    imshow(sobel)
    subplot(2,2,3)
    prewitt = edge(img,  "prewitt");
    imshow(prewitt)

end




%To-Dos:


