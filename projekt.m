
%Step 1: Datenimport

%Constants
PATHXCLTABLE = "C:\Users\Tim\Documents\MATLAB\patients_25.xlsx";
PATHCASES = "C:\Users\Tim\Documents\MATLAB\Cases\case_";

%Read Data
data = readtable(PATHXCLTABLE);
patientData = data(strcmp(data.DatensatzVerwenden, "Y"),:);
patientIDs = patientData.CaseID;
% Ich habe die Datensätze limitiert, damit das Skript schneller läuft
amountPatients = 3 %length(patientIDs)

cases = cell(1,amountPatients);
masks = cell(1,amountPatients);
for i = 1:amountPatients
    id = patientIDs(i,1);
    id_str = string(id);
    id_len = strlength(id_str);
    adding_zeros = 5 - id_len;

    for n = 1:adding_zeros
        id_str = append("0", id_str);
    end

    finalPathImages = append(PATHCASES, id_str, "\imaging.nii.gz");
    finalPathMasks = append(PATHCASES, id_str, "\segmentation.nii.gz"); 

    cases{i} = niftiread(finalPathImages);
    masks{i} = niftiread(finalPathMasks);

    %imfuse
end

%To-Dos:
%Funktion zur Visualisierung einer Schicht der Aufnahme (2D)
%sowie Überlargerung der (semitransparenten) Maske über der Aufnahme (2D)

