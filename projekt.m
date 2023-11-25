
%Constants
PATHXCLTABLE = "C:\Users\Tim\Documents\MATLAB\patients_25.xlsx";
PATHCASES = "C:\Users\Tim\Documents\MATLAB\Cases\case_";

%Read Data
data = readtable(PATHXCLTABLE);
patientData = data(strcmp(data.DatensatzVerwenden, "Y"),:);
patientIDs = patientData.CaseID;

cases = cell(1,length(patientIDs));

for i = 1:length(patientIDs)
    id = patientIDs(i,1);
    id_str = string(id);
    id_len = strlength(id_str);
    adding_zeros = 5 - id_len;

    for n = 1:adding_zeros
        id_str = append("0", id_str);
    end

    finalPath = append(PATHCASES, id_str, "\imaging.nii.gz")
    cases{i} = niftiread(finalPath);

    %imfuse
end

test_Case = niftiread("Cases\case_00003\imaging.nii.gz");

