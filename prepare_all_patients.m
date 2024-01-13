% funktion, um alle MRT Volumen und Masken vorzubereiten
% Reinladen und auf die gewünschten Maße zuschneiden
% Inputs:
%    patient_table: Tabelle der Patienten, die reingeladen werden sollen
%    path_dir:      Pfad, des Ordners, in dem alle Fälle gespeichert sind
%    path_tar:      Pfad auf dem die Matrizen gespeichert werden sollen

function [] = prepare_all_patients(patient_table,path_dir,path_tar)

    for i = 1:size(patient_table,1)
        cd(path_dir)
        patient = patient_table(i,:);

        adding_zeros = 5 - numel(num2str(patient.CaseID));

        id_str = string(patient.CaseID);
        for n = 1:adding_zeros
            id_str = append("0", id_str);
        end
        %path_dir = append(path_dir,'/case_' ,id_str, "/imaging.nii.gz");

        V = niftiread(append('case_',id_str,'/imaging.nii.gz'));
        mask = niftiread(append('case_',id_str,'/segmentation.nii.gz'));

        V = interpolate_nifti_z(V,patient);
        mask = interpolate_nifti_z(mask,patient);
        
        cd(path_tar)

        save(id_str,"V","mask");

    end

end