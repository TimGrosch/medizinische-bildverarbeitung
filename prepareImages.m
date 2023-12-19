% Input:
% slice:            to prepare slice of volume
% data:             entry of the patient in the patienttable (one row)

function [img_l, img_r] = prepareImages(slice,data)

    slice = slice(data.gew_hlteAxialeSchichten_z_:data.Var11,:);
    img_l = rescale(slice(:,1:size(slice,2)));
    img_r = rescale(slice(:,size(slice,2)+1:end));

end