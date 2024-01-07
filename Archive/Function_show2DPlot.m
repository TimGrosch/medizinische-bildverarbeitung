function [img] = show2DPlot(V, mask, chosen_layer, y_min, y_max)
    figure;
    img = squeeze(V(y_min:y_max,chosen_layer,:));
    imshow(img);
    mask = squeeze(mask(y_min:y_max,chosen_layer,:));
    fused = imfuse(img, mask);
    diffused = imdiffusefilt(fused);
    %imshow(diffused);
end