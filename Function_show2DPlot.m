function show2DPlot(chosen_case, chosen_layer, cases_cell, masks_cell, y_min, y_max)
    figure;
    img = squeeze(cases_cell{chosen_case}(y_min:y_max,chosen_layer,:));
    mask = squeeze(masks_cell{chosen_case}(y_min:y_max,chosen_layer,:));
    fused = imfuse(img, mask);
    diffused = imdiffusefilt(fused);
    imshow(diffused);
end