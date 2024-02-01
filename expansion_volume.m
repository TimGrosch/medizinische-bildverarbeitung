function[new_mask] = expansion_volume(seg,starting_slice,V,side)

if strcmp(side,'left')
    V = V(:,:,257:end);
else
    V = V(:,:,1:256);
end

% seg = kidney_seg;
% slice = example_coronal_layer;
% V = V;

% 229 - 326

new_mask = zeros(size(V));
new_seg = seg;
current_slice = starting_slice;

while true
    current_slice = current_slice-1;
    inbetween = activecontour(squeeze(V(:,current_slice,:)),new_seg,'Chan-vese');

    % cutoff metric TODO
    if nnz(inbetween) > 1.5*nnz(new_seg)
        fprintf('Kidney stops at slice: %d\n',current_slice)
        break;
    else
        new_seg = inbetween;
    end

    % deleting every object except for the largest
    max_size = 0;
    list_objects = bwconncomp(new_seg);
    list_objects = list_objects.PixelIdxList;
    for j = 1:length(list_objects)
        arr = list_objects(j);
        curr = length(arr{1});
        if curr > max_size
            max_size = curr;
        end
    end
    new_seg = bwareaopen(new_seg, max_size-1);

    imshow(new_seg,[])
    new_mask(:,current_slice,:) = new_seg;
end

new_seg = seg;
current_slice = starting_slice;

while true
    current_slice = current_slice+1;
    inbetween = activecontour(squeeze(V(:,current_slice,:)),new_seg,'Chan-vese');

    % cutoff metric TODO
    if nnz(inbetween) > 1.5*nnz(new_seg)
        fprintf('Kidney stops at slice: %d\n',current_slice)
        break;
    else
        new_seg = inbetween;
    end

    % delete every object except for the largest in the mask
    max_size = 0;
    list_objects = bwconncomp(new_seg);
    list_objects = list_objects.PixelIdxList;
    for j = 1:length(list_objects)
        arr = list_objects(j);
        curr = length(arr{1});
        if curr > max_size
            max_size = curr;
        end
    end
    new_seg = bwareaopen(new_seg, max_size-1);


    imshow(new_seg,[])
    new_mask(:,current_slice,:) = new_seg;
end

end