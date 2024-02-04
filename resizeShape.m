function resizedMask = resizeShape(originalMask, originalShape, x, y, resizeFactor)
    % Input:
    % originalMask: Binary mask containing the shape
    % resizeFactor: Scaling factor for resizing the shape

    % Ensure resizeFactor is greater than 1 for increasing the size
    if resizeFactor <= 1
        error('Resize factor should be greater than 1 for increasing the size.');
    end

    Shape = imresize(originalShape,resizeFactor);

    area = floor([y-size(Shape,1), y+size(Shape,1);
            x-size(Shape,2), x+size(Shape,2)]/2);

    area = max(area,1);
    area(1,:) = min(area(1,:),size(area,1));
    area(2,:) = min(area(2,:),size(area,2));

    resizedMask = zeros(size(originalMask));
    resizedMask(area(1,1):area(1,2),area(2,1):area(2,2)) = Shape;

    % Display the original and resized masks for verification (optional)
    figure;
    imshow(originalMask), title('Original Mask');
    hold on
    visboundaries(resizedMask,'Color','r')
end