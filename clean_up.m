function[img] = clean_up(img,size)

    if nargin == 1
        size = 5;
    end

    %img = bwmorph(img,'fill');
    img = imclose(img,strel('disk',size));
    %img = bwareaopen(img,50);
    

end