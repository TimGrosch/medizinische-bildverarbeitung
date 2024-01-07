function [img_combined] = set2(img, coords, v, yRef, xRef)
%{
DESCRIPTION:
Combine an edge image (described by its coordinates) to an other image:
Set value "v" to coordinates "coords" of image "img"
with initial position yRef and xRef.

INPUT:
img: original image
coords: coordinates of edge image to be combined
v: value to be set (defines the RGB color of the edge, e.g. 1=white
yRef: offset y
xRef: offset x

OUTPUT:
img_combined: combined image
%}

if nargin<4
    yRef=0;
    xRef=0;
end
if xRef<0
    xRef=0;
end
if yRef<0
    yRef=0;
end

s = size(coords,1);
img_combined = img;
for n=1:1:s
    img_combined(yRef+coords(n,1), xRef+coords(n,2)) = v;
end
end