function [score,y,x,acc] = generalized_hough_transform(target, reference, refY, refX)
%{
DESCRIPTION:
modified from: http://mbcoder.com/generalized-hough-transform/
Find reference image (edge image) in target image (edge image)
using generalized hough transform (GHT).

INPUT:
target: target image as binary edge image (edges = 1, rest = 0)
reference: reference image as binary edge image (edges = 1, rest = 0)
refY, refX: reference point (y,x coordinates) of the reference image,
e.g. middle point of the reference or top-left corner etc.

OUTPUT
score: maximum score of the best match
y, x: coordinates of the reference point of the best match
acc: accumulator (hough space)

%}

%--- Get the reference edge points ----------------------------------------
% find all y,x cordinates equal 1 in reference image (=edges)
[y,x] = find(reference>0); 
maxPoints = size(x,1); % number of points in the reference image

if (maxPoints<1)
    disp('Error: Cannot determine shape of reference (no data points equal 1)');
    quit();
end

% --- Create gradient map of reference ------------------------------------
gradient_ref = gradient_direction(reference);

% --- Create R-Table ------------------------------------------------------
maxAngles = 180; % devide the angel space to maxAngles uniformed space bins
binCounter = zeros(maxAngles); % counter for the amount of edge points associated with each angle gradient
Rtable = zeros(maxAngles, maxPoints, 2);
for k=1:1:maxPoints
    bin = round((gradient_ref(y(k),x(k))/pi)*(maxAngles-1))+1; % transform continuous gradient angles to discrete angle bins
    binCounter(bin) = binCounter(bin) + 1; % increase binCounter at that position
    Rtable(bin, binCounter(bin), 1) = y(k) - refY; % offset y directiom
    Rtable(bin, binCounter(bin), 2) = x(k) - refX; % offset x direction
end

%--- Get the target edge points--------------------------------------------
% find all y,x cordinates equal 1 in target image (=edges)
[y, x] = find(target>0);
maxPoints_tar = size(x,1);
if (maxPoints_tar<1)
    disp('Error: Cannot determine edges of target (no data points equal 1)');
    quit();
end

% --- Create gradient map of target ---------------------------------------
gradient_tar = gradient_direction(target);
size_tar = size(target);

% --- Create and populate hough space (accumulator) -----------------------
houghspace = zeros(size_tar);

for k = 1:1:maxPoints_tar
    bin = round((gradient_tar(y(k), x(k))/pi)*(maxAngles-1))+1; % transform continuous gradient angles to discrete angle bins
    for j = 1:1:binCounter(bin)
        ty = y(k) - Rtable(bin, j, 1);
        tx = x(k) - Rtable(bin, j, 2);
        if (ty>0) && (ty<size_tar(1)) && (tx>0) && (tx<size_tar(2))
            houghspace(ty, tx) =  houghspace(ty, tx)+1; % increase hough space at that position
        end
    end
end

% --- Find best match in hough space --------------------------------------
acc = houghspace;
acc = acc./sqrt(sum(sum(reference))); % normalize acc by the size of the reference to avoid bias towards large reference
score = max(max(acc)); % use maximum of accumulator as score
[y, x] = find(acc == max(max(acc)), 1, 'first'); % find maximum position

%{
Optional: Find a better metric (score) to describe the quality of the best
match that can be used in "find_object.m" to determine if a match was
found at all. Currently max(acc) is used for that purpose.
%}

end

