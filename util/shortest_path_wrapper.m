function [ spRows , spCols ] = shortest_path_wrapper( bwEdgeImage, startRow, startCol, ... 
                           stopRow, stopCol,  colSpan, weight, imagePadding)
%
% Cuts out a "relevant" region in a black and white image, 
% and finds the shortest path through that region. The shortest path
% is then back-translated to the original coordinates

%
% Input parameters:
%
% bwEdgeImage =  input black and white image (matrix) 
%                containing edges as '1's 
% startRow = start row for shortest path algorithm in original image
% startCol = start column for shortest path in original image
% stopRow =  stop row for shortest path in original image
% stopCol =  stop column for shortest path in original image
% colSpan =  connection column-wise span: every matrix element is connected to 
%            elements in the row below within a "span"
%            [col - colSpan , col - colSpan]              
% weight = cost for jumping one horizontal step to the left or right. 
%           The cost for n steps to left or right is weight^n
% imagePadding = how many columns the the left and right of start/stop
%                column which are added to the cut-out region
%
% Output parameters:
%
% spRows = row numbers of the shortest path in original image
% spCols = column numbers of the shortest path in original image
%
% % Dependencies: shortest_path_bw_image.m
%
 
[numberOfRows ,numberOfCols] = size(bwEdgeImage);

% Cut out relevant part of kymograph
upperCornerX = min(startCol,stopCol) - imagePadding; 
upperCornerX = max(1,upperCornerX);
                 % x-coordinate for upper left edge in cut-out region
upperCornerY = startRow;
                 % y-coordinate for upper left edge in cut-out region
windowWidth = max(startCol,stopCol) + ....
                imagePadding - upperCornerX + 1;
                 % width of cut-out region
windowWidth = min(numberOfCols + 1 - upperCornerX , windowWidth);             
windowHeight = stopRow - startRow + 1;
                 % height of cut-out region                 
% create a "cut-out" image                 
bwEdgeImageCutOut = bwEdgeImage(upperCornerY:upperCornerY + windowHeight - 1 , ...
             upperCornerX:upperCornerX + windowWidth - 1);

% figure out position of keypoints in cut-out image (bwImage)
startXCutOut = startCol - upperCornerX + 1;
startYCutOut = startRow - upperCornerY + 1;
stopXCutOut =  stopCol - upperCornerX + 1;
stopYCutOut = stopRow - upperCornerY + 1;
% find shortest path in cut-out image
[ spRowsDsDNALeftCutOut  spColsDsDNALeftCutOut ] = shortest_path_bw_image( bwEdgeImageCutOut, startYCutOut, ... 
                 startXCutOut, stopYCutOut , stopXCutOut, colSpan, weight);
% translate back to original image coordinates
spRows =  spRowsDsDNALeftCutOut + upperCornerY - 1;
spCols =  spColsDsDNALeftCutOut + upperCornerX - 1;

end
    