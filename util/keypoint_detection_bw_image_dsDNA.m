function  [atFirstTimeLeftEdgeX  , atFirstTimeLeftEdgeY , ...
           atFirstTimeRightEdgeX , atFirstTimeRightEdgeY ,...
           atCutLeftHairpinLeftEdgeX  , atCutLeftHairpinRightEdgeX   ,...
           atCutRightHairpinLeftEdgeX , atCutRightHairpinRightEdgeX  , ...
           atCutY, ...
           atEndLeftHairpinLeftEdgeX  , atEndLeftHairpinRightEdgeX  , ...
           atEndLeftHairpinY , ...
           atEndRightHairpinLeftEdgeX , atEndRightHairpinRightEdgeX , ...   
           atEndRightHairpinY ] = ...
           keypoint_detection_bw_image_dsDNA(I_mo , ...
               medFilterforBWEdgeDetection, minHairpinTimeSpan)
    

    % Detecting keypoints in a kymograph for an unfolding ds-DNA molecule
    % before unfolding
    %
    % Output:
    %
    % atFirstTimeLeftEdgeX = x-coordinate for left edge in first time frame
    % atFirstTimeLeftEdgeY = y-coordinate (=1) for left edge in first time frame
    %
    % atFirstTimeRightEdgeX = x-coordinate for left edge in first time frame
    % atFirstTimeRightEdgeY = y-coordinate (=1) for left edge in first time frame
    %
    % atCut#X#Hairpin#Y#EdgeX = x-coordinate for left/right edge position 
    %                           of the left/right hairpin at the time of the cut
    % atCutY = y-coordinate at the time of the cut
    %
    % atEnd#X#Hairpin#Y#EdgeX = x-coordinate for left/right edge position (x or y coordinates)
    %                            of the left/right hairpin at the time it
    %                           disappeared
    %
    % atEndLeftHairpinY = y-coordinate at which the left hairpin
    %                     disappeared
    % atEndRightHairpinY = y-coordinate at which the left hairpin
    %                     disappeared    
    %
    % medFilterforBWEdgeDetection = if set to true then a median filter is
    %                          applied to black and white image 
    %
    % minHairpinTimeSpan = if a hairpin persists for a shorter number 
    %                      time frames than minHairpinTimeSpan it does 
    %                      not "count". A larger value gives more robust 
    %                      hairpin detection. NOTE: minHairpinTimeSpan>=0
    %
    % Input:
    %
    % I_mo = black-and white image (multi-Otsued image)
    %

    
    if nargin == 0                        

        I_mo_test=[0 0 1 1 1 1 1 1 0 0;
                   0 0 1 1 1 1 1 1 0 0;
                   0 0 1 1 0 1 1 1 0 0;
                   0 0 1 1 0 0 1 1 0 0; 
                   0 0 0 0 0 0 0 1 0 0; 
                   0 0 0 0 0 0 0 1 0 0; 
                   0 0 0 0 0 0 0 0 0 0; 
                   0 0 0 0 0 0 0 0 0 0; 
                   0 0 0 0 0 0 0 0 0 0];     
               

        I_mo=I_mo_test;


    end

    [M, N]=size(I_mo);       % N = number of columns (x) in image 
                             % M = number of rows (time, t) in image

    % --------------------------------------------
    % --- Filter to deal with too rugged edges ---
    % ---------------------------------------------
    if medFilterforBWEdgeDetection
        I_mo = medfilt2(I_mo); 
    end
    
    % ----------------------
    % -- Detect all edges --
    % ----------------------
    
    LeftEdges = zeros(M,N);
    RightEdges = zeros(M,N);
    numberOfLeftEdges = zeros(M,1);
    numberOfRightEdges = zeros(M,1);
    for rowIdx = 1:M    % loop over rows (time frames)
        x=1:1:N-1;
        Q = I_mo(rowIdx,x)-I_mo(rowIdx,x+1); % Q is non-zero only if an edge is detected  
        leftEdgeTemp = find(Q == -1);  %   detect all left edges
        if ~isempty(leftEdgeTemp)
           leftEdges(rowIdx,1:length(leftEdgeTemp)) = leftEdgeTemp ;
           numberOfLeftEdges(rowIdx) = length(leftEdgeTemp);
        end
        rightEdgeTemp = find(Q == 1);   % detect all right edges
        if ~isempty(rightEdgeTemp)
           rightEdges(rowIdx,1:length(rightEdgeTemp)) = rightEdgeTemp + 1;
           numberOfRightEdges(rowIdx) = length(rightEdgeTemp);
        end
    end
    % check so that the number of left edges = number of right edges
    if sum(numberOfRightEdges - numberOfRightEdges)>0
        disp('number of right edges is not = number of left edges')
    end
    
    
    % ----------------------------------
    % --- Return left and right       ----
    % --- edges at first time frame ----
    % ----------------------------------
    atFirstTimeLeftEdgeX = leftEdges(1,1);
    finalIdx = numberOfRightEdges(1);
    atFirstTimeRightEdgeX = rightEdges(1,finalIdx);
    atFirstTimeLeftEdgeY = 1;
    atFirstTimeRightEdgeY = 1;
    

    % -----------------------------------------
    %  --- Find and return cut coordinates ----
    % -----------------------------------------
    
    rowIdx = 1; 
    foundCut = false;
    while (rowIdx <= M & ~foundCut) 
       rowIdx = rowIdx + 1;
       if (numberOfLeftEdges(rowIdx) == 2 & numberOfLeftEdges(rowIdx-1)==1)
          foundCut = true;
       end
    end
    % Return output
    atCutY = rowIdx ;
    atCutLeftHairpinLeftEdgeX = leftEdges(atCutY,1);
    atCutLeftHairpinRightEdgeX = rightEdges(atCutY,1);
    atCutRightHairpinLeftEdgeX = leftEdges(atCutY,2);
    atCutRightHairpinRightEdgeX = rightEdges(atCutY,2);
    
    
    % ----------------------------------------------
    % --- Identify x and y-coordinates at end   ---- 
    % --- time for the two hairpins and return  ----
    % ----------------------------------------------
    
    
    % In the image, set to zero all rows < cutY. Then find two connected 
    % components (hairpin 1 and hairpin 2). 
    endRow = zeros(1,2);
    endCol = zeros(1,2);
    %
    I_mo(1:atCutY + minHairpinTimeSpan - 1,:)=0;   
    [M N] = size(I_mo);
    cc = bwconncomp(I_mo);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [values idx] = sort(numPixels);
    
     
    % Loop through the largest connected components.
    % Start with largest, then second largest.
    allIdx = zeros(2,2);
    endCols = zeros(2,2);
    meanEndCol = zeros(2,1); 
    for counter = 1:2    
        ccIdx = idx(length(idx)-counter+1);
        cc.PixelIdxList{ccIdx};
        colsCC = floor((cc.PixelIdxList{ccIdx}-1)/M)+1;
        rowsCC = mod((cc.PixelIdxList{ccIdx}-1),M)+1;
        endRow(counter) = max(rowsCC);
        tempIdx = find(rowsCC == endRow(counter));
        endCols(counter,1) = min(colsCC(tempIdx'));
        endCols(counter,2) = max(colsCC(tempIdx'));
        meanEndCol(counter) = mean(endCols(counter,:));
    end
    % "sort" CCs in order to determine "left" and "right" hairpin.  
    [values idx] = sort(meanEndCol);
    
    % Return output 
    atEndLeftHairpinY= endRow(idx(1));
    atEndRightHairpinY = endRow(idx(2));
    atEndLeftHairpinLeftEdgeX = min(endCols(idx(1),:)) - 1;
    atEndLeftHairpinRightEdgeX = max(endCols(idx(1),:)) + 1; 
    atEndRightHairpinLeftEdgeX = min(endCols(idx(2),:)) - 1;
    atEndRightHairpinRightEdgeX = max(endCols(idx(2),:)) + 1; 
   
  
 
end




