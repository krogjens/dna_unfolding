function  [atInputTimeLeftEdgeX , atInputTimeRightEdgeX , ...
           atLastTimeY , atLastTimeLeftEdgeX , atLastTimeRightEdgeX ] = ... 
           keypoint_detection_bw_image_ssDNA(I_mo , tLeftInput , ...
           tRightInput , medFilterforBWEdgeDetection )
    

    % Detecting keypoints in a kymograph for an unfolding DNA molecule
    %
    % Output:
    %
    % atInputTimeLeftEdgeX = x-coordinate for left edge 
    %                        at input time frame, tLeftInput
    %
    % atInputTimeRightEdgeX = x-coordinate for right edge 
    %                         at input time frame, tRightInput
    %
    % atLastTimeY = y-coordinate for last time frame
    % 
    % atLastTimeLeftEdgeX = x-coordinate for left edge in last time frame
    %
    % atLastTimeRightEdgeX = x-coordinate for right edge in last time frame
    %
    % Input:
    %
    % I_mo = black and white image (multi-Otsued image)
    %
    % tLeftInput = input time at which the left edge are to be determined
    %
    % tRightInput = input time at which the right edge are to be determined
    %
    % medFilterforBWEdgeDetection = if set to true then a median filter is
    %                          applied to black and white image 
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
    
    
    % ----------------------------------------------
    % ----- Return edges at input time frame -------
    % ----------------------------------------------
    atInputTimeLeftEdgeX = leftEdges(tLeftInput,1);
    finalIdx = numberOfRightEdges(tRightInput);
    atInputTimeRightEdgeX = rightEdges(tRightInput,finalIdx);
    
    % ----------------------------------------------
    % --- Return coordinates at last time frame ----
    % ----------------------------------------------
    atLastTimeY = M;
    atLastTimeLeftEdgeX = leftEdges(M,1);
    finalIdx = numberOfRightEdges(M);
    atLastTimeRightEdgeX = rightEdges(M,finalIdx);
   
 
end




