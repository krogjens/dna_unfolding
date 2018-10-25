function FoundEdges = hairpin_tracing(inputFileName, inputDirName , ... ,
    outputDirName, displayImages, saveImages , saveTrajectories)

    % Detects and traces out edges in a kymograph with 
    % circular-to-linear unfolding event
    % 
    % inputDirName = name of directory with input files
    % inputFilename = name of iput tif-file with kymograph
    % outputDirName = name of output directory
    % displayImages = true if one wants to display images
    % saveImages = true if one want to save the output images
    % saveTrajectories = true if one want to save the trajectories
    %                    as text-files
    %


    % CHOSE SMOOTHING METHOD

    method = 5; % For the published dna unfolding analysis we use method 5 which has 
                % no inherent smoothing, since we want all the noise dynamics to be
                % present for the analysis. Be mindful that smoothing may reduce variance
                % and thus induce artificially large friction and force estimates.   
    
    % -----------------------------
    % ---- "Tweak" parameters -----
    % -----------------------------
    sigmaSmooth = 4;         % width of smoothing gaussian (units of pixels)
                             % used for pre-processing of images 
    weight = 3;            % cost for jumping to a white node in 
                             % shortest path algorithm (the cost is minimized)
                             % is weight^n, where n is the horizontal distance
                             % between white pixels in neighbouring rows
                             % note: choose w>=1
    colSpan = 17;           % column-wise span: every white element in 
                             % multi-Otsued black and white image 
                             % is connected to white  
                             % elements in the row below within a "span"
                             %  [col - colSpan , col - colSpan] 
                             % when doing shortest path algorithm to find
                             % edge trajactories.
    imagePadding = 20;       % when "cutting out" a relevant part of the image,
                             % this is the extra amount of padding added to the
                             % left and right
    minHairpinTimeSpan = 0;  % if a hairpin persists for a shorter number 
                             % time frames than minHairpinTimeSpan it does 
                             % not "count". A larger value gives more robust 
                             % hairpin detection. NOTE: minHairpinTimeSpan>=0
    medFilterforBWEdgeDetection = false;
                             % if = true then a median filter is applied
                             % before the derivative filter is applied.
                             % This avoids "hook" effect, i.e., that
                             % more than 2 edges are detected before cut
                             % occured. Must be used if no Pre-processing
                             % is used (case 5. below)
    
    
    % ----------------------------------------
    % ---- Load the image from .tif file. ----
    % ----------------------------------------
    
    %fullInputFileName = strcat(inputDirName,'/',inputFileName);
    fullInputFileName = inputFileName;
    imgArr = GetImgArrFromFile(fullInputFileName);    
    imgArr = imgArr ./ max(imgArr(:));
    origImgArr = imgArr;
    
    
    % --------------------------------------------------------
    % ---- Pre-processing (filtering) for noise reduction ----
    % --------------------------------------------------------

if method == 1;
      % --------------------------------------
      % ---            1.              -------
      % --- OLD method: 2D Gaussian filter ---
      % --------------------------------------
     sigmaSmoothOld = 8;
     f = fspecial('gaussian', [7 3], sigmaSmoothOld);
     imgArrSmooth = conv2(origImgArr, f, 'valid');

elseif method == 2

      % ------------------------------------------
      % ---            2.                  -------
      % --- 1D Gaussian filter + Median filter ---
      % --------------------------------------
    
     % [optimal value is sigma_smooth = sigma_psf = 2 pixels, see 
     % P.C. Torche, Noise reduction in single time frame optical DNA maps,PLOS
     % One (2017), Suppl Fig. S5.] 
     f = fspecial('gaussian', [round(4*sigmaSmooth) 1], sigmaSmooth);
     imgArrSmooth = conv2(origImgArr, f, 'valid');
     imgArrSmooth = medfilt2(imgArrSmooth); 

elseif method == 3
      % ---------------------------------------------
      % ---            3.                      ------
      % --- 1D Window-Sinc filter + Median filter ---
      % ---------------------------------------------
     numberOfRows = size(origImgArr,1);
     for rowIdx = 1:numberOfRows
        sigmaPSF = 2;      % width of optical point spread function (units of pixels)
        fc = 1/(pi*sigmaPSF); % optimal value, see P.C. Torche article above
        imgArrSmooth(rowIdx,:) = window_sinc_filter(origImgArr(rowIdx,:),fc);
     end;   
     imgArrSmooth = medfilt2(imgArrSmooth); 

elseif method == 4
     % ------------------------------------------
     % ----           4.                      ---
     % ---         Only Median filter         ---
     % ------------------------------------------    
     imgArrSmooth = medfilt2(origImgArr); 
  
else
     % ------------------------------------------
     % ----           5.                      ---
     % ---         No filter                  ---
     % ------------------------------------------ 
%      medFilterforBWEdgeDetection = false;
      medFilterforBWEdgeDetection = true;
      imgArrSmooth = origImgArr;
     
end     
     
    % -------------------------------------------------
    % ---- Thresholding via the multi-Otsu method -----
    % -------------------------------------------------
    
    % --- Multi-Otsu with three levels -----
    % --- (background,ss-DNA,ds-DNA)   -----
    
    thresh = multithresh(imgArrSmooth,2);
    seg_image = imquantize(imgArrSmooth, thresh);
    % region with double-stranded (intact) DNA
    seg_dsDNA=seg_image;
    seg_dsDNA(seg_dsDNA~=3) = 0;   
    seg_dsDNA(seg_dsDNA==3) = 1;
    % region with single-stranded DNA
    seg_ssDNA=seg_image;
    seg_ssDNA(seg_ssDNA~=2) = 0;   
    seg_ssDNA(seg_ssDNA==2) = 1;
    
    % ds-DNA: Retain only largest connected component,
    % and fill in "holes".
    CC_dsDNA = bwconncomp(seg_dsDNA);
    numPixels = cellfun(@numel,CC_dsDNA.PixelIdxList);
    [~,idx] = max(numPixels);
    seg_dsDNA(CC_dsDNA.PixelIdxList{idx}) = 2;
    seg_dsDNA(seg_dsDNA ~= 2) = 0;
    seg_dsDNA(seg_dsDNA == 2) = 1;
    seg_dsDNA = imfill(seg_dsDNA,'holes');   
        
    % ss-DNA region: retain only largest connected component,
    % and fill in "holes" 
    CC_ssDNA = bwconncomp(seg_ssDNA);
    numPixels = cellfun(@numel,CC_ssDNA.PixelIdxList);
    [~,idx] = max(numPixels);
    seg_ssDNA(CC_ssDNA.PixelIdxList{idx}) = 2;
    seg_ssDNA(seg_ssDNA ~= 2) = 0;
    seg_ssDNA(seg_ssDNA == 2) = 1;
    seg_ssDNA = imfill(seg_ssDNA,'holes');
    
        
    % --- Apply derivative filter to b-w ----
    % --- images to identify edges  ----

    % --- ds-DNA ---
    
    % Identify left edges
    f = [1 -1];
    if medFilterforBWEdgeDetection
        seg_dsDNA_ForEdgeDetection = medfilt2(seg_dsDNA); 
    else
        seg_dsDNA_ForEdgeDetection = seg_dsDNA;
    end
    % Note: whatever filter is used here must also be used in
    % keypoint_detection_bw_image_dsDNA.m 
    imgLeftEdges_dsDNA = conv2(seg_dsDNA_ForEdgeDetection, f, 'valid');
    imgLeftEdges_dsDNA(imgLeftEdges_dsDNA < 0) = 0;
                        % set all '-1' (right edges) to 0 

    % Identify right edges
    f = [-1 1];
    if medFilterforBWEdgeDetection
        seg_dsDNA_ForEdgeDetection = medfilt2(seg_dsDNA);
    else
        seg_dsDNA_ForEdgeDetection = seg_dsDNA;
    end
    % Note: whatever filter is used here must also be used in
    % keypoint_detection_bw_image_dsDNA.m 
    imgRightEdges_dsDNA = conv2(seg_dsDNA_ForEdgeDetection , f, 'valid');
    imgRightEdges_dsDNA(imgRightEdges_dsDNA < 0) = 0;
                        % set all '-1' (right edges) to 0 
    imgRightEdges_dsDNA = circshift(imgRightEdges_dsDNA,[0 1]);
                        % shift by one column to retain symmetry
                        % between left and right hairpins
   
    % --- ss-DNA ---
    
    % Identify left edges
    f = [1 -1];
    if medFilterforBWEdgeDetection
        seg_ssDNA_ForEdgeDetection = medfilt2(seg_ssDNA); 
    else
        seg_ssDNA_ForEdgeDetection = seg_ssDNA;
    end
    % Note: whatever filter is used here must also be used in
    % keypoint_detection_bw_image_ssDNA.m    
    imgLeftEdges_ssDNA = conv2(seg_ssDNA_ForEdgeDetection, f, 'valid');
    imgLeftEdges_ssDNA(imgLeftEdges_ssDNA < 0) = 0;
                        % set all '-1' (right edges) to 0 

    % Identify right edges
    f = [-1 1];
    if medFilterforBWEdgeDetection
        seg_ssDNA_ForEdgeDetection = medfilt2(seg_ssDNA);
    else
        seg_ssDNA_ForEdgeDetection = seg_ssDNA;
    end
     % Note: whatever filter is used here must also be used in
    % keypoint_detection_bw_image_ssDNA.m    
    imgRightEdges_ssDNA = conv2(seg_ssDNA_ForEdgeDetection, f, 'valid');
    imgRightEdges_ssDNA(imgRightEdges_ssDNA < 0) = 0;
                        % set all '-1' (right edges) to 0 
    imgRightEdges_ssDNA = circshift(imgRightEdges_ssDNA,[0 1]);
                        % shift by one column to retain symmetry
                        % between left and right hairpins
   
    
                        
    % -------------------------------------------------------------------
    % ----- Detect keypoints in the black-and-white dsDNA kymograph -----
    % -------------------------------------------------------------------

    
   [atFirstTimeLeftEdgeX  , atFirstTimeLeftEdgeY , atFirstTimeRightEdgeX , ...
   atFirstTimeRightEdgeY , atCutLeftHairpinLeftEdgeX  , ... 
   atCutLeftHairpinRightEdgeX   , atCutRightHairpinLeftEdgeX , ... 
   atCutRightHairpinRightEdgeX  , atCutY, atEndLeftHairpinLeftEdgeX , ... 
   atEndLeftHairpinRightEdgeX  , atEndLeftHairpinY , ...
   atEndRightHairpinLeftEdgeX , atEndRightHairpinRightEdgeX , ...
   atEndRightHairpinY ...
   ] = keypoint_detection_bw_image_dsDNA(seg_dsDNA , ...
               medFilterforBWEdgeDetection , minHairpinTimeSpan);  

                     
    % -------------------------------------------------------------------
    % ----- Detect keypoints in the black-and-white ssDNA kymograph -----
    % -------------------------------------------------------------------

    [atInputTimeLeftEdgeX , atInputTimeRightEdgeX , ...
           atLastTimeY , atLastTimeLeftEdgeX , atLastTimeRightEdgeX ] = ... 
           keypoint_detection_bw_image_ssDNA(seg_ssDNA , ... 
              atEndLeftHairpinY , atEndRightHairpinY , ...
               medFilterforBWEdgeDetection );
    
    
    % ----------------------------------------------------
    % -------               ds-DNA:                -------
    % ------- Link the six edges into trajectories -------
    % ------- (shortest path algorithm)            -------
    % ------- This sensitive to "good" segmentation ------
    % ------- so catch error if they occur         -------
    % ----------------------------------------------------
     
    errorInShortestPath = false;
    try 
        
        % ---  1. Intact (ds-DNA) DNA, left edge trajactory ---

        [ spRowsDsDNALeft , spColsDsDNALeft ] = ... 
          shortest_path_wrapper(imgLeftEdges_dsDNA, atFirstTimeLeftEdgeY , ... 
          atFirstTimeLeftEdgeX  , atCutY,  atCutLeftHairpinLeftEdgeX, ... 
          colSpan, weight, imagePadding);


        % ---  2. Intact (ds-DNA) DNA, right edge trajactory ---

         [ spRowsDsDNARight , spColsDsDNARight  ] = ...
           shortest_path_wrapper(imgRightEdges_dsDNA, atFirstTimeRightEdgeY , ... 
           atFirstTimeRightEdgeX  , atCutY,  atCutRightHairpinRightEdgeX, ... 
           colSpan, weight, imagePadding);


        % --- 3. Left edge of left-most hairpin  -------

        [ spRowsLeftHairpinLeftEdge, spColsLeftHairpinLeftEdge] = ... 
          shortest_path_wrapper(imgLeftEdges_dsDNA, atCutY , ... 
          atCutLeftHairpinLeftEdgeX , atEndLeftHairpinY, atEndLeftHairpinLeftEdgeX ,   ... 
          colSpan, weight, imagePadding);


         % ---  4. Right edge of left-most hairpin  -------

         [ spRowsLeftHairpinRightEdge, spColsLeftHairpinRightEdge] = ...
           shortest_path_wrapper(imgRightEdges_dsDNA, atCutY , ... 
           atCutLeftHairpinRightEdgeX , atEndLeftHairpinY, atEndLeftHairpinRightEdgeX ,   ... 
           colSpan, weight, imagePadding);


         % --- 5. Left edge of right-most hairpin  -------
         [ spRowsRightHairpinLeftEdge, spColsRightHairpinLeftEdge] = ... 
           shortest_path_wrapper(imgLeftEdges_dsDNA, atCutY , ... 
           atCutRightHairpinLeftEdgeX , atEndRightHairpinY, atEndRightHairpinLeftEdgeX ,   ... 
           colSpan, weight, imagePadding);


         % --- 6. Right edge of right-most hairpin  -------
        [ spRowsRightHairpinRightEdge, spColsRightHairpinRightEdge] = ...
          shortest_path_wrapper(imgRightEdges_dsDNA, atCutY , ... 
          atCutRightHairpinRightEdgeX , atEndRightHairpinY, atEndRightHairpinRightEdgeX ,   ... 
          colSpan, weight, imagePadding);
     
      % Indicate whether hairpins were detected with length > 2 (as needed
      % by bayesian analysis)
      if length(spRowsLeftHairpinLeftEdge) > 2 && length(spRowsRightHairpinLeftEdge) > 2
          FoundEdges = true;
      else
          FoundEdges = false;
      end


        % ----------------------------------------------------
        % -------               ss-DNA:                -------
        % ------- Link the two edges into trajectories -------
        % ------- (shortest path algorithm)            -------
        % ----------------------------------------------------

        % ---  1. ss-DNA, left edge trajactory ---

        [ spRowsSsDNALeft , spColsSsDNALeft ] = ... 
          shortest_path_wrapper(imgLeftEdges_ssDNA, atEndLeftHairpinY , ... 
          atInputTimeLeftEdgeX  , atLastTimeY , atLastTimeLeftEdgeX , ... 
          colSpan, weight, imagePadding);


        % ---  2. ss-DNA, right edge trajactory ---

        [ spRowsSsDNARight , spColsSsDNARight ] = ... 
          shortest_path_wrapper(imgRightEdges_ssDNA, atEndRightHairpinY , ... 
          atInputTimeRightEdgeX  , atLastTimeY , atLastTimeRightEdgeX , ... 
          colSpan, weight, imagePadding);

    catch MExp
        
        disp(strcat('Error in shortest path algorithm in file: ',inputFileName))
        disp(MExp.message)
        disp(char(10)) 
        errorInShortestPath = true;
        
    end 
  
    % -----------------------------------
    % ---------------- Plot -------------
    % -----------------------------------
      
    % --- Side-by-side montage --- 
% figure(42)
%         subimage(origImgArr), xlabel('raw')  
%             hold on, set(gca,'xtick',[]), set(gca,'ytick',[]) 
%             set(gca,'FontSize',10)
%              % 1. Intact DNA, left edge trajectory
%              plot(spColsDsDNALeft,spRowsDsDNALeft,'g.')
% 
%             % 2. Intact DNA, left edge trajectory
%             plot(spColsDsDNARight,spRowsDsDNARight,'g.')
% 
%             % 3. Left hairpin, left edge trajectory
%             plot(spColsLeftHairpinLeftEdge,spRowsLeftHairpinLeftEdge,'r.')
% 
%             % 4. Left hairpin, right edge trajectory
%             plot(spColsLeftHairpinRightEdge,spRowsLeftHairpinRightEdge,'r.')
% 
%             % 5. Left hairpin, left edge trajectory
%             plot(spColsRightHairpinLeftEdge,spRowsRightHairpinLeftEdge,'m.')
% 
%             % 6. Left hairpin, left edge trajectory
%             plot(spColsRightHairpinRightEdge,spRowsRightHairpinRightEdge,'m.')
%             % --- ss-DNA ---
% 
%             % 1. Intact DNA, left edge trajectory
%              plot(spColsSsDNALeft,spRowsSsDNALeft,'b.')
% 
%             % 2. Intact DNA, left edge trajectory
%             plot(spColsSsDNARight,spRowsSsDNARight,'b.')


    
   
    if displayImages
        
        figure
        h = suptitle(inputFileName);
        set(h,'Interpreter','none')
        set(h,'Fontsize',10)

        subplot(2,3,1), subimage(origImgArr), xlabel('raw')  
            hold on, set(gca,'xtick',[]), set(gca,'ytick',[]) 
            set(gca,'FontSize',10)

        subplot(2,3,2), subimage(imgArrSmooth)  
           xlabel('smoothed image') , hold on, 
           set(gca,'xtick',[]), set(gca,'ytick',[])
           set(gca,'FontSize',10)

        subplot(2,3,3), subimage(seg_dsDNA), xlabel('ds-DNA')   
            hold on, set(gca,'xtick',[]), set(gca,'ytick',[])
            set(gca,'FontSize',10)         

        subplot(2,3,4), subimage(seg_ssDNA),  xlabel('ss-DNA') , 
           hold on, set(gca,'xtick',[]), set(gca,'ytick',[])  
           set(gca,'FontSize',10)

        subplot(2,3,5), subimage(imgLeftEdges_dsDNA +imgRightEdges_dsDNA )  
           xlabel('derivative ds-DNA') , hold on, 
           set(gca,'xtick',[]), set(gca,'ytick',[])

         subplot(2,3,6), subimage(imgLeftEdges_ssDNA +imgRightEdges_ssDNA )  
           xlabel('derivative ss-DNA') , hold on, 
           set(gca,'xtick',[]), set(gca,'ytick',[])

        colormap(gray)  
        hold on;


        % --- Mark "keypoints" in smoothed image ---

        subplot(2,3,2)
        plot(atFirstTimeLeftEdgeX,atFirstTimeLeftEdgeY,'go','linewidth',2);
        plot(atFirstTimeRightEdgeX ,atFirstTimeRightEdgeY , 'go', 'linewidth', 2);
        %
        plot(atCutLeftHairpinLeftEdgeX, atCutY, 'ro', 'linewidth', 2);
        plot(atCutLeftHairpinRightEdgeX, atCutY, 'ro', 'linewidth', 2);
        plot(atCutRightHairpinLeftEdgeX, atCutY, 'mo', 'linewidth', 2);
        plot(atCutRightHairpinRightEdgeX, atCutY, 'mo', 'linewidth', 2);
        %
        plot(atEndLeftHairpinLeftEdgeX, atEndLeftHairpinY, 'ro', 'linewidth', 2);
        plot(atEndLeftHairpinRightEdgeX, atEndLeftHairpinY, 'ro', 'linewidth', 2);
        plot(atEndRightHairpinLeftEdgeX, atEndRightHairpinY, 'mo', 'linewidth', 2);
        plot(atEndRightHairpinRightEdgeX, atEndRightHairpinY, 'mo', 'linewidth', 2);



        % --- Mark "keypoints" smoothed image ---

        subplot(2,3,2)
        plot(atInputTimeLeftEdgeX, atEndLeftHairpinY, 'bo', 'linewidth', 2);
        plot(atLastTimeLeftEdgeX, atLastTimeY, 'co', 'linewidth', 2);
        %
        plot(atInputTimeRightEdgeX, atEndRightHairpinY, 'bo', 'linewidth', 2);
        plot(atLastTimeRightEdgeX, atLastTimeY, 'co', 'linewidth', 2);


        % --- Mark edge trajectories if they were corrected identified ---
        
        if ~errorInShortestPath 
            % --- ds-DNA ---

             % 1. Intact DNA, left edge trajectory
             subplot(2,3,1)
             plot(spColsDsDNALeft,spRowsDsDNALeft,'g.')

            % 2. Intact DNA, left edge trajectory
            subplot(2,3,1)
            plot(spColsDsDNARight,spRowsDsDNARight,'g.')

            % 3. Left hairpin, left edge trajectory
            subplot(2,3,1)
            plot(spColsLeftHairpinLeftEdge,spRowsLeftHairpinLeftEdge,'r.')

            % 4. Left hairpin, right edge trajectory
            subplot(2,3,1)
            plot(spColsLeftHairpinRightEdge,spRowsLeftHairpinRightEdge,'r.')

            % 5. Left hairpin, left edge trajectory
            subplot(2,3,1)
            plot(spColsRightHairpinLeftEdge,spRowsRightHairpinLeftEdge,'m.')

            % 6. Left hairpin, left edge trajectory
            subplot(2,3,1)
            plot(spColsRightHairpinRightEdge,spRowsRightHairpinRightEdge,'m.')


            % --- ss-DNA ---

            % 1. Intact DNA, left edge trajectory
             subplot(2,3,1)
             plot(spColsSsDNALeft,spRowsSsDNALeft,'b.')

            % 2. Intact DNA, left edge trajectory
            subplot(2,3,1)
            plot(spColsSsDNARight,spRowsSsDNARight,'b.')
        end

     end
    
    
     % ----------------------------------
     % --------------- Save -------------
     % ----------------------------------
     
    % --- Save images ---
    
    if saveImages
        
        fullOutputImageFileName = strcat(outputDirName,'/', ...
            inputFileName(1:end-4), '_images.png');
        print(gcf, '-dpng', fullOutputImageFileName);
        disp(char(10))
        disp(['images saved to:' 9 9 fullOutputImageFileName])
        
    end
    
   
    % --- Save the edge positions to text-files ---
    if saveTrajectories && ~errorInShortestPath
        
        % Intact DNA (kymograph before the cut)

        fullOutputFileName = strcat(outputDirName,'/',...
            inputFileName(1:end-4), '_intactDsDNA.txt');
        fileID = fopen(fullOutputFileName,'w');
        fprintf(fileID,'%22s %25s %20s\n','time','dsDNA left-edge','dsDNA right-edge');
        for counter = 1:length(spRowsDsDNALeft)
            fprintf(fileID,'%20.0f %20.0f %20.0f\n',...
              spRowsDsDNALeft(counter),spColsDsDNALeft(counter), ... 
              spColsDsDNARight(counter));
        end
        fclose(fileID);   
        disp(['Data saved to:' 9 9 fullOutputFileName])


         % Left hairpin 

        fullOutputFileName = strcat(outputDirName,'/',... 
            inputFileName(1:end-4), '_leftHairpin.txt');
        fileID = fopen(fullOutputFileName,'w');
        fprintf(fileID,'%20s %27s %20s\n','time','Hairpin left-edge','Hairpin right-edge');
        for counter = 1:length(spRowsLeftHairpinLeftEdge)
            fprintf(fileID,'%20.0f %20.0f %20.0f\n',...
              spRowsLeftHairpinLeftEdge(counter), spColsLeftHairpinLeftEdge(counter), ...
              spColsLeftHairpinRightEdge(counter));
        end
        fclose(fileID);   
        disp(['Data saved to:' 9 9 fullOutputFileName])

        % Right hairpin 

        fullOutputFileName = strcat(outputDirName,'/',...
            inputFileName(1:end-4), '_rightHairpin.txt');
        fileID = fopen(fullOutputFileName,'w');
        fprintf(fileID,'%20s %27s %20s\n','time','Hairpin left-edge','Hairpin right-edge');
        for counter = 1:length(spRowsRightHairpinLeftEdge)
            fprintf(fileID,'%20.0f %20.0f %20.0f\n',...
              spRowsRightHairpinLeftEdge(counter), spColsRightHairpinLeftEdge(counter), ...
              spColsRightHairpinRightEdge(counter));
        end
        fclose(fileID);   
        disp(['Data saved to:' 9 9 fullOutputFileName])

        % ss-DNA regions (after unfolding of first/second hairpin)

        % Left edge

        fullOutputFileName = strcat(outputDirName,'/',...
            inputFileName(1:end-4), '_ssDNALeftEdge.txt');
        fileID = fopen(fullOutputFileName,'w');
        fprintf(fileID,'%21s %25s\n','time','ssDNA left-edge');
        for counter = 1:length(spRowsSsDNALeft)
            fprintf(fileID,'%20.0f %20.0f\n',...
              spRowsSsDNALeft(counter),spColsSsDNALeft(counter));
        end
        fclose(fileID);   
        disp(['Data saved to:' 9 9 fullOutputFileName])

        % Right edge

        fullOutputFileName = strcat(outputDirName,'/',...
            inputFileName(1:end-4), '_ssDNARightEdge.txt');
        fileID = fopen(fullOutputFileName,'w');
        fprintf(fileID,'%21s %25s\n','time','ssDNA right-edge');
        for counter = 1:length(spRowsSsDNARight)
            fprintf(fileID,'%20.0f %20.0f\n',...
              spRowsSsDNARight(counter),spColsSsDNARight(counter));
        end
        fclose(fileID);   
        disp(['Data saved to:' 9 9 fullOutputFileName])

%        printspecs(inputDirname,outputDirname,inputFilename(1:end-4));
    end
    
    
    
end

%============================ Subfunctions ===============================%

%=========================================================================%
function [ imgArr ] = GetImgArrFromFile(fullInputFileName)

    tiffInfo = imfinfo(fullInputFileName);   % Get the TIF file information
    no_frame = numel(tiffInfo);     % Get the number of images in the file
    traceCell = cell(no_frame,1);   % Preallocate the cell array
    for frame = 1:no_frame
      traceCell{frame} = imread(fullInputFileName,'Index',frame,'Info',tiffInfo);
    end
    data = traceCell{1,1};
    
    imgArr = zeros(size(data));
    
    for i = 1:size(data,1)
        imgArr(i,:) = data(i,:);
    end
    
end

