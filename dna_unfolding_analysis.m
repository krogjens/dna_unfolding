function dna_unfolding_analysis(kymo,pixelSize,frameRate)

if nargin == 0 % If no arguments are given perform the default example analysis
   % Target kymograph folder and perform Otsu segmentation to obtain coordinates
   kymo = 'p3 and l3-ZVI Export-03_molecule_1_kymograph.tif'; % Path to kymograph
   expName = kymo(1:end-14); % Snap off "_kymograph.tif" to make a folder for analysis files
   mkdir(expName) % Folder where analysis files and results will be placed

   pixelSize = 1.59e-7; % Manually add the pixelsize (meter)
   frameRate = 9; % Manually add the framerate for the used camera (Hz)
end
addpath('util')
addpath('ns')

% Trace out edges in kymographs
FoundEdges = dna_edge_tracing_skel(kymo,expName);

if FoundEdges,
    fprintf('Edges detected in %s. \nStarting analysis \n',kymo);
    
    % Collect coordinates into x_l, x_r and L data file along resolution and framerate
    DataFileName = dna_coord_collector(expName,pixelSize,frameRate);

    % Run hp_skeleton and create bayesout.mat and results.txt files
    dna_bayes_skel(expName,DataFileName)
    fprintf('Bayesian analysis completed. \n')
    
    % Plot edges on top of kymograph
    figure(1)
    dna_unfolding_plotter(kymo,expName);
    % Create figure showing likelihood function
    figure(2)
    dna_likelihood_plotter(kymo,expName)
    
else
    fprintf('Was not able to find edges in kymograph: %s \n',kymo)
end
