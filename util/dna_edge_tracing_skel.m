function FoundEdges = dna_edge_tracing_skel(kymo,expdir)
% Analyze all files
outputDirName_0p05X = expdir;   % must exist
displayImages = false;
saveImages = false;
saveTrajectories = true;       

try   
   FoundEdges = hairpin_tracing(kymo,expdir,outputDirName_0p05X, displayImages, saveImages , saveTrajectories );
catch MExp
   FoundEdges = false;
   disp('Error which is not shortest path related.')
   disp(MExp.message)
   disp(char(10)) 
end 
