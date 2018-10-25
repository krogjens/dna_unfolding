function [ spRows , spCols ] = shortest_path_bw_image( bwImage, startRow, startCol, ... 
                           stopRow, stopCol, colSpan, weight)
%
% Find the shortest path through a black and white image, 
% where all 1s are connected across rows.
%
% Input parameters:
%
% bwImage =  input black and white image (matrix)
% startRow = start row for shortest path algorithm
% startCol = start column for shortest path 
% stopRow =  stop row for shortest path
% stopCol =  stop column for shortest path 
% colSpan =  connection column-wise span: every matrix element is connected to 
%            elements in the row below within a "span"
%            [col - colSpan , col - colSpan]              
% weight = cost for jumping one horizontal step to the left or right. 
%           The cost for n steps to left or right is weight^n
%
% Output parameters:
%
% spRows = row numbers of the shortest path 
% spCols = column numbers of the shortest path
%
% NOTE: this function has not been optimized for speed.
%
% Initialize
numberOfRows = size(bwImage,1);
numberOfCols = size(bwImage,2);               
numberOfNodes = length(find(bwImage == 1 ));

% Make "translation tables" between
% node index to (rows,columns)
nodeCol = zeros(numberOfNodes,1);  
nodeRow = zeros(numberOfNodes,1);  
nodeIdx = zeros(numberOfRows,numberOfCols);  
nodeCounter = 0;
for rowNumber = 1:numberOfRows 
   for colNumber = 1:numberOfCols    
            
       if bwImage(rowNumber,colNumber) == 1
           nodeCounter = nodeCounter + 1;
           nodeIdx(rowNumber,colNumber) = nodeCounter; 
           nodeRow(nodeCounter) = rowNumber;    
           nodeCol(nodeCounter) = colNumber; 
       end
       
   end
end


% Create graph and add edges (links between nodes)
G = digraph;
for nodeCounter = 1:numberOfNodes  
       
       currentRow = nodeRow(nodeCounter);
       currentCol = nodeCol(nodeCounter);      
       if currentRow <= numberOfRows - 1 
           % find the nodes which are connected to the current node
           connRow = currentRow + 1;
           connCols = (currentCol - colSpan):(currentCol + colSpan);
           idx = find(connCols >=1 & connCols <= numberOfCols );
           connCols = connCols(idx);
           nonZeroIdx = find(bwImage(connRow,connCols) == 1);      
           connNodes = nodeIdx(connRow,connCols(nonZeroIdx));
            % create connections in the graph
           edgeWeights = weight.^(abs(connCols(nonZeroIdx)-currentCol));
           G = addedge(G,nodeCounter,connNodes,edgeWeights); 
  
       end
       
end


%  % Shortest path algorithm (choose path which minimizes total cost)
startNode = nodeIdx(startRow,startCol);
stopNode = nodeIdx(stopRow,stopCol);
if startNode == 0 | stopNode == 0
    disp('Warning: start or stop positions must be = 1 in black and white matrix')
    spRows = [];
    spCols = [];   
else
    nodesOptimalPath = shortestpath(G,startNode,stopNode,'Method','acyclic');
    spRows = nodeRow(nodesOptimalPath);
    spCols = nodeCol(nodesOptimalPath);   
end

end

