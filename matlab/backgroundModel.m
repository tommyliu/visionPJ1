% Generate background model
%
% model = backgroundModel(background_edges)
%   edges_folder     - path to background edges folder
%   model       1x16 - background model
%
% [NOTE]
% The model has 16 elements because we are considering 4 edge orientations.
% Refer to getSparseRepresentation for more information.
function model = backgroundModel(edges_folder)

model = zeros(1,16);

edges = dir(edges_folder);
num_edges = 0;
for file_id=1:size(edges,1)
   
    % don't process folders ni hidden files
    if (edges(file_id).isdir), continue; end;
    edge_mat = edges(file_id).name;
    if (edge_mat(1)=='.'), continue; end;
    % display(sprintf('Processing %s',edge_mat));
    
    % get edges mat
    edges_image = readEdges(fullfile(edges_folder,edge_mat));

    for e=1:16
       model(e) = sum(sum(edges_image == e-1));
    end
    
    num_edges = num_edges+1;
    
end

display(sprintf('Processed %d edge maps',num_edges));

% turn the edges distribution into a probability function
model = model/sum(model);
