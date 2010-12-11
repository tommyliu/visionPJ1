% Generate foreground appereance model
%
% model = backgroundModel(training_file, edges_folder, t)
%   training_file            - File with parts location per training image
%   edges_folder             - path to background edges folder
%   t                        - size of template
%   model       {txtx16}x6 - foreground model
%
% [NOTE]
% The model is a cell of 6 elements (6 parts). Each element is a matrix of
% txtx16 elements, where for each template we have probabilities for all
% possible edge orientations (from 0 to 15)
function model = foregroundModel(training_file, edges_folder, t)

model = {};
num_parts = 6;
mid_feat = t/2;

% get part locations
[part_locations, edges_file] = getPartLocations(training_file);

% initialize model
for p=1:num_parts, model{p} = zeros(t,t,16); end

% process each image
imid = 0;
for file_id=1:length(edges_file)    
    
    % read edges
    edges_path = fullfile(edges_file{file_id});
    edges = readEdges(edges_path);
    [M,N] = size(edges);    
    
    for p=1:num_parts               % work with each part independently
        
        location = part_locations(p,:,file_id);
        
        % extract feature window
        t_row = max(location(1)-mid_feat+1,1);    % top row
        b_row = min(location(1)+mid_feat,M);  % bottom row
        l_col = max(location(2)-mid_feat+1,1);      % left column
        r_col = min(location(2)+mid_feat,N);    % right column
        feature_window = edges(t_row:b_row, l_col:r_col);
        
%         imagesc(feature_window);
%         pause
        
        full_feature = zeros(t,t);
        full_t_row = t_row-(location(1)-mid_feat+1)+1;
        full_b_row = t-(location(1)+mid_feat-b_row);
        full_l_col = l_col-(location(2)-mid_feat+1)+1;
        full_r_col = t-(location(2)+mid_feat-r_col);
        full_feature(full_t_row:full_b_row,full_l_col:full_r_col) = ...
            feature_window;
                
%         imagesc(full_feature);
%         pause
        
        % count edges distribution
        for e=1:16
            model{p}(:,:,e) = ...
                model{p}(:,:,e) + (full_feature == e-1);
        end
        
    end
    
    imid = imid+1;
end

% turn the edges distribution into a probability function
for p=1:num_parts 
    feat_dist = model{p}(:,:,:);
    for r=1:t
       for c=1:t
           edge_dist = feat_dist(r,c,:);
           model{p}(r,c,:) = edge_dist./sum(edge_dist);
       end
    end
end

display(sprintf('\tProcessed %d images for the foreground model', imid));


