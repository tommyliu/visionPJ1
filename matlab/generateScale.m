% Compute scale for images (edges) in learning data set
edges_folder = '/tmpVisionPJ/test_edges/';
motorbikes_folder = '/tmpVisionPJ/motorbikes/';
output_folder = 'kadir_features_test';
training_file = '../motorbike_training.dat';

[partLocations,edgesPath] = getPartLocations(training_file);

edges = dir(edges_folder);

scales = zeros(15,3);
% scales = zeros(length(imagesPath),3);

% for file_id=1:size(edges,1)

for file_id=1:15  

    image_name = edgesPath{file_id};
    image_basename = basename(image_name);
    display(sprintf('%s',image_basename));
    
    edge = readEdges(image_name);
    image = imread(fullfile(motorbikes_folder,...
                    sprintf('%s.jpg',image_basename(2:end))));
                
    [h,w] = size(edge);
    [M,N,C] = size(image);
    image = imresize(image, [h w]);
    
   
    s1=5; 			%Start Scale
    s2=20;			%Stop Scale
    AA=0;			%Anti-aliased sampling (not available with Parzen windowing).
    nbins=16;		%number of bins (set to 0 for Parzen window PDF estimation)
    gsigma=1;		%sigma for Parzen window (if nbins=0. Only available on 1D)
    wt=0.6;                 %threshold on Saliency values
    yt=0;                   %threshold on inter-scale saliency

    %1D example
    div=(255/(nbins-1));	%quantisation of image.
    Y=CalcScaleSaliency(uint8(double(image)./div),s1,s2, nbins, gsigma,AA);
    fprintf('Clustering features...');
    C=GreedyCluster(Y, wt,yt);
    fprintf('done\n');
    num_features=size(C,2);
    
    % find closes point to reference
    reference = partLocations(1,1:2,file_id);
    
    tmp = C(1:2,:)'-repmat(reference,num_features,1);
    distances = sum(tmp.^2,2);
    [dummy,i] = min(distances);
    scales(file_id,:) = C(1:3,i)';
    display(sprintf('%.2f %.2f %.2f', ...
        scales(file_id,1), scales(file_id,2), scales(file_id,3)));
    
    myfigure=figure(1);
    imshow(image); hold on;
    for i=1:num_features
      drawcircle(C(1,i)+1,C(2,i)+1,C(3,i)*2+1,'g');
    end;
    
    output_file = fullfile(output_folder,image_basename);
    saveas(myfigure,output_file,'jpg');
    
end

save(fullfile(output_folder,'scales_test.dat'),'scales');