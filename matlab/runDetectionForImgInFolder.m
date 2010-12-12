% images path
images_folder = '/tmpVisionPJ/data_10/';
output_folder = 'results_10';

images = dir(images_folder);
for file_id=1:size(images,1)
    
    % don't process folders ni hidden files
    if (images(file_id).isdir), continue; end;
    image_path = images(file_id).name;
    if (image_path(1)=='.'), continue; end;
    display(sprintf('%s',image_path));
    
    image = readEdges(fullfile(images_folder,image_path));
   
    timerStart = tic;
    mycost = costMap(image,bg_model,fg_model);
    [best_cost,best_locations,dt_images] = kfan_matching(mycost, sigma, mu);
    timeElapsed = toc(timerStart);
    
    close all;
    h1 = figure(1);
    h2 = figure(2);
    for i=1:6
        a = mycost(:,:,i)';
        figure(1);subplot(2,3,i);imagesc(a); caxis([-5000 2000]); colorbar; axis equal;
        figure(2);subplot(2,3,i);imagesc(dt_images(:,:,i)); colorbar; axis equal;
    end

    saveas(h1,fullfile(output_folder,sprintf('%s_cost_map',basename(image_path))),'jpg');
    saveas(h2,fullfile(output_folder,sprintf('%s_distance',basename(image_path))),'jpg');

    h3 = figure(3);
    imagesc(image);colormap(gray);
    hold on;
    for i=1:6
        rectangle('Position',[best_locations(i,2)-5,best_locations(i,1)-5,5,5],'FaceColor','red')
    end
    
    saveas(h3,fullfile(output_folder,sprintf('%s_location',basename(image_path))),'jpg');
    results_file = sprintf('%s_result',basename(image_path));
    
    save(fullfile(output_folder,results_file),'image_path','mycost','best_cost', ...
        'best_locations','dt_images','timeElapsed');

end