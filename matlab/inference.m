% load bike_model
% image_path = '/tmpVisionPJ/test_edges/10641.dat';
% %image_path = '/tmpVisionPJ/test_edges/770800.dat';
% %image_path = '/tmpVisionPJ/test_edges/10001.dat';
% %image_path = '/tmpVisionPJ/test_edges/10003.dat';
% image = readEdges(image_path);
% fprintf('Loaded %s (%dx%d)\n',image_path, size(image,1), size(image,2)); 
% mycost = costMap(image,bg_model,fg_model);


[best_cost,best_locations,dt_images] = kfan_matching(mycost, sigma, mu);

close all;
figure(1);
figure(2);
for i=1:6
    a = mycost(:,:,i)';
    figure(1);subplot(2,3,i);imagesc(a); caxis([-5000 2000]); colorbar;

    figure(2);subplot(2,3,i);imagesc(dt_images(:,:,i)); colorbar;
%     figure(2);subplot(2,3,i);imagesc(cp(:,:,i)); caxis([-5000 2000]); colorbar;
end


figure(3);
imagesc(image);colormap(gray);
hold on;
for i=1:6
    rectangle('Position',[best_locations(i,2),best_locations(i,1),2,2],'FaceColor','red')
end