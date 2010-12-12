hard_bg_model =  [ 0.4907, 0.0753, 0.0223, 0.0275, 0.0927,...
0.0152,    0.0384,    0.0272,...
0.0229,    0.0272,    0.0070,    0.0228,...
0.0380, 0.0270, 	0.0283,    0.0376 ];

mycost = costMap(image_001,hard_bg_model,fg_model);

close all;
figure(1);

for i=1:6
    a = mycost(:,:,i)';
    figure(1);subplot(2,3,i);imagesc(a); caxis([-5000 2000]); colorbar;
    figure(2);subplot(2,3,i);imagesc(cp(:,:,i)); caxis([-5000 2000]); colorbar;
end

