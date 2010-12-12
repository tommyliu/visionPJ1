function costmaps = distance_transform(imgs,sigmas)

imgs = double(imgs);
[h,w,k] = size(imgs);
BSIZE =10;
costmaps = zeros(h,w,k);
for i=1:k
    costmaps(:,:,i) = conv2(imgs(:,:,i),fspecial('gaussian',BSIZE,3),'same');
end


% img = imresize(img,0.2);
% figure(1);
% subplot(2,2,1);imagesc(img);
% dist_img_total = zeros(h,w);
% for i=1:h
%     for j=1:w
%         dist_img_cur = zeros(h,w);
%         dist_img_cur(i,j) = 1;
%         img_cur = img(max(i-BSIZE,1):min(h,i+BSIZE),max(1,j-BSIZE):min(w,j+BSIZE));        
%         dist_img_cur = dist_img_cur(max(i-BSIZE,1):min(h,i+BSIZE),max(1,j-BSIZE):min(w,j+BSIZE));        
%         dist_img_cur = bwdist(dist_img_cur);
%         dist_img_cur = double(img_cur)./((dist_img_cur+1).^2);
%         dist_img_total(i,j) =sum(sum(dist_img_cur));
%     end
% end
% 
% 
% dist_img = img>(max(max(img))-min(min(img)))*0.6;
% dist_img_total = bwdist(dist_img);
% dist_img_total = 1-dist_img_total;
% display(sprintf('(elapsed_time: %.2f)', toc(tStartS)));


% subplot(2,2,2);imagesc(dist_img_total);
% subplot(2,2,3);imagesc(gaussian_blur);