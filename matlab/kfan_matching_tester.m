load bike_model


rrow = 100;
rcol = 50;
HSIZE = 30;

costmaps = zeros(200,300,6);
t = fspecial('gaussian',HSIZE,10);


for i=1:6
    r = floor(rand*80);
    costmaps(r+floor(rrow+mu(i,1)):r+floor(rrow+HSIZE-1+mu(i,1)),r+floor(rcol+mu(i,2)):r+floor(rcol+HSIZE-1+mu(i,2)),i) = t;
    
    costmaps(:,:,i) = costmaps(:,:,i)/max(max(costmaps(:,:,i)));
end


[best_cost,best_locations] = kfan_matching(costmaps, sigma, mu)

image = sum(costmaps,3);
figure();
imagesc(image);
hold on;
for i=1:6
    rectangle('Position',[best_locations(i,2),best_locations(i,1),2,2],'FaceColor','white')
end