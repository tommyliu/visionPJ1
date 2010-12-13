%runKadir

edge = '/tmpVisionPJ/test_edges/10001.dat';
edge = readEdges(edge);
realim = '/tmpVisionPJ/motorbikes/0001.jpg';
image = imread(realim);

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
Y=CalcScaleSaliency(uint8(double(image)./div),s1,s2, nbins, gsigma,AA)
fprintf('Clustering features...');
C=GreedyCluster(Y, wt,yt);
fprintf('done\n');
num_features=size(C,2);

figure(1);
imshow(image); hold on;
for i=1:num_features
  drawcircle(C(1,i)+1,C(2,i)+1,C(3,i)*2+1,'g');
end;