%Copyright 1998-2004 Timor Kadir.
%Kadir/Brady Feature detector (Scale Saliency) Code 
%CalcEntropyScaleOpt.c 
%For non-commericial use only.
s1=3; 			%Start Scale
s2=33;			%Stop Scale
AA=0;			%Anti-aliased sampling (not available with Parzen windowing).
nbins=16;		%number of bins (set to 0 for Parzen window PDF estimation)
gsigma=1;		%sigma for Parzen window (if nbins=0. Only available on 1D)
wt=0.5;                 %threshold on Saliency values
yt=0;                   %threshold on inter-scale saliency

%1D example
div=(255/(nbins-1));	%quantisation of image.
im=imread('image_0001.jpg');
Y=CalcScaleSaliency(uint8(double(im)./div),s1,s2, nbins, gsigma,AA);
fprintf('Clustering features...');
C=GreedyCluster(Y, wt,yt);
fprintf('done\n');
figure;
imshow(im);
length=size(C,2);
for i=1:length
  drawcircle(C(1,i)+1,C(2,i)+1,C(3,i)*2+1,'g');
end;




