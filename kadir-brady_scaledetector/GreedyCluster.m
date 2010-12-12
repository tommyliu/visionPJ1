%Copyright 1998-2004 Timor Kadir.
%Kadir/Brady Feature detector (Scale Saliency) Code  
%CalcEntropyScaleOpt.c 
%For non-commericial use only..

%Greedy clusterer for salient feature post-processing.

function C=GreedyCluster(Y,wt,yt)
C=[];
j=1;
[a b]=sort(-Y(6,:)); %Sort in ascending saliency
Y=Y(:,b);
MaxVal=Y(6,1);
ind=find(Y(6,:)>yt);
Y=Y(:,ind);
ind=find(Y(5,:)>wt);
Y=Y(:,ind);
while(size(Y,2)>0)
  Yi=Y(:,1);
  C(:,j)=Yi;
  j=j+1;
  a=repmat(Yi,[1 size(Y,2)]);
  dists=sqrt(sum((Y([1:2],:)-a([1:2],:)).^2));
  ind=find(dists>Yi(3));
  Y=Y(:,ind);
  %size(Y,2)
end;
