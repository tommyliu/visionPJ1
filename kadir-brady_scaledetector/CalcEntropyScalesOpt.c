/*
Copyright 1998-2004 Timor Kadir.

Kadir/Brady Feature detector (Scale Saliency) Code 

For non-commericial use only.
*/

#include <math.h>
#include <string.h>
#include "mex.h"
#include "mat.h"
#include "General.h"
#include "Hists.h"
#ifdef FASTPLOG
#define PLOGPRES 10000
double PLOGP[PLOGPRES];
#endif

void PrintCopyright()
{
   mexPrintf("Kadir/Brady Feature Detector (Scale Saliency) Version 1.5. Copyright Timor Kadir 2004.\n");
   mexPrintf("Matlab(TM) binary for Non Commercial use only.\n");
}


mxArray ** CirclePix2(int StartS, int StopS, int m)
{
   int roi_dims[2], centre, i,j,s, count=0, scale;
   double radius2, Prevradius2=-1, z;
   uint32_T *pScalePix;
   mxArray **pROIPix;

   int MaxScale=(StopS*2)+1;
   roi_dims[0]=MaxScale*MaxScale;
   centre=(MaxScale+1)/2;

   pROIPix=calloc(StopS+1, sizeof(mxArray *) );

   for (s=StartS;s<=StopS;s++){
      pROIPix[s]=mxCreateNumericArray(1,roi_dims,mxUINT32_CLASS,mxREAL);
      pScalePix=(UINT32_T *) mxGetPr(pROIPix[s]);
      scale=(s*2)+1;
      radius2=pow((scale-1)/2,2);
      for (j=1;j<=MaxScale;j++)
         for (i=1;i<=MaxScale;i++){
            z=pow((i-centre),2) +pow((j-centre),2);
            if (z<=radius2 && z>Prevradius2) {
               pScalePix[count++]=i-1+(j-1)*m;
            }
         }
         mxSetM(pROIPix[s],(count));
         count=0;
         Prevradius2=radius2;
   }
   return(pROIPix);
}

void AACirclePix(int StartS, int StopS, int m, mxArray ***pROIP, mxArray ***pROIW)
{
   int i,j,s, ind=0, range,dims[2],*pP,*pW;
   double z,w;

   range=ceil(StopS)+1;
   *pROIP=calloc(StopS+1, sizeof(mxArray *));
   *pROIW=calloc(StopS+1, sizeof(mxArray *));
   dims[0]=2*range+1;
   dims[1]=2*range+1;

   for (s=StartS;s<=StopS;s++){
      ind=0;
      (*pROIP)[s]=mxCreateNumericArray (2,dims,mxUINT32_CLASS,mxREAL);
      (*pROIW)[s]=mxCreateNumericArray (2,dims,mxUINT32_CLASS,mxREAL);
      pP=(int *) mxGetPr((*pROIP)[s]);
      pW=(int *) mxGetPr((*pROIW)[s]);
      for (j=-range;j<=range;j++)
         for (i=-range;i<=range;i++){
            z=sqrt(pow(i,2) +pow(j,2));
            w=((1/(1.0 + pow(z/s,42)))-(1/(1.0 + pow(z/(s-1),42))));
            if (w>0.0001){
               pP[ind]=i+j*m;
               pW[ind++]=(int)rint(w*1000);
            }
         }
         mxSetM((*pROIP)[s],(ind));
         mxSetM((*pROIW)[s],(ind));
   }
}


inline double plogp (double p)
{
   if (p==0)
      return(0);
   else
      return(p*log(p));
}

void CreateParzen(double GSigma, double *GMean, double **FGauss)
{
   int size,i;
   size=GSigma*3;
   *GMean=size/2;
   mexPrintf("size %d\n",size);
   *FGauss=(double *)calloc(size, sizeof(double));

   for(i=0;i<size;i++){
      (*FGauss)[i]=(1./(GSigma*sqrt(2*3.142)))*exp(-pow(((double)i-(*GMean)+1)/(GSigma),2));
   }
}

void inline ParzenSmooth(int *Histo1, double *HistoSmooth, double *FGauss, double GSigma, double GMean)
{
   int istart,istop,i,b;
   double temp;

   for(b=0;b<256;b++){
      temp=0;
      istart=-MIN(0,b-GMean);
      istop=MIN(GMean*2,256-b+GMean-1);
      for(i=istart;i<istop;i++){
         temp+=(double)Histo1[b+i-(int)GMean]*FGauss[i];
      }
      HistoSmooth[b]=temp;
   }
}

mxArray *do_CalcSalScale1D(const mxArray *Image, int s_start,int s_stop, int nbins)
{
   int Pointx, Pointy, scale,x,y,j, *Histo1, s,hs,i,m,n, counts, *peaks, *lengths, sum, *pPNZ;
   int x_start, x_stop, y_start, y_stop, NoOfScales=0, ind=0;
   double tempNormBin, *Histo2, *tempPr, progress, *pEntropyArray, *pDistArray, *pBestSaliency, Norm, *pImageDouble;
   uint8_T *pImage, *pROI_icon;
   mxArray *rhs[2], *lhs, **ROIPix, *ProgressBar, *EntropyArray, *DistArray, *BestSaliency;
   uint32_T **pROIPix;

   pImage = (uint8_T *) mxGetPr(Image); /* the input image*/
   m=mxGetM(Image);
   n=mxGetN(Image);

   Histo1=(int *) calloc(nbins, sizeof(int));
   Histo2=(double *) calloc(nbins, sizeof(double));
   pPNZ=(int *) calloc(nbins, sizeof(int));

   s=(s_stop*2)+1;
   hs=(s-1)/2;
   x_start=hs;
   x_stop=m-hs;

   y_start=x_start;
   y_stop=n-hs;

   /* setup entropy output array */
   EntropyArray= mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pEntropyArray= (double *) mxGetPr(EntropyArray);

   /*setup distance output array*/
   DistArray=mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pDistArray= (double *) mxGetPr(DistArray);

   /* Create output matrix of bestscales and their saliency values */
   BestSaliency=mxCreateDoubleMatrix(6,100000,mxREAL);
   pBestSaliency=( double *) mxGetPr(BestSaliency);

   peaks=calloc(s_stop, sizeof(int));
   ROIPix=CirclePix2(s_start, s_stop,m);
   lengths=calloc(s_stop+1, sizeof(int *));
   pROIPix=calloc(s_stop+1, sizeof(uint32_T **));

   for (scale=s_start;scale<=s_stop;scale++){
      lengths[scale]=mxGetM(ROIPix[scale]);
      pROIPix[scale]=(uint32_T *) mxGetPr(ROIPix[scale]);
   }
#ifdef SHOWPROG
   // Just the progress bar stuff
   rhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
   tempPr=(double *) mxGetPr(rhs[0]);
   *tempPr=0;
   rhs[1]=mxCreateString("Calculating..");
   mexCallMATLAB(1, &lhs, 2, rhs, "waitbar");
   rhs[1]=lhs;
   progress=1/(double)(x_stop-x_start);
#endif

   for (j=0;j<nbins;j++){
      Histo1[j]=0; 
      Histo2[j]=0;
   }
   for (Pointx=x_start;Pointx<x_stop;Pointx++){
#ifdef SHOWPROG
      *tempPr+=progress;
      mexCallMATLAB(0, NULL, 2, rhs, "waitbar");
#endif
      for (Pointy=y_start;Pointy<y_stop;Pointy++){
         x=Pointx-hs;
         y=Pointy-hs;
         pROI_icon=pImage+x+(y)*m;
         sum=0;
         ind=0;
         for (scale=s_start;scale<=s_stop;scale++){
            HistUint8_1DOpt4(pROI_icon, Histo1, pROIPix[scale], lengths[scale]);
            sum+=lengths[scale];
            pEntropyArray[scale]=0;
            pDistArray[scale]=0;
            for (j=0;j<nbins;j++){
               if (Histo1[j]>0){
                  tempNormBin=(double)Histo1[j]/sum;
#ifdef FASTPLOG
                  pEntropyArray[scale]-=PLOGP[(int)(tempNormBin*PLOGPRES)];
#else
                  pEntropyArray[scale]-=plogp(tempNormBin);
#endif
                  pDistArray[scale]+=fabs((tempNormBin)-(Histo2[j]));
                  Histo2[j]=tempNormBin;
               }
               else{
                  pDistArray[scale]+=Histo2[j];
                  Histo2[j]=0;      
               }
            }
         }
         counts=0;
         pDistArray[s_start]=pDistArray[s_start+1];

         /* Smooth the DistArray */
         for (scale=s_start+1;scale<s_stop;scale++)
            pDistArray[scale]=(pDistArray[scale-1]+pDistArray[scale]+pDistArray[scale+1])/3;

         for (scale=s_start+1;scale<s_stop;scale++){
            if ((pEntropyArray[scale]>pEntropyArray[scale-1]) &&
               (pEntropyArray[scale]>pEntropyArray[scale+1])){
                  peaks[counts]=scale;
                  counts++;
               }
         }
         if (NoOfScales+counts>=mxGetN(BestSaliency)){
            mxSetN(BestSaliency,NoOfScales+counts+100000);
            pBestSaliency=mxRealloc(pBestSaliency,6*(NoOfScales+counts+100000)*sizeof(double));
            mxSetPr(BestSaliency,pBestSaliency);
         }
         /*  assign the best (peaks) scales and global saliency for this x,y location */
         if (counts>0)
            for (i=NoOfScales;i<NoOfScales+counts;i++){
               pBestSaliency[i*6+0]=Pointx;  /* the x location */
               pBestSaliency[i*6+1]=Pointy;  /* the y location */
               pBestSaliency[i*6+2]=peaks[i-NoOfScales];	  /*the best scale*/
               Norm=pow(pBestSaliency[i*6+2],2)/(2*pBestSaliency[i*6+2]-1);
               pBestSaliency[i*6+3]=pEntropyArray[(peaks[i-NoOfScales])];	  /*HD*/
               pBestSaliency[i*6+4]=pDistArray[(peaks[i-NoOfScales])]*Norm;	  /*WD*/
               pBestSaliency[i*6+5]=Norm*pEntropyArray[(peaks[i-NoOfScales])]*(pDistArray[(peaks[i-NoOfScales])]); /* the saliency */
            }

            NoOfScales+=counts;
            for (j=0;j<nbins;j++){
               Histo1[j]=0;
               Histo2[j]=0;
            }
      }/*end of y loop*/
   }/*end of x loop*/
   mxSetN(BestSaliency,NoOfScales);
#ifdef SHOWPROG
   mexCallMATLAB(0, NULL, 1, &lhs, "close");
#endif
   free(peaks);
   free(ROIPix);
   free(lengths);
   free(Histo1);
   free(Histo2);
   free(pPNZ);
   return(BestSaliency);
}

mxArray *do_CalcSalScale1DAA(const mxArray *Image, int s_start,int s_stop, int nbins)
{
   int Pointx, Pointy, scale,x,y,j,s,hs,i,m,n, counts, *peaks, *lengths, *sums, *pPNZ,*Histo1;
   int x_start, x_stop, y_start, y_stop, NoOfScales=0;
   double tempNormBin, *Histo2, *tempPr, progress, *pEntropyArray, *pDistArray, *pBestSaliency, Norm,sum;
   uint8_T *pImage, *pROI_icon;
   mxArray *rhs[2], *lhs, **ROIP, **ROIW, *ProgressBar, *EntropyArray, *DistArray, *BestSaliency;
   int **pROIP,**pROIW;
   pImage = (uint8_T *) mxGetPr(Image); /* the input image*/
   m=mxGetM(Image);
   n=mxGetN(Image);
   Histo1=(int *) calloc(nbins, sizeof(int));
   Histo2=(double *) calloc(nbins, sizeof(double));
   pPNZ=(int *) calloc(nbins, sizeof(int));

   s=(s_stop+1)*2+1;
   hs=(s-1)/2;
   x_start=hs;
   x_stop=m-hs;

   y_start=x_start;
   y_stop=n-hs;

   /* setup entropy output array */
   EntropyArray= mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pEntropyArray= (double *) mxGetPr(EntropyArray);

   /*setup distance output array*/
   DistArray=mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pDistArray= (double *) mxGetPr(DistArray);

   /* Create output matrix of bestscales and their saliency values */
   BestSaliency=mxCreateDoubleMatrix(6,100000,mxREAL);
   pBestSaliency=( double *) mxGetPr(BestSaliency);

   peaks=calloc(s_stop, sizeof(int));
   AACirclePix(s_start, s_stop,m, &ROIP, &ROIW);
   lengths=calloc(s_stop+1, sizeof(int *));
   sums=calloc(s_stop+1, sizeof(int *));
   pROIP=calloc(s_stop+1, sizeof(uint32_T **));
   pROIW=calloc(s_stop+1, sizeof(uint32_T **));

   for (scale=s_start;scale<=s_stop;scale++){
      lengths[scale]=mxGetM(ROIP[scale]);
      pROIP[scale]=(int *) mxGetPr(ROIP[scale]);
      pROIW[scale]=(int *) mxGetPr(ROIW[scale]);
      sums[scale]=0;
      for (j=0;j<lengths[scale];j++){
         sums[scale]+=*(pROIW[scale]+j);
      }
   }
#ifdef SHOWPROG
   // Just the progress bar stuff
   rhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
   tempPr=(double *) mxGetPr(rhs[0]);
   *tempPr=0;
   rhs[1]=mxCreateString("Calculating..");
   mexCallMATLAB(1, &lhs, 2, rhs, "waitbar");
   rhs[1]=lhs;
   progress=1/(double)(x_stop-x_start);
#endif

   for (j=0;j<nbins;j++){
      Histo1[j]=0; 
      Histo2[j]=0;
   }
   for (Pointx=x_start;Pointx<x_stop;Pointx++){
#ifdef SHOWPROG
      *tempPr+=progress;
      mexCallMATLAB(0, NULL, 2, rhs, "waitbar");
#endif
      for (Pointy=y_start;Pointy<y_stop;Pointy++){
         x=Pointx;
         y=Pointy;
         pROI_icon=pImage+x+(y)*m;
         sum=0;
         for (scale=s_start;scale<=s_stop;scale++){
            HistUint8_1DOpt4AA(pROI_icon, Histo1, pROIP[scale],pROIW[scale], lengths[scale]-1);
            sum+=sums[scale];
            pEntropyArray[scale]=0;
            pDistArray[scale]=0;
            for (j=0;j<nbins;j++){
               if (Histo1[j]>0){
                  tempNormBin=(double)Histo1[j]/sum;
#ifdef FASTPLOG
                  pEntropyArray[scale]-=PLOGP[(int)(tempNormBin*PLOGPRES)];
#else
                  pEntropyArray[scale]-=plogp(tempNormBin);
#endif
                  pDistArray[scale]+=fabs((tempNormBin)-(Histo2[j]));
                  Histo2[j]=tempNormBin;
               }
               else{
                  pDistArray[scale]+=Histo2[j];
                  Histo2[j]=0;      
               }
            }
         }
         counts=0;
         pDistArray[s_start]=pDistArray[s_start+1];

         for (scale=s_start+1;scale<s_stop;scale++){
            pDistArray[scale]=(pDistArray[scale-1]+pDistArray[scale]+pDistArray[scale+1])/3;
            if ((pEntropyArray[scale]>pEntropyArray[scale-1]) &&
               (pEntropyArray[scale]>pEntropyArray[scale+1])){
                  peaks[counts]=scale;
                  counts++;
               }
         }
         if (NoOfScales+counts>=mxGetN(BestSaliency)){
            mxSetN(BestSaliency,NoOfScales+counts+100000);
            pBestSaliency=mxRealloc(pBestSaliency,6*(NoOfScales+counts+100000)*sizeof(double));
            mxSetPr(BestSaliency,pBestSaliency);
         }
         /*  assign the best (peaks) scales and global saliency for this x,y location */
         if (counts>0)
            for (i=NoOfScales;i<NoOfScales+counts;i++){
               pBestSaliency[i*6+0]=Pointx;  /* the x location */
               pBestSaliency[i*6+1]=Pointy;  /* the y location */
               pBestSaliency[i*6+2]=peaks[i-NoOfScales];	  /*the best scale*/
               // TODO: fix scale normalisation
               Norm=peaks[i-NoOfScales];//pow(pBestSaliency[i*6+2],2)/(2*pBestSaliency[i*6+2]-1);
               pBestSaliency[i*6+3]=pEntropyArray[(peaks[i-NoOfScales])];	  /*HD*/
               pBestSaliency[i*6+4]=pDistArray[(peaks[i-NoOfScales])]*Norm;	  /*WD*/
               pBestSaliency[i*6+5]=Norm*pEntropyArray[(peaks[i-NoOfScales])]*(pDistArray[(peaks[i-NoOfScales])]); /* the saliency */
            }

            NoOfScales+=counts;
            for (j=0;j<nbins;j++){
               Histo1[j]=0;
               Histo2[j]=0;
            }
      }/*end of y loop*/
   }/*end of x loop*/
   mxSetN(BestSaliency,NoOfScales);
#ifdef SHOWPROG
   mexCallMATLAB(0, NULL, 1, &lhs, "close");
#endif
   free(peaks);
   free(ROIP);
   free(ROIW);
   free(lengths);
   free(Histo1);
   free(Histo2);
   free(pPNZ);
   free(pROIP);
   free(pROIW);
   free(sums);
   return(BestSaliency);
}

mxArray *do_CalcSalScale1DParzen(const mxArray *Image, int s_start,int s_stop, double GSigma)
{
   int Pointx, Pointy, scale,x,y,k, *Histo1, s,hs,i,m,n, counts, *peaks, *lengths, sum;
   int x_start, x_stop, y_start, y_stop, NoOfScales=0;
   double tempNormBin, *Histo2, *HistoSmooth, *tempPr, progress, *pEntropyArray, *pDistArray, *pBestSaliency, GMean, *FGauss, Norm;
   uint8_T *pImage, *pROI_icon;
   mxArray *rhs[2], *lhs, **ROIPix, *ProgressBar, *EntropyArray, *DistArray, *BestSaliency;
   uint32_T **pROIPix;
   int nbins=256;
   pImage = (uint8_T *) mxGetPr(Image); /* the input image*/
   m=mxGetM(Image);
   n=mxGetN(Image);
   Histo1=(int *) calloc(nbins, sizeof(int));
   Histo2=(double *) calloc(nbins, sizeof(double));
   HistoSmooth=(double *) calloc(nbins, sizeof(double));

   s=(s_stop*2)+1;
   hs=(s-1)/2;
   x_start=hs;
   x_stop=m-hs;

   y_start=x_start;
   y_stop=n-hs;

   /*Create Smoothing filter for Parzen Post Smoothing*/
   CreateParzen(GSigma, &GMean, &FGauss);

   /* setup entropy output array */
   EntropyArray= mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pEntropyArray= (double *) mxGetPr(EntropyArray);

   /*setup distance output array*/
   DistArray=mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pDistArray= (double *) mxGetPr(DistArray);

   /* Create output matrix of bestscales and their saliency values */
   BestSaliency=mxCreateDoubleMatrix(6,100000,mxREAL);
   pBestSaliency=( double *) mxGetPr(BestSaliency);

   peaks=calloc(s_stop, sizeof(int));
   ROIPix=CirclePix2(s_start, s_stop,m);
   lengths=calloc(s_stop+1, sizeof(int *));
   pROIPix=calloc(s_stop+1, sizeof(uint32_T **));

   for (scale=s_start;scale<=s_stop;scale++){
      lengths[scale]=mxGetM(ROIPix[scale]);
      pROIPix[scale]=(uint32_T *) mxGetPr(ROIPix[scale]);
   }

#ifdef SHOWPROG
   // Just the progress bar stuff
   rhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
   tempPr=(double *) mxGetPr(rhs[0]);
   *tempPr=0;
   rhs[1]=mxCreateString("Calculating..");
   mexCallMATLAB(1, &lhs, 2, rhs, "waitbar");
   rhs[1]=lhs;
   progress=1/(double)(x_stop-x_start);
#endif

   for (k=0;k<nbins;k++){
      Histo1[k]=0; 
      Histo2[k]=0;
   }
   for (Pointx=x_start;Pointx<x_stop;Pointx++){
#ifdef SHOWPROG
      *tempPr+=progress;
      mexCallMATLAB(0, NULL, 2, rhs, "waitbar");
#endif
      for (Pointy=y_start;Pointy<y_stop;Pointy++){
         x=Pointx-hs;
         y=Pointy-hs;
         pROI_icon=pImage+x+(y)*m;
         sum=0;
         for (scale=s_start;scale<=s_stop;scale++){
            int sind=scale-1;
            HistUint8_1DOpt3Parz(pROI_icon, Histo1, pROIPix[scale], lengths[scale]);
            ParzenSmooth(Histo1, HistoSmooth, FGauss, GSigma, GMean);
            sum=0;
            for (k=0;k<nbins;k++) sum+=Histo1[k];
            pEntropyArray[sind]=0;
            pDistArray[sind]=0;
            for (k=0;k<nbins;k++){
               tempNormBin=(double)HistoSmooth[k]/sum;
#ifdef FASTPLOG				
               pEntropyArray[sind]-=PLOGP[(int)(tempNormBin*PLOGPRES)];
#else
               pEntropyArray[sind]-=plogp(tempNormBin);
#endif
               pDistArray[sind]+=fabs((tempNormBin)-(Histo2[k]));
               Histo2[k]=tempNormBin;
            }
         }
         counts=0;
         /* Smooth the DistArray */
         for (scale=s_start+1;scale<s_stop;scale++)
            pDistArray[scale]=(pDistArray[scale-1]+pDistArray[scale]+pDistArray[scale+1])/3;

         for (scale=s_start;scale<=s_stop-2;scale++){
            if ((pEntropyArray[scale]<pEntropyArray[scale+1]) &&
               (pEntropyArray[scale+1]>pEntropyArray[scale+2])){
                  peaks[counts]=scale+1;
                  counts++;
               }
         }
         if (NoOfScales+counts>=mxGetN(BestSaliency)){
            mxSetN(BestSaliency,NoOfScales+counts+100000);
            pBestSaliency=mxRealloc(pBestSaliency,6*(NoOfScales+counts+100000)*sizeof(double));
            mxSetPr(BestSaliency,pBestSaliency);
         }
         /*  assign the best (peaks) scales and global saliency for this x,y location */
         if (counts>0)
            for (i=NoOfScales;i<NoOfScales+counts;i++){
               pBestSaliency[i*6+0]=Pointx;  /* the x location */
               pBestSaliency[i*6+1]=Pointy;  /* the y location */
               pBestSaliency[i*6+2]=peaks[i-NoOfScales]+1;	  /*the best scale*/
               Norm=pow(pBestSaliency[i*6+2],2)/(2*pBestSaliency[i*6+2]-1);
               pBestSaliency[i*6+3]=pEntropyArray[(peaks[i-NoOfScales])];	  /*HD*/
               pBestSaliency[i*6+4]=pDistArray[(peaks[i-NoOfScales])]*Norm;	  /*WD*/
               pBestSaliency[i*6+5]=Norm*pEntropyArray[(peaks[i-NoOfScales])]*(pDistArray[(peaks[i-NoOfScales])]); /* the saliency */
            }

            NoOfScales+=counts;
            for (k=0;k<nbins;k++){
               Histo1[k]=0;
               HistoSmooth[k]=0;
               Histo2[k]=0;
            }
      }/*end of y loop*/
   }/*end of x loop*/
   mxSetN(BestSaliency,NoOfScales);
#ifdef SHOWPROG
   mexCallMATLAB(0, NULL, 1, &lhs, "close");
#endif
   free(peaks);
   free(ROIPix);
   free(lengths);
   return(BestSaliency);
}


mxArray *do_CalcSalScale2D(const mxArray *Image, int s_start,int s_stop, int nbins)
{
   int Pointx, Pointy, scale,x,y,j, *Histo1, s,hs,i,m,n, counts, *peaks, *lengths, sum, ind=0, *pPNZ;
   int x_start, x_stop, y_start, y_stop, NoOfScales=0, binmax=0;
   double tempNormBin, *Histo2, *tempPr, progress, *pEntropyArray, *pDistArray, *pBestSaliency, Norm;
   uint8_T *pImage1D, *pImage2D, *pROI_icon1D, *pROI_icon2D;
   mxArray *rhs[2], *lhs, **ROIPix, *ProgressBar, *EntropyArray, *DistArray, *BestSaliency;
   uint32_T **pROIPix;

   m=mxGetM(Image);
   n=mxGetN(Image)/2;
   pImage1D = (uint8_T *) mxGetPr(Image); /* the input image*/
   pImage2D = pImage1D+m*n;
   binmax=nbins*nbins;
   Histo1=(int *) calloc(binmax, sizeof(int));
   Histo2=(double *) calloc(binmax, sizeof(double));
   pPNZ=(int *) calloc(binmax, sizeof(int));

   s=(s_stop*2)+1;
   hs=(s-1)/2;
   x_start=hs;
   x_stop=m-hs;

   y_start=x_start;
   y_stop=n-hs;

   /* setup entropy output array */
   EntropyArray= mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pEntropyArray= (double *) mxGetPr(EntropyArray);

   /*setup distance output array*/
   DistArray=mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pDistArray= (double *) mxGetPr(DistArray);

   /* Create output matrix of bestscales and their saliency values */
   BestSaliency=mxCreateDoubleMatrix(6,100000,mxREAL);
   pBestSaliency=( double *) mxGetPr(BestSaliency);

   peaks=calloc(s_stop, sizeof(int));
   ROIPix=CirclePix2(s_start, s_stop,m);
   lengths=calloc(s_stop+1, sizeof(int *));
   pROIPix=calloc(s_stop+1, sizeof(uint32_T **));

   for (j=0;j<binmax;j++){
      Histo1[j]=0; 
      Histo2[j]=0;
   }
   for (scale=s_start;scale<=s_stop;scale++){
      lengths[scale]=mxGetM(ROIPix[scale]);
      pROIPix[scale]=(uint32_T *) mxGetPr(ROIPix[scale]);
   }

#ifdef SHOWPROG
   // Just the progress bar stuff
   rhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
   tempPr=(double *) mxGetPr(rhs[0]);
   *tempPr=0;
   rhs[1]=mxCreateString("Calculating..");
   mexCallMATLAB(1, &lhs, 2, rhs, "waitbar");
   rhs[1]=lhs;
   progress=1/(double)(x_stop-x_start);
#endif

   for (Pointx=x_start;Pointx<x_stop;Pointx++){
#ifdef SHOWPROG
      *tempPr+=progress;
      mexCallMATLAB(0, NULL, 2, rhs, "waitbar");
#endif
      for (Pointy=y_start;Pointy<y_stop;Pointy++){
         x=Pointx-hs;
         y=Pointy-hs;
         pROI_icon1D=pImage1D+x+(y)*m;
         pROI_icon2D=pImage2D+x+(y)*m;
         sum=0;
         for (scale=s_start;scale<=s_stop;scale++){
            int sind=scale-1;
            double isum;
            HistUint8_2DOpt4(pROI_icon1D, pROI_icon2D, Histo1, pROIPix[scale], lengths[scale], nbins);
            sum+=lengths[scale];
            pEntropyArray[sind]=0;
            pDistArray[sind]=0;
            isum=1./sum; 
            for (j=0;j<ind;j++)
               if (Histo1[pPNZ[j]]==0){ 	//this optimised so that it only does this when Histo2>0 (a bin is in PNZ)
                  double r=Histo2[pPNZ[j]];
                  r+=pDistArray[sind];
                  Histo2[pPNZ[j]]=0;
                  pDistArray[sind]=r;
               }
               ind=0;
               for (j=0;j<binmax;j++){
                  if (Histo1[j]!=0){
                     double r=Histo1[j];
                     pPNZ[ind++]=j;
                     tempNormBin=r*isum;
#ifdef FASTPLOG				
                     pEntropyArray[sind]-=PLOGP[(int)(tempNormBin*PLOGPRES)];
#else
                     pEntropyArray[sind]-=plogp(tempNormBin);
#endif
                     pDistArray[sind]+=fabs((tempNormBin)-(Histo2[j]));
                     Histo2[j]=tempNormBin;
                  }
               }
         }
         counts=0;
         /* Smooth the DistArray */
         for (scale=s_start+1;scale<s_stop;scale++)
            pDistArray[scale]=(pDistArray[scale-1]+pDistArray[scale]+pDistArray[scale+1])/3;

         for (scale=s_start;scale<=s_stop-2;scale++){
            if ((pEntropyArray[scale]<pEntropyArray[scale+1]) &&
               (pEntropyArray[scale+1]>pEntropyArray[scale+2])){
                  peaks[counts]=scale+1;
                  counts++;
               }
         }
         if (NoOfScales+counts>=mxGetN(BestSaliency)){
            mxSetN(BestSaliency,NoOfScales+counts+100000);
            pBestSaliency=mxRealloc(pBestSaliency,6*(NoOfScales+counts+100000)*sizeof(double));
            mxSetPr(BestSaliency,pBestSaliency);
         }
         /*  assign the best (peaks) scales and global saliency for this x,y location */
         if (counts>0)
            for (i=NoOfScales;i<NoOfScales+counts;i++){
               pBestSaliency[i*6+0]=Pointx;  /* the x location */
               pBestSaliency[i*6+1]=Pointy;  /* the y location */
               pBestSaliency[i*6+2]=peaks[i-NoOfScales]+1;	  /*the best scale*/
               Norm=pow(pBestSaliency[i*6+2],2)/(2*pBestSaliency[i*6+2]-1);
               pBestSaliency[i*6+3]=pEntropyArray[(peaks[i-NoOfScales])];	  /*HD*/
               pBestSaliency[i*6+4]=pDistArray[(peaks[i-NoOfScales])]*Norm;	  /*WD*/
               pBestSaliency[i*6+5]=Norm*pEntropyArray[(peaks[i-NoOfScales])]*(pDistArray[(peaks[i-NoOfScales])]); /* the saliency */
            }

            NoOfScales+=counts;
            for (j=0;j<binmax;j++){
               Histo1[j]=0;
               Histo2[j]=0;
            }
      }/*end of y loop*/
   }/*end of x loop*/
   mxSetN(BestSaliency,NoOfScales);
#ifdef SHOWPROG
   mexCallMATLAB(0, NULL, 1, &lhs, "close");
#endif
   free(peaks);
   free(ROIPix);
   free(lengths);
   free(Histo1);
   free(Histo2);
   free(pPNZ);
   return(BestSaliency);
}


mxArray *do_CalcSalScale3D(const mxArray *Image, int s_start,int s_stop, int nbins)
{
   int Pointx, Pointy, scale,x,y,j,o, *Histo1, s,hs,i,m,n, counts, *peaks, *lengths, *pPNZ, sum;
   int x_start, x_stop, y_start, y_stop, NoOfScales=0, binmax=0, ind=0;
   double tempNormBin=0, *Histo2, *tempPr, progress, *pEntropyArray, *pDistArray, *pBestSaliency, Norm, isum;
   uint8_T *pImage1D, *pImage2D, *pImage3D, *pROI_icon1D, *pROI_icon2D, *pROI_icon3D;
   mxArray *rhs[2], *lhs, **ROIPix, *ProgressBar, *EntropyArray, *DistArray, *BestSaliency;
   uint32_T **pROIPix;


   m=mxGetM(Image);
   n=mxGetN(Image)/3;
   pImage1D = (uint8_T *) mxGetPr(Image); /* the input image*/
   pImage2D = pImage1D+m*n;
   pImage3D = pImage1D+(m*n*2);
   binmax=nbins*nbins*nbins;
   Histo1=(int *) calloc(binmax, sizeof(int));
   Histo2=(double *) calloc(binmax, sizeof(double));
   pPNZ=(int *) calloc(binmax, sizeof(int));

   s=(s_stop*2)+1;
   hs=(s-1)/2;
   x_start=hs;
   x_stop=m-hs;

   y_start=x_start;
   y_stop=n-hs;

   /* setup entropy output array */
   EntropyArray= mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pEntropyArray= (double *) mxGetPr(EntropyArray);

   /*setup distance output array*/
   DistArray=mxCreateDoubleMatrix (1,s_stop+1,mxREAL);
   pDistArray= (double *) mxGetPr(DistArray);

   /* Create output matrix of bestscales and their saliency values */
   BestSaliency=mxCreateDoubleMatrix(6,100000,mxREAL);
   pBestSaliency=( double *) mxGetPr(BestSaliency);

   peaks=calloc(s_stop, sizeof(int));
   ROIPix=CirclePix2(s_start, s_stop,m);
   lengths=calloc(s_stop+1, sizeof(int *));
   pROIPix=calloc(s_stop+1, sizeof(uint32_T **));

   for (scale=s_start;scale<=s_stop;scale++){
      lengths[scale]=mxGetM(ROIPix[scale]);
      pROIPix[scale]=(uint32_T *) mxGetPr(ROIPix[scale]);
   }

#ifdef SHOWPROG
   // Just the progress bar stuff
   rhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
   tempPr=(double *) mxGetPr(rhs[0]);
   *tempPr=0;
   rhs[1]=mxCreateString("Calculating..");
   mexCallMATLAB(1, &lhs, 2, rhs, "waitbar");
   rhs[1]=lhs;
   progress=1/(double)(x_stop-x_start);
#endif

   for (j=0;j<binmax;j++){
      Histo1[j]=0; 
      Histo2[j]=0;
   }
   for (Pointx=x_start;Pointx<x_stop;Pointx++){
#ifdef SHOWPROG
      *tempPr+=progress;
      mexCallMATLAB(0, NULL, 2, rhs, "waitbar");
#endif
      for (Pointy=y_start;Pointy<y_stop;Pointy++){
         x=Pointx-hs;
         y=Pointy-hs;
         o=x+y*m;
         pROI_icon1D=pImage1D+o;
         pROI_icon2D=pImage2D+o;
         pROI_icon3D=pImage3D+o;
         sum=0;
         ind=0;
         for (scale=s_start;scale<=s_stop;scale++){
            int sind=scale-1;
            HistUint8_3DOpt4(pROI_icon1D, pROI_icon2D, pROI_icon3D, Histo1, pROIPix[scale], lengths[scale],nbins);
            sum+=lengths[scale];
            pEntropyArray[sind]=0;
            pDistArray[sind]=0;
            isum=1./sum; 
            //gettimeofday(&tv1, NULL);
            for (j=0;j<ind;j++)
               if (Histo1[pPNZ[j]]==0){ 	//this optimised so that it only does this when Histo2>0 (a bin is in PNZ)
                  double r=Histo2[pPNZ[j]];
                  r+=pDistArray[sind];
                  Histo2[pPNZ[j]]=0;
                  pDistArray[sind]=r;
               }
               ind=0;
               for (j=0;j<binmax;j++){
                  if (Histo1[j]!=0){
                     double r=Histo1[j];
                     pPNZ[ind++]=j;
                     tempNormBin=r/sum;
#ifdef FASTPLOG				
                     pEntropyArray[sind]-=PLOGP[(int)(tempNormBin*PLOGPRES)];
#else
                     pEntropyArray[sind]-=plogp(tempNormBin);
#endif
                     pDistArray[sind]+=fabs((tempNormBin)-(Histo2[j]));
                     Histo2[j]=tempNormBin;
                  }
               }
               //gettimeofday(&tv2, NULL);
               //tv3.tv_sec+=tv2.tv_sec-tv1.tv_sec;
         }
         counts=0;
         //Smooth the DistArray
         for (scale=s_start+1;scale<s_stop;scale++)
            pDistArray[scale]=(pDistArray[scale-1]+pDistArray[scale]+pDistArray[scale+1])/3;

         for (scale=s_start;scale<=s_stop-2;scale++){
            if ((pEntropyArray[scale]<pEntropyArray[scale+1]) &&
               (pEntropyArray[scale+1]>pEntropyArray[scale+2])){
                  peaks[counts]=scale+1;
                  counts++;
               }
         }
         if (NoOfScales+counts>=mxGetN(BestSaliency)){
            mxSetN(BestSaliency,NoOfScales+counts+100000);
            pBestSaliency=mxRealloc(pBestSaliency,6*(NoOfScales+counts+100000)*sizeof(double));
            mxSetPr(BestSaliency,pBestSaliency);
         }
         //assign the best (peaks) scales and global saliency for this x,y location
         if (counts>0)
            for (i=NoOfScales;i<NoOfScales+counts;i++){
               pBestSaliency[i*6+0]=Pointx;  /* the x location */
               pBestSaliency[i*6+1]=Pointy;  /* the y location */
               pBestSaliency[i*6+2]=peaks[i-NoOfScales]+1;	  /*the best scale*/
               Norm=pow(pBestSaliency[i*6+2],2)/(2*pBestSaliency[i*6+2]-1);
               pBestSaliency[i*6+3]=pEntropyArray[(peaks[i-NoOfScales])];	  /*HD*/
               pBestSaliency[i*6+4]=pDistArray[(peaks[i-NoOfScales])]*Norm;	  /*WD*/
               pBestSaliency[i*6+5]=Norm*pEntropyArray[(peaks[i-NoOfScales])]*(pDistArray[(peaks[i-NoOfScales])]); /* the saliency */
            }

            NoOfScales+=counts;
            for (j=0;j<binmax;j++){
               Histo1[j]=0;
               Histo2[j]=0;
            }
      }/*end of y loop*/
   }/*end of x loop*/
   mxSetN(BestSaliency,NoOfScales);
#ifdef SHOWPROG
   mexCallMATLAB(0, NULL, 1, &lhs, "close");
#endif
   free(peaks);
   free(ROIPix);
   free(lengths);
   free(Histo1);
   free(Histo2);
   free(pPNZ);
   return(BestSaliency);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   const mxArray *Image;
   int s_start, s_stop, Dims, i, nbins=0,AA=0;
   mxArray *BestSaliency;
   const int *Dims2;
   double sigma;

#ifdef PUBLICRELEASE
   PrintCopyright();
#endif
#ifdef FASTPLOG
   for (i=0;i<PLOGPRES;i++) PLOGP[i]=plogp((double)i/PLOGPRES);
#endif

   /*Validateinputs(nlhs, prhs, nrhs);*/

   /* get the input args */
   Image = prhs[0];
   s_start = *mxGetPr(prhs[1]);
   s_stop =  *mxGetPr(prhs[2]);
   nbins=*mxGetPr(prhs[3]);
   sigma= *mxGetPr(prhs[4]);
   AA= *mxGetPr(prhs[5]);

   Dims=mxGetNumberOfDimensions(Image);
   Dims2=mxGetDimensions(Image);

   if (Dims==2){
      mexPrintf("Input is scalar image. ");
      if (nbins){
         mexPrintf("Histogram with %d bins. ", nbins);
         mexPrintf( "Calculating Scale Saliency...");
         if (AA) {
            mexPrintf("using anti-aliased sampling...");
            BestSaliency=do_CalcSalScale1DAA(Image, s_start,s_stop, nbins);
         }
         else
            BestSaliency=do_CalcSalScale1D(Image, s_start,s_stop, nbins);	
         mexPrintf( "done.\n");
      }
      else
      {
         mexPrintf("Parzen window with sigma=%2.2f.\n", sigma);
         mexPrintf( "Calculating Scale Saliency...");
         BestSaliency=do_CalcSalScale1DParzen(Image, s_start,s_stop, sigma);
         mexPrintf( "done.\n");
      }
   }
   else if (Dims==3 && Dims2[2]==2){
      mexPrintf("Input is a 2D vector image. ");
      mexPrintf("Histogram with %d bins.\n", nbins);
      mexPrintf( "Calculating Scale Saliency...");
      BestSaliency=do_CalcSalScale2D(Image, s_start,s_stop, nbins);
      mexPrintf( "done.\n");
   }
   else if (Dims==3 && Dims2[2]==3){
      mexPrintf("Input is a 3D vector image. ");
      mexPrintf("Histogram with %d bins.\n", nbins);
      mexPrintf( "Calculating Scale Saliency...");
      BestSaliency=do_CalcSalScale3D(Image, s_start,s_stop, nbins);
      mexPrintf( "done.\n");
   }
   else {
      printf("Sorry cannot handle greater than 3 dimensional input images.\n");
      return;
   }

   plhs[0]=BestSaliency;
}

void ValidateInputs(int nlhs, const mxArray *prhs, int nrhs)
{
   /*Not yet implemented*/
}
