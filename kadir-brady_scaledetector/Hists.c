/*
Copyright 1998-2004 Timor Kadir.

Kadir/Brady Feature detector (Scale Saliency) Code 

For non-commericial use only.
*/


#include "Hists.h"

// All these rely on y being cleared if its a new histogram
// In this way we can incrementally add to bigger scales rather than start with totally new histograms 
// each time.
// 1 Dimensional histogramming

inline void HistUint8_1DOpt3Parz(uint8_T *a, int *y, uint32_T *pROIPix, int length) 
{ 
int i;

for (i = 0; i < length; i++)
	y[a[pROIPix[i]]]++;	
} 

inline void HistUint8_1DOpt4(uint8_T *a, int *y, uint32_T *pROIPix, int length) 
{ 
int i;

for (i = 0; i < length; i++)
	y[a[pROIPix[i]]]++;	
} 

inline void HistUint8_1DOpt4AA(uint8_T *a, int *y, int *pROIP, int *pROIW, int length) 
{ 
int i;

for (i = 0; i < length; i++)
	y[a[pROIP[i]]]+=pROIW[i];	
} 

inline void HistUint8_2DOpt4(uint8_T *a, uint8_T *b,int *y, uint32_T *pROIPix, int length, int nbins) 
{ 
int i;

for (i = 0; i < length; i++)
	y[a[pROIPix[i]] + b[pROIPix[i]]*nbins]++;	
} 

inline void HistUint8_3DOpt4(uint8_T *a, uint8_T *b, uint8_T *c, int *y, uint32_T *pROIPix, int length, int nbins) 
{ 
int i,o;

for (i = 0; i < length; i++){
  o=a[pROIPix[i]];
  o+=(b[pROIPix[i]]+c[pROIPix[i]]*nbins)*nbins;
  y[o]++;
  }
}

