/*
Copyright 1998-2004 Timor Kadir.

Kadir/Brady Feature detector (Scale Saliency) Code 

For non-commericial use only.
*/

#include <math.h>
#include <string.h>
#include "mex.h"
#include "mat.h"
#ifndef _included_General
#include "General.h"
#endif

inline void HistUint8_1DOpt3Parz(uint8_T *a, int *y, uint32_T *pROIPix, int length);
inline void HistUint8_1DOpt4(uint8_T *a, int *y, uint32_T *pROIPix, int length);
inline void HistUint8_2DOpt4(uint8_T *a, uint8_T *b, int *y, uint32_T *pROIPix, int length, int nbins);
inline void HistUint8_3DOpt4(uint8_T *a, uint8_T *b, uint8_T *c, int *y, uint32_T *pROIPix, int length, int nbins);
inline void HistUint8_1DOpt4AA(uint8_T *a, int *y, int *pROIP, int *pROIW, int length);
