/*
Copyright 1998-2004 Timor Kadir.

Kadir/Brady Feature detector (Scale Saliency) Code 

For non-commericial use only.
*/

#define _included_General

#ifdef _WIN32
#define inline 
#define rint(x) (( (x-floor(x))>=0.5)?(ceil(x)):(floor(x)))
#endif
#define MAX(x,y) ((x)<(y)?(y):(x)) 
#define MIN(x,y) ((x)>(y)?(y):(x)) 
