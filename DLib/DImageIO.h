//
//  DLib: A simple image processing library.
//
//  David Crandall, 2003-2005
//  crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
//
//
#include "DImage.h"

DImage ReadPPMImage(const char *filename);
DImage ReadPBMImage(const char *filename);
void WritePPMImage(const DImage &img, const char *filename);
void WritePPMImage(const DPlane &img, const char *filename);
