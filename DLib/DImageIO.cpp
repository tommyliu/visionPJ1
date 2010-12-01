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
#include "stdio.h"
#include <fstream>
#include <iostream>
#include "DImage.h"
#include <assert.h>

using namespace std;

DImage ReadPPMImage(const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  char temp[1024];
  assert(fp);
  fgets(temp, 1024, fp);
  assert(temp[0] == 'P' && temp[1] == '6');

  int width, height, maxval;
  do {
    fgets(temp, 1024, fp);
  } while(temp[0] == '#');

  sscanf(temp, "%d %d", &width, &height);

  fgets(temp, 1024, fp);
  sscanf(temp, "%d", &maxval);

  DImage image(3, height, width);

  for(int i=0; i<height; i++)
    for(int j=0; j<width; j++)
      for(int k=0; k<3; k++)
      {
        image[k][i][j] = fgetc(fp);
      } 

  fclose(fp);

  return image;
}


DImage ReadPBMImage(const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  char temp[1024];

  fgets(temp, 1024, fp);
  assert(temp[0] == 'P' && temp[1] == '4');

  int width, height, maxval;
  do {
    fgets(temp, 1024, fp);
  } while(temp[0] == '#' || temp[0]=='\n');

  sscanf(temp, "%d %d", &width, &height);

  DImage image(1, height, width);

  for(int i=0; i<height; i++)
    for(int j=0; j<ceil(width/8.0); j++)
      {
         int tmp = fgetc(fp);
         for(int k=128, l=0; k>=1 && l<8; k=k>>1, l++)
         {
           if(j*8+l < width)
              image[0][i][j*8+l] = (tmp & k)?255:0;
         }
      } 

  fclose(fp);

  return image;
}


void WritePPMImage(const DImage &img, const char *filename)
{
  FILE *fp = fopen(filename, "wb");
  char temp[1024];

  // write magic number
  fprintf(fp, "P6\n");

  // write dimensions
  fprintf(fp, "%d %d\n", img.cols(), img.rows());
  
  // write max pixel value
  fprintf(fp, "255\n");

  for(int i=0; i<img.rows(); i++)
    for(int j=0; j<img.cols(); j++)
      for(int k=0; k<3; k++)
	fputc(int(img[k][i][j]), fp);

  fclose(fp);

  return;
}

void WritePPMImage(const DPlane &img, const char *filename)
{
  FILE *fp = fopen(filename, "wb");
  char temp[1024];

  // write magic number
  fprintf(fp, "P6\n");

  // write dimensions
  fprintf(fp, "%d %d\n", img.cols(), img.rows());
  
  // write max pixel value
  fprintf(fp, "255\n");

  for(int i=0; i<img.rows(); i++)
    for(int j=0; j<img.cols(); j++)
      for(int k=0; k<3; k++)
	fputc(int(img[i][j]), fp);

  fclose(fp);

  return;
}
