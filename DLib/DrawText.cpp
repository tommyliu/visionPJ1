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
#include <math.h>
#include "DrawText.h"


/*  GIMP header image file format (RGB): /Users/dave/maya/test_imgs/letters2.h  */

static unsigned int width = 564;
static unsigned int height = 13;

/*  Call this macro repeatedly.  After each use, the pixel data can be extracted  */

#define HEADER_PIXEL(data,pixel) {\
  pixel[0] = (((data[0] - 33) << 2) | ((data[1] - 33) >> 4)); \
  pixel[1] = ((((data[1] - 33) & 0xF) << 4) | ((data[2] - 33) >> 2)); \
  pixel[2] = ((((data[2] - 33) & 0x3) << 6) | ((data[3] - 33))); \
  data += 4; \
}
static char *header_data =
	"````````````````````````````````````````````````````````````````"
	"````````````````!!!!````````````````````````````````````````````"
	"````````````````````````````````````````````````````````!!!!````"
	"````!!!!````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````!!!!````````````````````````````````````````````"
	"````````````````````````!!!!````````````````!!!!````!!!!````````"
	"````````!!!!````!!!!````````!!!!!!!!!!!!````````````!!!!````````"
	"````!!!!````!!!!!!!!````````````````````````!!!!````````````````"
	"````!!!!````````````````!!!!````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````!!!!"
	"!!!!!!!!````````````````!!!!````````````````!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!````````````````````!!!!````````!!!!!!!!!!!!!!!!"
	"!!!!````````````!!!!!!!!````````!!!!!!!!!!!!!!!!!!!!````````!!!!"
	"!!!!!!!!````````````!!!!!!!!!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````!!!!!!!!!!!!````````````````"
	"!!!!!!!!````````````````!!!!````````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!````````!!!!!!!!!!!!````````````!!!!!!!!!!!!!!!!"
	"!!!!````!!!!!!!!!!!!!!!!!!!!````````!!!!!!!!!!!!````````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!````````````!!!!!!!!!!!!!!!!````"
	"!!!!````````````!!!!````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````!!!!!!!!!!!!!!!!!!!!````"
	"````````!!!!!!!!!!!!````````!!!!````````````````````!!!!!!!!!!!!"
	"````````````````!!!!````````````````````````````````````````````"
	"````````````````````````````````````````!!!!````````````````````"
	"````````````````````````````````````````!!!!````````````````````"
	"````````````````!!!!!!!!````````````````````````````````!!!!````"
	"````````````````````````!!!!````````````````````````!!!!````````"
	"!!!!````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````!!!!!!!!````````````!!!!````````````!!!!!!!!````````"
	"````````````````````````````````````````!!!!````````````````!!!!"
	"````!!!!````````````````!!!!````!!!!````!!!!````!!!!````!!!!````"
	"!!!!````!!!!````````!!!!!!!!````````!!!!````````````````````!!!!"
	"````````````````!!!!````````````````````````!!!!````````!!!!````"
	"!!!!````!!!!````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"!!!!````!!!!````````````!!!!````````!!!!!!!!````````````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````````!!!!!!!!````````"
	"!!!!````````````````````````!!!!````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````````"
	"!!!!````````!!!!````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````!!!!````````````!!!!````!!!!````````!!!!````````"
	"!!!!````````````````````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````!!!!````````````````````"
	"````````!!!!````!!!!````````!!!!````````!!!!````````````````````"
	"!!!!!!!!````!!!!!!!!````!!!!!!!!````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````````!!!!````````````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````````"
	"````````!!!!````````````!!!!````````````````!!!!````````````````"
	"````````````!!!!````````````````!!!!````````````````````````````"
	"````````````````!!!!````````````````````````````````````!!!!````"
	"````````````````````````````````````````````````````````!!!!````"
	"````````````````````````````!!!!````````````````````````````````"
	"````````!!!!````````````````````````````````````````````````````"
	"````````````````!!!!````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````!!!!````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````!!!!````````````````````!!!!````````````"
	"````````!!!!````````````````````````````````````````````!!!!````"
	"````````````!!!!````!!!!````````!!!!!!!!!!!!!!!!!!!!!!!!!!!!````"
	"!!!!````````````!!!!````!!!!````!!!!````!!!!````````````````````"
	"````````````!!!!````````````````!!!!````````````````````````!!!!"
	"````````````!!!!!!!!!!!!````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````!!!!````````!!!!````````!!!!!!!!````!!!!````!!!!````"
	"````````````````````````!!!!````````````````````!!!!````````!!!!"
	"````!!!!````````!!!!````````````````````!!!!````````````````````"
	"````````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````````````!!!!!!!!````````````````!!!!!!!!````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````!!!!````!!!!````````!!!!!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````!!!!````````````````````!!!!````"
	"````````!!!!````!!!!````````````````````!!!!````````````````````"
	"!!!!````````````````````!!!!````````````!!!!````````````!!!!````"
	"````````````````````````!!!!````!!!!````!!!!````````````!!!!````"
	"````````````````!!!!!!!!````!!!!!!!!````!!!!!!!!````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````````````````````"
	"!!!!````````````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````!!!!````!!!!````````!!!!````!!!!````````````!!!!````!!!!"
	"````````````````````!!!!````````````````!!!!````````````````````"
	"!!!!````````````````````````!!!!````````````!!!!````!!!!````````"
	"````````````````````````````````````!!!!````````````!!!!!!!!!!!!"
	"!!!!````!!!!````!!!!!!!!````````````!!!!!!!!!!!!````````````!!!!"
	"!!!!!!!!!!!!````````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!!!!!!!!!````!!!!````!!!!!!!!````````````!!!!!!!!````"
	"````````````!!!!!!!!!!!!````````!!!!````````!!!!````````````````"
	"!!!!````````````!!!!!!!!````!!!!````````!!!!````!!!!!!!!````````"
	"````!!!!!!!!!!!!````````!!!!````!!!!!!!!````````````!!!!!!!!!!!!"
	"!!!!````!!!!!!!!````!!!!!!!!````````!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!!!!!!!!!!!!!!!!!````````````!!!!````````````````````"
	"!!!!````````````````````!!!!````````````````````````````````````"
	"````````!!!!````````````````````````````````````````````!!!!````"
	"!!!!````````!!!!!!!!````````````````!!!!````!!!!````````````!!!!"
	"````````````````````````````````````````````!!!!````````````````"
	"````````````````!!!!````````!!!!!!!!!!!!````````````````!!!!````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````!!!!````````!!!!````!!!!````!!!!````"
	"````````!!!!````````````````````````!!!!````````````````!!!!!!!!"
	"````````!!!!````````!!!!````````!!!!!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````````````````!!!!````````````!!!!!!!!!!!!````````"
	"!!!!````````````!!!!````````````!!!!!!!!````````````````!!!!!!!!"
	"````````````````````````!!!!!!!!````````````````````````!!!!!!!!"
	"````````````````````````````!!!!````````!!!!````!!!!````!!!!````"
	"````!!!!````!!!!````````!!!!!!!!!!!!!!!!````````!!!!````````````"
	"````````!!!!````````````!!!!````!!!!!!!!!!!!!!!!````````!!!!!!!!"
	"!!!!!!!!````````!!!!````````````````````!!!!!!!!!!!!!!!!!!!!````"
	"````````!!!!````````````````````````````!!!!````!!!!!!!!````````"
	"````````!!!!````````````````````!!!!````!!!!````!!!!````!!!!````"
	"!!!!````!!!!````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````````!!!!!!!!!!!!"
	"````````````````!!!!````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````!!!!````!!!!````````````!!!!````````````"
	"````!!!!````!!!!````````````````!!!!````````````````````!!!!````"
	"````````````````!!!!````````````````````````!!!!````````````!!!!"
	"````!!!!````````````````````````````````````````````````````````"
	"!!!!````````````!!!!````!!!!!!!!````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````!!!!"
	"````````````````!!!!````````````!!!!````!!!!!!!!````````!!!!````"
	"````````!!!!````````````````````````!!!!````````!!!!````!!!!````"
	"````````````````!!!!````````````!!!!````!!!!````!!!!````!!!!!!!!"
	"````````!!!!````!!!!````````````!!!!````!!!!!!!!````````!!!!````"
	"!!!!````````````!!!!````````!!!!!!!!````````````!!!!````````````"
	"!!!!````````!!!!````````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````!!!!````!!!!````````"
	"!!!!````````````!!!!````````````````!!!!````````````````!!!!````"
	"````````````````!!!!````````````````````!!!!````````````````!!!!"
	"!!!!````````!!!!````````!!!!````````````````````````````````````"
	"````````!!!!````!!!!````````````!!!!!!!!````````````````!!!!````"
	"!!!!````!!!!````!!!!````````````````````````````````````````!!!!"
	"````````````````````````````````!!!!````!!!!````!!!!````!!!!````"
	"````````!!!!````````````````````````````````````````````````````"
	"````````````````````````````````````````!!!!````````````!!!!!!!!"
	"````````!!!!````````````!!!!````````````````````!!!!````````````"
	"````````````````!!!!````!!!!!!!!!!!!!!!!!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````````!!!!````````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!!!!!````````````````````````````"
	"````````````````````````````````!!!!!!!!````````!!!!!!!!!!!!!!!!"
	"!!!!````````````!!!!!!!!````````````````!!!!````````````!!!!````"
	"!!!!````!!!!````````!!!!````!!!!````````!!!!````````````!!!!````"
	"!!!!````````````````````!!!!````````````!!!!````!!!!````````````"
	"````````!!!!````````````````````!!!!````````!!!!!!!!````!!!!````"
	"````````!!!!````````````!!!!````````````````````````````!!!!````"
	"!!!!!!!!````````````````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````!!!!````!!!!````!!!!````````````!!!!````!!!!!!!!"
	"!!!!!!!!````````!!!!````````````!!!!````!!!!!!!!!!!!!!!!````````"
	"````````````````!!!!````````````!!!!````````````!!!!````````````"
	"!!!!````````!!!!````!!!!````````!!!!````!!!!````!!!!````````````"
	"!!!!````````````````````!!!!````````````````````!!!!````````````"
	"````````!!!!````````````````````````!!!!````````````````````!!!!"
	"````````!!!!````````````!!!!````````````````````````````````````"
	"````````````````!!!!````````````!!!!````!!!!````````````!!!!````"
	"!!!!````````````````````!!!!````````````!!!!````!!!!!!!!!!!!!!!!"
	"!!!!````````!!!!````````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````````````!!!!````````````````````````!!!!````````"
	"!!!!!!!!````````````````````````!!!!````````````!!!!````!!!!````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````!!!!````````!!!!````````````````"
	"````!!!!!!!!````````````````!!!!````````````````!!!!````````````"
	"!!!!````````!!!!````!!!!````````!!!!````!!!!````!!!!````````````"
	"!!!!````````````````!!!!````!!!!````````````````!!!!````````````"
	"````````!!!!````````````````````!!!!````````````````````!!!!````"
	"````````!!!!````````!!!!!!!!````````````!!!!````````````````````"
	"````````````````!!!!!!!!!!!!!!!!!!!!!!!!````````!!!!````!!!!````"
	"````!!!!````!!!!````!!!!!!!!````````!!!!````!!!!````````````````"
	"````````````!!!!````````````````````````````````!!!!````````````"
	"!!!!````````````!!!!!!!!!!!!!!!!!!!!````````````````````````````"
	"!!!!!!!!!!!!!!!!!!!!````````````````````````````````!!!!````````"
	"````````!!!!````````````!!!!````````````!!!!````````````````!!!!"
	"````````````````````````````````!!!!````````````````!!!!````````"
	"````````````````!!!!````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````````````````````!!!!````````````"
	"````````````````````````````````````````!!!!!!!!````````````````"
	"````````````````````````````````````````!!!!!!!!````````!!!!````"
	"````````!!!!````````!!!!!!!!````!!!!!!!!!!!!!!!!!!!!````!!!!````"
	"````````!!!!````!!!!````````````````````!!!!````````````!!!!````"
	"!!!!````````````````````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````!!!!````!!!!````````````!!!!````````````````````"
	"!!!!````````````!!!!````!!!!````````!!!!!!!!````!!!!````````````"
	"!!!!````!!!!````````````````````!!!!````````````!!!!````!!!!````"
	"!!!!````````````````````````````!!!!````````````!!!!````````````"
	"!!!!````````````!!!!````````!!!!````!!!!````````````!!!!````!!!!"
	"````````````!!!!````!!!!````````````````!!!!````````````````!!!!"
	"````````````````````````!!!!````````````````````````!!!!````````"
	"````````````!!!!````````````````````````````````````````````````"
	"````````````````````````````````!!!!````````````!!!!````!!!!````"
	"````````!!!!````!!!!````````````````````!!!!````````````!!!!````"
	"!!!!````````````````````````!!!!````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````````````!!!!````````````````````"
	"````!!!!````````!!!!````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````!!!!````!!!!````````````!!!!````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````````!!!!"
	"````````````````````````````!!!!````````````!!!!````````````````"
	"!!!!````````````!!!!````````!!!!````!!!!````````!!!!````!!!!````"
	"!!!!````````````!!!!````````````````!!!!````!!!!````````````!!!!"
	"````````````````!!!!!!!!````````````````````````!!!!````````````"
	"````````````!!!!!!!!````````````````````````````````````````````"
	"````````````````````````````````````!!!!````!!!!````````!!!!````"
	"!!!!````!!!!````!!!!````````!!!!````!!!!!!!!````````````!!!!````"
	"````````````````````````````!!!!````````````````````````````````"
	"!!!!````````````````````````````````````!!!!````````````````````"
	"!!!!!!!!````````````````````````````````````````!!!!!!!!````````"
	"````!!!!````````````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````````````!!!!````````````!!!!````````````"
	"````!!!!````````````````````````!!!!````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````````````!!!!````````````````!!!!"
	"````````````````!!!!!!!!````````````````!!!!!!!!````````````````"
	"!!!!!!!!````````!!!!!!!!!!!!!!!!!!!!````````````!!!!!!!!````````"
	"````````````````````````````!!!!````````````````!!!!````````````"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````!!!!````````!!!!````````````````````!!!!````````````````````"
	"!!!!````````````!!!!````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!````!!!!````````!!!!````````!!!!````"
	"````````````````!!!!````````````!!!!````!!!!````````!!!!!!!!````"
	"!!!!````````````!!!!````!!!!````````````````````!!!!````````````"
	"!!!!````!!!!````````!!!!````````!!!!````````````!!!!````````````"
	"!!!!````````````!!!!````````````!!!!````````````!!!!````````````"
	"````!!!!````!!!!````````!!!!````````````!!!!````````````!!!!````"
	"````````!!!!````````````````````````````!!!!````````````````````"
	"````````!!!!````````````````!!!!````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````!!!!"
	"!!!!````!!!!````````````!!!!````!!!!````````````!!!!````!!!!````"
	"````!!!!!!!!````!!!!````````````!!!!````````!!!!````````````````"
	"!!!!````````!!!!!!!!````!!!!````````````!!!!````````````!!!!````"
	"````````````````````!!!!````````!!!!````````!!!!````````````````"
	"!!!!````````````!!!!````!!!!````!!!!````!!!!````````````!!!!````"
	"!!!!````````````!!!!````!!!!````````````!!!!````!!!!````````!!!!"
	"!!!!````````!!!!````````````````!!!!````````````!!!!````````!!!!"
	"````````!!!!````!!!!````````!!!!!!!!````````````!!!!````````````"
	"````!!!!````!!!!````````````!!!!````!!!!````````````````!!!!````"
	"````````!!!!````````````````````````````!!!!````````````````````"
	"!!!!````````````````````!!!!````````````````````````````````````"
	"````````!!!!````````````````````````````````````````!!!!````!!!!"
	"````````````!!!!!!!!!!!!````````!!!!````````````!!!!````````!!!!"
	"!!!!!!!!````!!!!````````````````````````````````!!!!````````````"
	"````````````!!!!````````````````````````````````````````!!!!````"
	"````````````````!!!!!!!!````````````````````````````````````````"
	"!!!!!!!!````````!!!!````````````````````````!!!!!!!!!!!!````````"
	"!!!!!!!!!!!!!!!!!!!!````!!!!!!!!!!!!!!!!!!!!````````!!!!!!!!!!!!"
	"````````````````````!!!!````````!!!!!!!!!!!!!!!!````````````!!!!"
	"!!!!!!!!````````````````!!!!````````````````!!!!!!!!!!!!````````"
	"````!!!!!!!!````````````````````!!!!!!!!````````````````!!!!!!!!"
	"````````````````````````!!!!!!!!````````````````````````!!!!!!!!"
	"````````````````````````!!!!````````````````````!!!!!!!!!!!!````"
	"!!!!````````````!!!!````!!!!!!!!!!!!!!!!````````````!!!!!!!!!!!!"
	"````````!!!!!!!!!!!!````````````!!!!!!!!!!!!!!!!!!!!````!!!!````"
	"````````````````````!!!!!!!!!!!!````````!!!!````````````!!!!````"
	"````!!!!!!!!!!!!````````````!!!!!!!!!!!!````````!!!!````````````"
	"!!!!````!!!!!!!!!!!!!!!!!!!!````!!!!````````````!!!!````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!````````!!!!````````````````````"
	"````!!!!!!!!!!!!````````!!!!````````````!!!!````````!!!!!!!!!!!!"
	"````````````````!!!!````````````````!!!!!!!!!!!!````````````````"
	"!!!!````````````````!!!!````!!!!````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!!!!!!!!!!!!!!!!!````````````!!!!````"
	"````````````````````````!!!!````````````````!!!!````````````````"
	"````````````````````````````````````````````````````````````````"
	"````!!!!!!!!````!!!!````!!!!!!!!!!!!!!!!````````````!!!!!!!!!!!!"
	"````````````!!!!!!!!````!!!!````````!!!!!!!!!!!!````````````!!!!"
	"````````````````````!!!!!!!!````!!!!````!!!!````````````!!!!````"
	"````````!!!!!!!!````````````````````!!!!````````!!!!````````````"
	"!!!!````````````!!!!!!!!````````!!!!````!!!!````!!!!````!!!!````"
	"````````!!!!````````!!!!!!!!!!!!````````!!!!!!!!!!!!!!!!````````"
	"````!!!!!!!!````!!!!````!!!!!!!!!!!!````````````````!!!!!!!!!!!!"
	"````````````````!!!!!!!!````````````!!!!!!!!````!!!!````````````"
	"!!!!````````````````!!!!````!!!!````````!!!!````````````!!!!````"
	"````````!!!!````````````!!!!!!!!!!!!!!!!!!!!````````````!!!!````"
	"````````````````!!!!````````````````````!!!!````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"!!!!````````````````````````!!!!````````````````````````````````"
	"````````````````````````````````````!!!!````````````````````````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````!!!!````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````!!!!````````````````````````````````````````````````!!!!"
	"````````````````````````````````!!!!!!!!!!!!!!!!!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````!!!!````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````!!!!````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````!!!!````````````````````````````````````````"
	"````````!!!!````````````````````!!!!````````````````````!!!!````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````!!!!````````````````!!!!````````````````````"
	"````````````````````````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````!!!!````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````!!!!!!!!````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````!!!!````````````````````````````````````"
	"````````````!!!!````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````!!!!!!!!!!!!"
	"````````````````````````````````````````````````````````!!!!!!!!"
	"!!!!````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````!!!!````````````````````````````````````!!!!````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````!!!!````````````````````````````"
	"````````````````````````!!!!````````````````````````````````````"
	"````````!!!!````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````!!!!````````!!!!````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````!!!!!!!!!!!!````````````"
	"````````````````````!!!!!!!!!!!!````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````!!!!!!!!````````````"
	"````````````````!!!!!!!!````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````````````````````````````````````````````````````"
	"````````````````";



static const int numbers[10][120] = {
    {0, 0, 0, 1, 1, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 0, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 0, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 1, 1, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 0, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 0, 0, 0,
     0, 1, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 0, 1, 1, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 0, 0, 0, 1, 0, 0,
     0, 1, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 0, 0, 0, 1, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 1, 1, 1, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0},

    {0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 0, 0, 0, 0, 1, 0,
     0, 1, 1, 1, 1, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0}};



static const int chr_rows=13;
static const int chr_cols=6;

void draw_number(DPlane &img, int number, int digits, int row, int col,
	    int value)
{
    for(int d=digits; d>0; d--)
    {
        int num = (number % int(pow(10,d))) / int(pow(10,d-1));

        for(int i=0; i<chr_rows; i++)
            for(int j=0; j<chr_cols; j++)
                if(row+i >= 0 && col+j >=0 && row+i < img.rows() && col+j < img.cols())
                    img[row+i][col+j] = int(img[row+i][col+j])/2 | int(numbers[num][i*chr_cols + j] * value);

        col += chr_cols;
    }

    return;
}


static int first=1;

void draw_text(DPlane &img, char *str, int row, int col, int value)
{
  static DPlane characters(height, width);
  if(first)
    {
      char *data = header_data;
      for(int i=0; i<height; i++)
	for(int j=0; j<width; j++)
	  {
	    char pixel[4];

	    HEADER_PIXEL(data, pixel);

	    characters[i][j] = pixel[0];
	  }

      first=0;
    }

  DPoint pt(row, col);

  for(int i=0; i<strlen(str); i++)
    {
      int char_no = str[i] - 33;

      img.set_submatrix(pt, characters.extract(DRect(0, char_no*chr_cols, chr_rows-1, (char_no+1)*chr_cols-1)));
      pt.col(pt.col() + chr_cols);
    }
}
