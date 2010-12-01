SRCFILES = DColorCluster.cpp DImageIO.cpp DLabelComponents.cpp DPoint.cpp DDiagMatrix.cpp DIntegralImage.cpp DMatrix.cpp DrawText.cpp DGaussianPyramid.cpp DistanceTransform.cpp DPlane.cpp DRect.cpp DImage.cpp DBoxMinTransform.cpp DMultiDMatrix.cpp

OBJECTS = DColorCluster.o DImageIO.o DLabelComponents.o DPoint.o DDiagMatrix.o DIntegralImage.o DMatrix.o DrawText.o DGaussianPyramid.o DistanceTransform.o DPlane.o DRect.o DImage.o DBoxMinTransform.o DMultiDMatrix.o

host-type = $(shell arch)

ifeq ($(host-type),ppc)
	CCOPTS = -O3
else
	CCOPTS = -funroll-loops -march=native -mfpmath=sse,387 -msse2
endif


.SUFFIXES: .o .cpp

.cpp.o:
	g++ -DGSL_SUPPORT -I . -o $@ -c $< -I ../gsl-1.5 $(CCOPTS) -O3 -funroll-loops -ffast-math -march=native -mfpmath=sse,387 -msse3
	g++ -DGSL_SUPPORT -I . -o p3/$@ -c $< -I ../gsl-1.5 -O3 -funroll-loops -march=native -mfpmath=sse,387 -msse
	g++ -DGSL_SUPPORT -I . -o db/$@ -c $< -I ../gsl-1.5 -g
	g++ -DGSL_SUPPORT -I . -o pg/$@ -c $< -I ../gsl-1.5 -g -pg

all: libDLib.a 


libDLib.a: $(OBJECTS)
	ar rs libDLib.a $(OBJECTS)
	ar rs libDLib_p3.a p3/*.o
	ar rs libDLib_db.a db/*.o
	ar rs libDLib_pg.a pg/*.o


clean:
	rm *.o *.a
