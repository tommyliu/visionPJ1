
match_lite : match_lite.cpp Appear.cpp KFan.cpp Train_Appear.cpp Train_KFan.cpp
	g++ match_lite.cpp Appear.cpp KFan.cpp Train_Appear.cpp Train_KFan.cpp -IDLib -L. -Lgsl-1.5/.libs -Lgsl-1.5/cblas/.libs -lDLib -lgsl -lgslcblas -Igsl-1.5 -LDLib -o match_lite -I.  -O3 -funroll-loops -march=native -mfpmath=sse,387 -ffast-math  -msse3 -DGSL_SUPPORT


