export FP_TYPE ?= float

LIB_PREFIX=../..
OBJS = dio.o frqread.o matlabfunctions.o platinum.o star.o synthesis.o macresmain.o fft.o wavread.o
CXX = g++
CFLAGS = -O -msse3 -DFP_TYPE=$(FP_TYPE)
INCLUDES = -I.  -I$(LIB_PREFIX)/include
LIBS = $(LIB_PREFIX)/lib/libpyin.a $(LIB_PREFIX)/lib/libgvps.a

macres : $(OBJS) $(LIBS) fft.h frqread.h mmsystem.h windows.h wavread.h world.h
	$(CXX) -o $@ $(OBJS) $(LIBS) $(CFLAGS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm macres $(OBJS)
