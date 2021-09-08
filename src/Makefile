OBJS = dio.o frqread.o matlabfunctions.o platinum.o star.o synthesis.o macresmain.o fft.o wavread.o
CXX = g++
CFLAGS = -O -msse3
INCLUDES = -I.
LIBS = 

macres : $(OBJS) fft.h frqread.h mmsystem.h windows.h wavread.h world.h
	$(CXX) -o $@ $(OBJS) $(LIBS) $(CFLAGS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm macres $(OBJS)
