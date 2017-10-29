//tn_fnds v0.0.5   2012/3/17
//追加されているコメントには誤りがあるかもしれません。

// 音声分析合成法 WORLD by M. Morise
//
// FFTWを使うので，別途インストールが必要です．
//

// 現状で分かっているバグ
// decimateForF0 : 開始直後・終了間際4サンプルくらいに誤差が入ります．
//#include <fftsg.h>
//#include <fftw3.h>
#include <fft.h>

#include <stdlib.h>
#include <windows.h>
#include <math.h>

#define PI 3.1415926535897932384

// windowsならでは
#pragma warning( disable : 4996 )

#pragma comment(lib, "libfftw3-3.lib")
#pragma comment(lib, "libfftw3f-3.lib")
#pragma comment(lib, "libfftw3l-3.lib")

#define MAX_FFT_LENGTH 2048
#define FLOOR_F0 90.0//71.0    tn_fnds v0.0.4 低い方向へのF0誤検出防止。UTAUの原音はそんなに低くないだろうと楽観
#define DEFAULT_F0 500.0//150.0   tn_fnds v0.0.3 大きくしてみたら無声子音がきれいになった
#define LOW_LIMIT 65.0//EFB-GT

// 71は，fs: 44100においてFFT長を2048にできる下限．
// 70 Hzにすると4096点必要になる．
// DEFAULT_F0は，0.0.4での新機能．調整の余地はあるが，暫定的に決定する．

// F0推定法 DIO : Distributed Inline-filter Operation
void dio(double *x, int xLen, int fs, double framePeriod, 
		 double *timeAxis, double *f0);
int GetNumDIOSamples(int sample_rate, int num_samples, double frame_period);

// スペクトル包絡推定法 STAR : Synchronous Technique and Adroit Restoration
int getFFTLengthForStar(int fs);

//void star(double *x, int xLen, int fs, double *timeAxis, double *f0,
//		  double **specgram);
//void getMinimumPhaseSpectrum(double *inputSpec, fftw_complex *spectrum, fftw_complex *cepstrum, int fftl);

// 非周期性指標推定法 PLATINUM : 名称未定
void pt100(double *x, int xLen, int fs, double *timeAxis, double *f0,  
		 double **residualSpecgram);

//tn_fnds v0.0.3 にて追加
int pt101(double *x, int xLen, int fs, double *timeAxis, double *f0, 
		 double ***residualSpecgram, int **residualSpecgramLength, int *residualSpecgramIndex);

//tn_fnds v0.0.4 にて追加
void PulseResidualWindow(double **residualSpecgram, int *residualSpecgramLength, int pCount);

// WORLD Synthesis
//void synthesis(double *f0, int tLen, double **specgram, double **residualSpecgram, int fftl, double framePeriod, int fs, 
//			   double *synthesisOut, int xLen);
void synthesisPt100(double *f0, int tLen, double **residualSpecgram, int fftl, double framePeriod, int fs, 
			   double *synthesisOut, int xLen);
//void getMinimumPhaseSpectrum(double *inputSpec, fftw_complex *spectrum, fftw_complex *cepstrum, int fftl);

//tn_fnds v0.0.3 にて追加
void synthesisPt101(double fixedDefault_f0, double *f0, int tLen, double **aperiodicity, int *ResidualSpecgramLength,
					int *fixedResidualSpecgramIndex, double *volume,
					int fftl, double framePeriod, int fs, double *synthesisOut, int xLen);
//------------------------------------------------------------------------------------
// Matlab 関数の移植
double std2(double *x, int xLen);
void inv(double **r, int n, double **invr);
//void fftfilt(double *x, int xlen, double *h, int hlen, int fftl, double *y);
float randn(void);
void histc(double *x, int xLen, double *y, int yLen, int *index);
void interp1(double *t, double *y, int iLen, double *t1, int oLen, double *y1);
long decimateForF0(double *x, int xLen, double *y, int r);
void filterForDecimate(double *x, int xLen, double *y, int r);
int myround(double x);
void diff(double *x, int xLength, double *ans);
void interp1Q(double x, double shift, double *y, int xLength, double *xi, int xiLength, double *ans);