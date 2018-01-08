#include "world.h"

#include <stdio.h> // for debug
#include <stdlib.h>
#include <math.h>


void tandemSTRAIGHTGeneralBody(double *x, int xLen, int fs, double f0, double t, double f0LowLimit, 
							   double q1, double exponentControl, double correctionForBlackman, int fftl,
							   double * sliceSTRAIGHT);

int getFFTLengthForTandemStraight(int fs)
{
	double correctionForBlackman = 2.5;
	double f0LowLimit = 60.0;
	return (int)pow(2.0, 1.0+(int)(log(correctionForBlackman*fs/f0LowLimit+1) / log(2.0)));
}

// Spectrum extraction based on TANDEM-STRAIGHT
// x	: Input signal whose length is xLen [sample]
// xLen : Length of the input signal.
// f0	: Estimated F0 contour
void tandemStraight(double *x, int xLen, int fs, double *timeAxis, double *f0,
		 double **specgram)
{
	int i,j;
	double framePeriod = (timeAxis[1]-timeAxis[0])*1000.0;
	double f0LowLimit = 60.0; // F0が60 Hz以下の場合は無声音として扱う
	double defaultF0 = 300;
	double q1 = -0.315669; // よく分からん
	double exponentControl = 1.0/6.0;
	double correctionForBlackman = 2.5; // ここを変更してはいけない．

	int	fftl = (int)pow(2.0, 1.0+(int)(log(correctionForBlackman*fs/f0LowLimit+1) / log(2.0)));
	int tLen = getSamplesForDIO(fs, xLen, framePeriod);

	double currentTime;
	double currentF0;
	double *sliceSTRAIGHT;
	sliceSTRAIGHT = (double *)malloc(sizeof(double) * fftl);
	for(i = 0;i < tLen;i++)
	{
		currentTime = timeAxis[i];
		currentF0 = f0[i] < f0LowLimit ? defaultF0 : f0[i]; // 将来的に60 Hz以下だったらに変える．
		tandemSTRAIGHTGeneralBody(x, xLen, fs, currentF0, currentTime, f0LowLimit, 
			q1, exponentControl, correctionForBlackman, fftl, 
			sliceSTRAIGHT);
		for(j = 0;j <= fftl/2;j++)
			specgram[i][j] = sliceSTRAIGHT[j];
	}

	free(sliceSTRAIGHT);


}

// matlabに順ずる丸め
static int myround2(double x)
{
	if(x > 0)
		return (int)(x+0.5);
	else
		return (int)(x-0.5);
}

void cumsum(double *x, int xLength, double *ans)
{
	ans[0] = x[0];
	for(int i = 1;i < xLength;i++)
	{
		ans[i] = ans[i-1] + x[i];
	}
}

void diff(double *x, int xLength, double *ans)
{
	for(int i = 0;i < xLength-1;i++)
	{
		ans[i] = x[i+1] - x[i];
	}
	return;
}


// 河原先生作成のMatlab関数interp1Hを移植．
// 基本的には同じだが，配列の要素数を明示的に指定する必要がある．
void interp1Q(double x, double shift, double *y, int xLength, double *xi, int xiLength, double *ans)
{
	double deltaX;
	double *xiFraction, *deltaY;
	int *xiBase;
	int i;

	xiFraction = (double *)malloc(xiLength*sizeof(double));
	deltaY = (double *)malloc(xLength*sizeof(double));
	xiBase = (int *)malloc(xiLength*sizeof(int));

	deltaX = shift;
	for(i = 0;i < xiLength;i++)
	{
		xiBase[i] = (int)floor((xi[i] - x) / deltaX);
		xiFraction[i] = (double)(xi[i]-x)/deltaX - (double)xiBase[i];
	}
	diff(y, xLength, deltaY);
	deltaY[xLength-1] = 0.0;

	for(i = 0;i < xiLength;i++)
	{
		ans[i] = y[xiBase[i]] + deltaY[xiBase[i]]*xiFraction[i]; 
	}

	free(xiFraction);
	free(xiBase);
	free(deltaY);
}

void tandemSTRAIGHTGeneralBody(double *x, int xLen, int fs, double f0, double t, double f0LowLimit, 
							   double q1, double exponentControl, double correctionForBlackman, int fftl,
							   double * sliceSTRAIGHT)
{
	int i,j;
	double t0 = 1.0 / f0;

	int *baseIndex, *preIndex, *postIndex; // i付きのも包含する (Matlab版参照)
	int nFragment = (int)(0.5 + correctionForBlackman*(double)fs/f0/2.0);

	baseIndex = (int *)malloc(sizeof(int) * (nFragment*2+1));
	preIndex  = (int *)malloc(sizeof(int) * (nFragment*2+1));
	postIndex = (int *)malloc(sizeof(int) * (nFragment*2+1));

	for(i = -nFragment, j = 0;i <= nFragment;i++, j++)
		baseIndex[j] = i;

	for(i = 0;i <= nFragment*2;i++)
	{
		preIndex[i]  = min(xLen, max(1, myround2((t-t0/4.0)*(double)fs+1+baseIndex[i]) ) ) - 1;
		postIndex[i] = min(xLen, max(1, myround2((t+t0/4.0)*(double)fs+1+baseIndex[i]) ) ) - 1;
	}

	double *preSegment, *postSegment;
	double *preWindow, *postWindow;
	double preTime, postTime, preAverage, postAverage;
	preSegment  = (double *)malloc(sizeof(double) * (nFragment*2+1));
	postSegment = (double *)malloc(sizeof(double) * (nFragment*2+1));
	preWindow   = (double *)malloc(sizeof(double) * (nFragment*2+1));
	postWindow  = (double *)malloc(sizeof(double) * (nFragment*2+1));
	preAverage  = postAverage = 0.0;
	for(i = 0;i <= nFragment*2;i++)
	{
		preSegment[i]  = x[preIndex[i]];
		postSegment[i] = x[postIndex[i]];
		preTime  = (double)baseIndex[i]/(double)fs/(correctionForBlackman/2.0) + 
			((t-t0/4.0)*(double)fs - (double)(myround2((t-t0/4.0)*(double)fs))) / (double)fs;
		postTime = (double)baseIndex[i]/(double)fs/(correctionForBlackman/2.0) + 
			((t+t0/4.0)*(double)fs - (double)(myround2((t+t0/4.0)*(double)fs))) / (double)fs;

		preWindow[i]  = 0.5*cos(PI*preTime*f0) +0.42+0.08*cos(2.0*PI*preTime*f0);
		postWindow[i] = 0.5*cos(PI*postTime*f0)+0.42+0.08*cos(2.0*PI*postTime*f0);
		preAverage  += preWindow[i]*preWindow[i];
		postAverage += postWindow[i]*postWindow[i];
	}
	preAverage  = sqrt(preAverage);
	postAverage = sqrt(postAverage);
	for(i = 0;i <= nFragment*2;i++)
	{
		preWindow[i]  /= preAverage;
		postWindow[i] /= postAverage;
	}

	// 波形のスペクトルを計算
	double				*waveform1, *waveform2;
	double				*tandemSpec;
	waveform1  = (double *)malloc(sizeof(double) * fftl);
	waveform2  = (double *)malloc(sizeof(double) * fftl);
	tandemSpec = (double *)malloc(sizeof(double) * fftl);

	fftw_plan			forwardFFT1;				// FFTセット
	fftw_plan			forwardFFT2;				// FFTセット
	fftw_complex		*ySpec1, *ySpec2;	// スペクトル
	ySpec1 = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);
	ySpec2 = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);
	forwardFFT1 = fftw_plan_dft_r2c_1d(fftl, waveform1, ySpec1, FFTW_ESTIMATE);
	forwardFFT2 = fftw_plan_dft_r2c_1d(fftl, waveform2, ySpec2, FFTW_ESTIMATE);

	double dc1, dc2; // DCは直流成分
	dc1 = dc2 = 0.0;
	for(i = 0;i <= nFragment*2;i++)
	{
		dc1 += preSegment[i];
		dc2 += postSegment[i];
	}
	dc1 /= (double)(nFragment*2+1);
	dc2 /= (double)(nFragment*2+1);
	// パワースペクトルの計算
	for(i = 0;i <= nFragment*2;i++)
	{
		waveform1[i] = (preSegment[i] -dc1) * preWindow[i];
		waveform2[i] = (postSegment[i]-dc2) * postWindow[i];
	}
	for(;i < fftl;i++)
	{
		waveform1[i] = 0.0;
		waveform2[i] = 0.0;
	}

	fftw_execute(forwardFFT1); // FFTの実行
	fftw_execute(forwardFFT2); // FFTの実行

	for(i = 0;i <= fftl/2;i++)
	{
		tandemSpec[i] = ySpec1[i][0]*ySpec1[i][0] + ySpec1[i][1]*ySpec1[i][1]
		              + ySpec2[i][0]*ySpec2[i][0] + ySpec2[i][1]*ySpec2[i][1];
	}
	free(ySpec1); free(ySpec2);
	fftw_destroy_plan(forwardFFT1);
	fftw_destroy_plan(forwardFFT2);
	free(preSegment); free(postSegment);
	free(preWindow); free(postWindow);
	free(baseIndex); free(preIndex); free(postIndex);

	// tandem2straight
	double *dSpectrum, dFrequencyAxis, dShift;
	dSpectrum		= (double *)malloc(sizeof(double) * fftl*2);
	double *shaper;
	shaper = (double *)malloc(sizeof(double) * ((int)(f0*fftl/(double)fs)*2 + 1) );

	// 計算コストを少しでも減らす
	dFrequencyAxis = -(double)fs;
	dShift = (double)fs/(double)fftl; 
	dSpectrum[0] = dSpectrum[fftl] = tandemSpec[0];
	for(i = 1;i < fftl/2;i++)
	{
		dSpectrum[i] = dSpectrum[i+fftl] = tandemSpec[i];
		dSpectrum[fftl/2+i] = dSpectrum[fftl/2+i+fftl] = tandemSpec[fftl/2-i];
	}
	dSpectrum[fftl/2] = dSpectrum[fftl/2+fftl] = tandemSpec[fftl/2];

	int tmp = (int)(f0*fftl/(double)fs);
	for(i = 0;i <= tmp*2;i++)
		shaper[i] = 0.5 + 0.5*cos(PI*(i-tmp)/(f0/(double)fs*(double)fftl));

	double *xi, *filler;
	xi		= (double *)malloc(sizeof(double) * ((int)(f0*fftl/(double)fs)*2 + 1) );
	filler  = (double *)malloc(sizeof(double) * ((int)(f0*fftl/(double)fs)*2 + 1) );
	tmp = (int)(f0*fftl/(double)fs);
	for(i = 0;i < ((int)(f0*fftl/(double)fs)*2 + 1);i++)
		xi[i] = f0 - fabs((double)i-tmp)/(double)fftl*(double)fs; 
	interp1Q(dFrequencyAxis, dShift, dSpectrum, fftl*2, xi, ((int)(f0*fftl/(double)fs)*2 + 1), filler);

	for(i = 0;i < ((int)(f0*fftl/(double)fs)*2 + 1);i++)
		dSpectrum[fftl+i-tmp] = dSpectrum[fftl+i-tmp]*(1-shaper[i]) + filler[i]*shaper[i];

	// ok Dec.25
	double parameter = 0.001;
	double *dSegment, *centers;
	dSegment	= (double *)malloc(sizeof(double) * fftl*2);
	centers			= (double *)malloc(sizeof(double) * (fftl/2 + 1) );

	dSegment[0]		= pow(dSpectrum[0], parameter)*(double)fs/(double)fftl;
	for(i = 1;i < fftl*2;i++)
		dSegment[i] = pow(dSpectrum[i], parameter)*(double)fs/(double)fftl + dSegment[i-1];

	for(i = 0;i <= fftl/2;i++)
		centers[i] = (double)i / (double)fftl * (double)fs - f0/2.0;

	double *lowLevels, *highLevels;
	lowLevels  = (double *)malloc(sizeof(double) * (fftl/2+1));
	highLevels = (double *)malloc(sizeof(double) * (fftl/2+1));
	interp1Q(dFrequencyAxis + (double)fs/(double)fftl/2.0, dShift, dSegment, fftl*2, centers, (fftl/2 + 1), lowLevels);
	for(i = 0;i <= fftl/2;i++)
		centers[i] += f0;
	interp1Q(dFrequencyAxis + (double)fs/(double)fftl/2.0, dShift, dSegment, fftl*2, centers, (fftl/2 + 1), highLevels);

	for(i = 0;i <= fftl/2;i++)
	{
		sliceSTRAIGHT[i] = max(0.00000000000001, 
			pow( (highLevels[i]-lowLevels[i])/f0, 1.0/parameter) );
	}

	free(lowLevels); free(highLevels);
	free(dSegment); free(centers);
	free(xi); free(filler);
	free(shaper);
	free(dSpectrum);
	free(waveform1); free(waveform2); 
	free(tandemSpec); 
}