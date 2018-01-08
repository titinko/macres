// MacRes v0.0.1   1/7/18
#include "world.h"

#include <stdio.h> // for debug
#include <stdlib.h>

// spectrum, cepstrumは毎回malloc, freeするのが面倒だから．
/*
int getOneFrameSegment(double *f0, int tLen, double **specgram, double **aperiodicity, int fftl, double framePeriod, double currentTime, int fs, double defaultF0,
						fftw_complex *spectrum, fftw_complex *cepstrum, 
						double *response, int xLen);
*/
/*
void getMinimumPhaseSpectrum(double *inputSpec, fftw_complex *spectrum, fftw_complex *cepstrum, int fftl)
{
	int i;
	double real, imag;
	fftw_plan forwardFFT, inverseFFT;
	forwardFFT = fftw_plan_dft_1d(fftl, spectrum, cepstrum, FFTW_FORWARD, FFTW_ESTIMATE);
	inverseFFT = fftw_plan_dft_1d(fftl, cepstrum, spectrum, FFTW_BACKWARD, FFTW_ESTIMATE);
	// 値を取り出す
	for(i = 0;i <= fftl/2;i++)	
	{
		spectrum[i][0] = log(inputSpec[i])/2.0;
		spectrum[i][1] = 0.0;
	}
	for(;i < fftl;i++)
	{
		spectrum[i][0] = spectrum[fftl-i][0];
		spectrum[i][1] = 0.0;
	}
	fftw_execute(forwardFFT);
	for(i = 1;i < fftl/2;i++)
	{
		cepstrum[i][0] = 0.0;
		cepstrum[i][1] = 0.0;
	}
	for(;i < fftl;i++)
	{
		cepstrum[i][0] *= 2.0;
		cepstrum[i][1] *= 2.0;
	}
	fftw_execute(inverseFFT);
	for(i = 0;i < fftl;i++)
	{
		real = exp(spectrum[i][0]/(double)fftl)*cos(spectrum[i][1]/(double)fftl);
		imag = exp(spectrum[i][0]/(double)fftl)*sin(spectrum[i][1]/(double)fftl);
		spectrum[i][0] = real;
		spectrum[i][1] = imag;
	}
}
*/
// 特定時刻の応答を取得する．
/*
void getOneFrameSegment(double *f0, int tLen, double **specgram, double **residualSpecgram, int fftl, double framePeriod, double currentTime, int fs, double defaultF0,
						fftw_complex *spectrum, fftw_complex *cepstrum, 
						double *response, int xLen)
{
	int i;
	double real, imag, tmp;
	fftw_plan	inverseFFT_RP;				// FFTセット
	int currentFrame, currentPosition;
	inverseFFT_RP = fftw_plan_dft_c2r_1d(fftl, spectrum, response ,  FFTW_ESTIMATE);
	currentFrame = (int)(currentTime/(framePeriod/1000.0) + 0.5);	
	currentPosition = (int)(currentTime*(double)fs);
	tmp = currentTime + 1.0/(f0[currentFrame] == 0.0 ? defaultF0 : f0[currentFrame]);
	// 値を取り出す
	getMinimumPhaseSpectrum(specgram[currentFrame], spectrum, cepstrum, fftl);
	spectrum[0][0] *= residualSpecgram[currentFrame][0];
	for(i = 1;i <= fftl/2;i++)
	{
		real = spectrum[i][0]*residualSpecgram[currentFrame][(i-1)*2+1] - spectrum[i][1]*residualSpecgram[currentFrame][i*2];
		imag = spectrum[i][0]*residualSpecgram[currentFrame][i*2] + spectrum[i][1]*residualSpecgram[currentFrame][(i-1)*2+1];
		spectrum[i][0] = real;
		spectrum[i][1] = imag;
	}
	fftw_execute(inverseFFT_RP);
	fftw_destroy_plan(inverseFFT_RP);
}
*/


void synthesisPt100(double *f0, int tLen, double **aperiodicity, int fftl, double framePeriod, int fs, 
			   double *synthesisOut, int xLen)
{
	int i,j;

	double currentTime = 0.0;
	int currentPosition = 0;//currentTime / framePeriod;
	int currentFrame = 0;
	for(i = 0;;i++)
	{
		currentPosition = (int)(currentTime*(double)fs);//周期単位で継ぎ足していく
		for(j = 0;j < fftl/2;j++)
		{
			if(j+currentPosition >= xLen) break;
			synthesisOut[j+currentPosition] += aperiodicity[currentFrame][j];
		}

		// 更新
		currentTime += 1.0/(f0[currentFrame] == 0.0 ? DEFAULT_F0 : f0[currentFrame]);//時刻を1周期分進める
		currentFrame = (int)(currentTime/(framePeriod/1000.0) + 0.5);//次に継ぎ足すデータ位置は次の時刻に最も近いフレーム
		currentPosition = (int)(currentTime*(double)fs);
		if(j+currentPosition >= xLen || currentFrame >= tLen) break;
	}

	return;
}

void synthesisPt101(double fixedDefault_f0, double *f0, int tLen, double **aperiodicity, int *ResidualSpecgramLength,
					int *fixedResidualSpecgramIndex, double *volume,
					int fftl, double framePeriod, int fs, double *synthesisOut, int xLen)
{
	int i,j;
	double currentTime = 0;
	int currentPosition = 0;//currentTime / framePeriod;
	int currentFrame = 0;

	for(i = 0;;i++)
	{
		for(j = 0;j < ResidualSpecgramLength[fixedResidualSpecgramIndex[currentFrame]];j++)
		{
			if(j+currentPosition >= xLen) break;
			synthesisOut[max(0, j+currentPosition)] += aperiodicity[fixedResidualSpecgramIndex[currentFrame]][j] * volume[currentFrame];
		}

		// 更新
		currentTime += 1.0/(f0[currentFrame] == 0.0 ? fixedDefault_f0 : f0[currentFrame]);//時刻を1周期分進める
		currentFrame = (int)(currentTime/(framePeriod/1000.0) + 0.5);//次に継ぎ足すデータ位置は次の時刻に最も近いフレーム
		currentPosition = (int)(currentTime*(double)fs);//周期単位で継ぎ足していく

		if(j+currentPosition >= xLen || currentFrame >= tLen) break;
	}
	return;
}
