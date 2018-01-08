#include "world.h"

#include <stdio.h> // for debug
#include <stdlib.h>
#include <math.h>

#define NOMINALCUTOFF 600.0

// QMFの設計．これが一番の鬼門だから妥協しよう．
void designQMFpairOfFilters(int fs, double *hHP, double *hLP);

void setH(double *x, int xLen, int segmentLength, int nMargin2, int indexBias, int currentPositionInSample, int t0InSamples, 
		  double **H)
{
	int i, j;

	int index;
	for(i = -1;i < 2;i++)
	{
		for(j = 0;j < segmentLength;j++)
		{
			index = max(0, min(xLen-1, i+currentPositionInSample-indexBias-t0InSamples+j));
			H[j][i+1] = x[index];
			index = max(0, min(xLen-1, i+currentPositionInSample-indexBias+t0InSamples+j));
			H[j][i+4] = x[index];
		}
	}
}

void calcHw(double **H, int segmentLength, int nMargin2, double **w,
			double **Hw)
{
	int i,j,k;
	double tmp;
	for(i = 0;i < nMargin2;i++)
	{
		for(j = 0;j < segmentLength;j++)
		{
			tmp = 0.0;
			for(k = 0;k < segmentLength;k++)
			{
				tmp += H[k][i] * w[k][j];
				
			}
			Hw[i][j] = tmp;
		}
	}
}

void calcR(double **Hw, int nMargin2, int segmentLength, double **H,
		   double **R)
{
	int i,j,k;
	double tmp;
	for(i = 0;i < nMargin2;i++)
	{
		for(j = 0;j < nMargin2;j++)
		{
			tmp = 0.0;
			for(k = 0;k < segmentLength;k++)
			{
				tmp += Hw[i][k] * H[k][j];
				
			}
			R[i][j] = tmp;
		}
	}
}

void calcHwx(double **Hw, int nMargin2, int segmentLength, double *x, int origin, 
			 double *Hwx)
{
	int i,j;
	double tmp;
	for(i = 0;i < nMargin2;i++)
	{
		tmp = 0.0;
		for(j = 0;j < segmentLength;j++)
		{
			tmp += Hw[i][j]*x[origin+j];
		}
		Hwx[i] = tmp;
	}
}

void calca(double **invR, int nMargin2, double *Hwx,
		   double *a)
{
	int i,j;
	double tmp;
	for(i = 0;i < nMargin2;i++)
	{
		tmp = 0.0;
		for(j = 0;j < nMargin2;j++)
		{
			tmp += invR[i][j]*Hwx[j];
		}
		a[i] = tmp;
	}
}

void calcHa(double **H, int segmentLength, int nMargin2, double *a,
			double *Ha)
{
	int i,j;
	double tmp;
	for(i = 0;i < segmentLength;i++)
	{
		tmp = 0.0;
		for(j = 0;j < nMargin2;j++)
		{
			tmp += H[i][j]*a[j];
		}
		Ha[i] = tmp;
	}
}

double calcStdwxHa(double *wsqrt, int segmentLength, double *x, int origin, double *Ha, double *wxHa)
{
	int i;
	for(i = 0;i < segmentLength;i++)
	{
		wxHa[i] = wsqrt[i]*(x[i+origin]-Ha[i]);
	}
	return std(wxHa, segmentLength);
}

double calcStdwx(double *wsqrt, int segmentLength, double *x, int origin, double *wx)
{
	int i;
	for(i = 0;i < segmentLength;i++)
	{
		wx[i] = wsqrt[i]*x[i+origin];
	}
	return std(wx, segmentLength);
}

void f0PredictionResidualFixSegmentW(double *x, int xLen, double fs, double targetF0, double *temporalPositions, double *vuv, int tLen, double initialTime, int durationMs,
									 double *rmsResidual)
{
	int i, j;
	int nMargin = 3;
	int segmentLength;

	segmentLength = (int)(fs*durationMs/2000.0+0.5)*2 + 1;

	//---------------------------------------------------------------------------------
	// メモリ確保
	double **w;
	w = (double **)malloc(sizeof(double *) * segmentLength);
	for(i = 0;i < segmentLength;i++) w[i] = (double *)malloc(sizeof(double) * segmentLength);

	for(i = 0;i < segmentLength;i++)
		for(j = 0;j < segmentLength;j++)
			w[i][j] = 0.0;
	for(i = 0;i < (segmentLength-1)/2;i++)
	{
		w[i][i] = 0.5 - 0.5*cos(((double)i+1.0)/(double)(segmentLength+1)*2.0*PI);
		w[segmentLength-i-1][segmentLength-i-1] = w[i][i];
	}
	w[(segmentLength-1)/2][(segmentLength-1)/2] = 1.0;

	double *wsqrt;
	wsqrt = (double *)malloc(sizeof(double)*segmentLength);
	for(i = 0;i < segmentLength;i++)
		wsqrt[i] = sqrt(w[i][i]);

	double **H, **Hw, **R, **invR;
	H    = (double **)malloc(sizeof(double *)*segmentLength);
	for(i = 0;i < segmentLength;i++)
		H[i] = (double *)malloc(sizeof(double)*(nMargin*2) );
	Hw   = (double **)malloc(sizeof(double *)*nMargin*2);
	for(i = 0;i < nMargin*2;i++)
		Hw[i] = (double *)malloc(sizeof(double)*segmentLength );
	R    = (double **)malloc(sizeof(double *)*nMargin*2);
	for(i = 0;i < nMargin*2;i++)
		R[i] = (double *)malloc(sizeof(double)*nMargin*2 );
	invR = (double **)malloc(sizeof(double *)*nMargin*2);
	for(i = 0;i < nMargin*2;i++)
		invR[i] = (double *)malloc(sizeof(double)*nMargin*2 );

	double *Hwx, *a, *Ha, *wxHa, *wx;
	Hwx   = (double *)malloc(sizeof(double) * nMargin*2);
	a     = (double *)malloc(sizeof(double) * nMargin*2);
	Ha	  = (double *)malloc(sizeof(double) * segmentLength);
	wxHa  = (double *)malloc(sizeof(double) * segmentLength);
	wx    = (double *)malloc(sizeof(double) * segmentLength);

	double *segmentIndex;
	segmentIndex = (double *)malloc(sizeof(double) * segmentLength);
	//---------------------------------------------------------------------------------

	int t0InSamples;
	int currentPositionInSample;
	int indexBias;
	int origin;

	t0InSamples = (int)(fs/targetF0 + 0.5);
	indexBias = (int)(fs/targetF0/2.0 + 0.5);

	//ok
	for(i = 0;i < tLen;i++)
	{
		currentPositionInSample = (int)(-initialTime+temporalPositions[i]*fs + 0.5)+1;

		if(vuv[i] != 0.0)
		{
			origin = max(0, min(xLen-1, currentPositionInSample-indexBias));
			setH(x, xLen, segmentLength, nMargin*2, indexBias, currentPositionInSample, t0InSamples, H);
			calcHw(H, segmentLength, nMargin*2, w, Hw);
			calcR(Hw, nMargin*2, segmentLength, H, R);
			calcHwx(Hw, nMargin*2, segmentLength, x, origin, Hwx);
			inv(R, nMargin*2, invR);
			calca(invR, nMargin*2, Hwx, a);
			calcHa(H, segmentLength, nMargin*2, a, Ha);
			rmsResidual[i] = calcStdwxHa(wsqrt, segmentLength, x, origin, Ha, wxHa) / 
				calcStdwx(wsqrt, segmentLength, x, origin, wx);
		}
		else
		{ // 無声音の場合，非周期性指標は意味を持たない．
			rmsResidual[i] = 0.0;
		}
	}

	// 使ったメモリの開放
	free(segmentIndex);
	free(wx); free(wxHa); free(Ha); free(a); free(Hwx);
	for(i = 0;i < nMargin*2;i++)	 free(invR[i]);  
	free(invR);
	for(i = 0;i < nMargin*2;i++)	 free(R[i]);  
	free(R);
	for(i = 0;i < nMargin*2;i++)	 free(Hw[i]); 
	free(Hw);
	for(i = 0;i < segmentLength;i++) free(H[i]);  
	free(H);
	free(wsqrt);
	for(i = 0;i < segmentLength;i++) free(w[i]); 
	free(w);
}

// 帯域毎の非周期性指標を計算する．
// この関数は，QMFによるサブバンド分解を主に担当．
// 戻り値未定なので，とりあえずの実装
void bandwiseAperiodicity(double *x, int xLen, int fs, double *f0, double *vuv, double *stretchedLocations, int tLen, int windowLengthMs,
						  double **aperiodicity)
{
	int i, j;
	double hHP[41], hLP[37];
	double nominalCutOff = NOMINALCUTOFF;

	designQMFpairOfFilters(fs, hHP, hLP);

	int cLen;
	cLen = (int)(log((double)fs/nominalCutOff)/log(2.0));

	double *cutOffList;
	cutOffList = (double *)malloc(sizeof(double) * cLen);

	for(i = 0;i < cLen;i++)
		cutOffList[i] = (double)fs / pow(2.0, (double)(i+2)); 

	int fftl;
	fftl = (int)pow(2.0, 1.0 + (int)(log((double)(xLen+82) ) / log(2.0))); // hHPの倍の長さ分余計に確保する必要がある

	double *wholeSignal, *highSignal, *lowSignal, *downSampledHighSignal;

	wholeSignal				= (double *)malloc(sizeof(double) * fftl);
	highSignal				= (double *)malloc(sizeof(double) * fftl);
	lowSignal				= (double *)malloc(sizeof(double) * fftl);
	downSampledHighSignal   = (double *)malloc(sizeof(double) * fftl);

	int wLen = xLen+82;
	double fsTmp;
	double *rmsResidual;
	rmsResidual = (double *)malloc(sizeof(double) * tLen);

	for(i = 0;i < xLen;i++) wholeSignal[i] = x[i];
	for(;i < fftl;i++) wholeSignal[i] = 0.0;
	for(i = 0;i < cLen-1;i++)
	{
		fsTmp = cutOffList[i]*2.0;
		fftfilt(wholeSignal, wLen, hHP, 41, fftl, highSignal);
		fftfilt(wholeSignal, wLen, hLP, 37, fftl, lowSignal );
		for(j = 0;j < wLen;j+=2) downSampledHighSignal[j/2] = highSignal[j];

		f0PredictionResidualFixSegmentW(downSampledHighSignal, (int)((double)wLen/2.0+0.5), fsTmp, f0[0], stretchedLocations, vuv, tLen, 41.0/2.0/fsTmp, windowLengthMs, rmsResidual);

		// ある帯域の非周期性指標をコピー
		for(j = 0;j < tLen;j++) aperiodicity[j][cLen-i-1] = rmsResidual[j];

		// サブバンド分解 (?)
		for(j = 0;j < wLen; j+=2) wholeSignal[j/2] = lowSignal[j];
		wLen = (int)((double)wLen/2.0 + 0.5) + 82;
		fftl = (int)pow(2.0, 1.0 + (int)(log((double)(wLen) ) / log(2.0)));
		for(j=j/2;j<fftl;j++) wholeSignal[j] = 0.0;

	}

	wLen = (wLen-82)*2;
	f0PredictionResidualFixSegmentW(wholeSignal, (int)((double)wLen/2.0+0.5), fsTmp, f0[0], stretchedLocations, vuv, tLen, 41.0/2.0/fsTmp, windowLengthMs, rmsResidual);
	for(j = 0;j < tLen;j++) aperiodicity[j][0] = rmsResidual[j];

	free(rmsResidual);
	free(downSampledHighSignal); free(lowSignal); free(highSignal); free(wholeSignal);
	free(cutOffList);
}

// 4倍にアップサンプリングするための関数 (1)
void setInterpolatedX(double *x, int xLen, double *interpolatedX)
{
	int i;

	interpolatedX[0] = x[0] * 0.14644660940672621;
	interpolatedX[1] = x[0] * 0.49999999999999994;
	interpolatedX[2] = x[0] * 0.85355339059327373;
	for(i = 0;i < xLen-1;i++)
	{
		interpolatedX[i*4 + 3] = x[i];
		interpolatedX[i*4 + 4] = x[i]*0.85355339059327373 + x[i+1]*0.14644660940672621;
		interpolatedX[i*4 + 5] = x[i]*0.49999999999999994 + x[i+1]*0.49999999999999994;
		interpolatedX[i*4 + 6] = x[i]*0.14644660940672621 + x[i+1]*0.85355339059327373;
	}
	interpolatedX[i*4 + 3] = x[xLen-1];
	interpolatedX[i*4 + 4] = x[xLen-1] * 0.85355339059327373;
	interpolatedX[i*4 + 5] = x[xLen-1] * 0.49999999999999994;
	interpolatedX[i*4 + 6] = x[xLen-1] * 0.14644660940672621;
	interpolatedX[i*4 + 7] = interpolatedX[i*4 + 8] = interpolatedX[i*4 + 9] = 0.0;
	return;
}

// 信号の見かけ上のF0を一定 (targetF0 Hz)にする．
void normalizeSignal(double *x, int xLen, int fs, double *f0, int tLen, double framePeriod, double targetF0, 
					 double **stretchedSignal, double **stretchedLocations, int *oLen)
{
	int i;
	int ixLen = xLen*4 + 6;

	double *interpolatedX;
	interpolatedX = (double *)malloc(sizeof(double) * (ixLen+16) );
	
	setInterpolatedX(x, xLen, interpolatedX);

	double *originalSignalTime;
	originalSignalTime = (double *)malloc(sizeof(double) * ixLen );

	for(i = 0;i < ixLen;i++)
		originalSignalTime[i] = (double)i/( (double)fs * 4.0);

	double *baseF0, *baseTime;
	baseF0   = (double *)malloc(sizeof(double) * (tLen+1));
	baseTime = (double *)malloc(sizeof(double) * (tLen+1));

	for(i = 0;i < tLen;i++)
	{
		baseF0[i]   = f0[i] == 0.0 ? targetF0 : f0[i];
		baseTime[i] = (double)i * framePeriod;
	}
	baseF0[tLen]   = targetF0;
	baseTime[tLen] = (double)i * framePeriod;

	double *interpolatedF0, *stretchedTime;
	interpolatedF0 = (double *)malloc(sizeof(double) * ixLen);
	stretchedTime  = (double *)malloc(sizeof(double) * ixLen);

	interp1(baseTime, baseF0, tLen+1, originalSignalTime, ixLen, interpolatedF0);

	double tmp;
	tmp = targetF0*( (double)fs*4.0);
	stretchedTime[0] = interpolatedF0[0]/tmp;
	for(i = 1;i < ixLen;i++)
		stretchedTime[i] = stretchedTime[i-1] + (interpolatedF0[i]/tmp);

	double *tmpTime, *stretchedSignal4;
	int sLen;
	sLen = (int)(stretchedTime[ixLen-1] * (double)fs * 4.0) + 1;
	tmpTime          = (double *)malloc(sizeof(double) * sLen);
	stretchedSignal4 = (double *)malloc(sizeof(double) * sLen);

	for(i = 0;i < sLen;i++)
		tmpTime[i] = (double)i / ( (double)fs*4.0);
	interp1(stretchedTime, interpolatedX, ixLen, tmpTime, sLen, stretchedSignal4);

	*stretchedLocations = (double *)malloc(sizeof(double) * tLen);
	interp1(originalSignalTime, stretchedTime, ixLen, baseTime, tLen, *stretchedLocations);

	*oLen = 1 + (int)((double)sLen/4.0);
	*stretchedSignal = (double *)malloc(sizeof(double) * (*oLen+16) );
	double *tmpSignal = *stretchedSignal;

	decimateForF0(stretchedSignal4, sLen, tmpSignal, 4);

	// メモリの開放
	free(stretchedSignal4);
	free(tmpTime);
	free(stretchedTime);
	free(baseF0);
	free(baseTime);
	free(interpolatedF0);
	free(originalSignalTime);
	free(interpolatedX);
}

int getBands_v3(int fs)
{
	double nominalCutOff = NOMINALCUTOFF;
	return (int)(log((double)fs/nominalCutOff)/log(2.0));
}

// TANDEM最新版で採用されている非周期性指標の移植
// 機能は一部限定的だが，品質はほぼ等価 (だと嬉しい)．
void aperiodicityRatio_v3(double *x, int xLen, int fs, double *f0, int tLen, double framePeriod, 
		 double **aperiodicity, double *targetF0_output)
{
	int i;
	int sLen;
	double maxF0, targetF0;
	maxF0 = 0.0;
	for(i = 0;i < tLen;i++)
		maxF0 = maxF0 > f0[i] ? maxF0 : f0[i];
	targetF0 = max(32.0, min(200.0, maxF0));
	*targetF0_output = targetF0;

	double *stretchedSignal, *stretchedLocations;
	stretchedSignal = stretchedLocations = NULL;

	normalizeSignal(x, xLen, fs, f0, tLen, framePeriod/1000.0, targetF0,
		&stretchedSignal, &stretchedLocations, &sLen);

	double *stretchedF0;
	stretchedF0 = (double *)malloc(sizeof(double) * tLen);

	for(i = 0;i < tLen;i++)
		stretchedF0[i] = targetF0;
	bandwiseAperiodicity(stretchedSignal, sLen, fs, stretchedF0, f0, stretchedLocations, tLen, (int)(2000.0/targetF0 + 0.5),
		aperiodicity);

	free(stretchedF0);
	free(stretchedSignal);
	free(stretchedLocations);
}

void calculateAperiodicity(double *aperiodicity, int cLen, int fftl, double f0, int fs, double targetF0, double *periodicSpec)
{
	int i;
	double *ap, *axis, *w, *tmpAP;
	ap		= (double *)malloc(sizeof(double) * (cLen+1) );
	axis	= (double *)malloc(sizeof(double) * (cLen+1) );
	w	    = (double *)malloc(sizeof(double) * (fftl/2+1) );
	tmpAP	= (double *)malloc(sizeof(double) * (fftl/2+1) );

	double *cutOffList;
	cutOffList = (double *)malloc(sizeof(double) * cLen);
	for(i = 0;i < cLen;i++)
		cutOffList[i] = (double)fs / pow(2.0, (double)(i+2));

	double f0Safe, stretchingFactor;
	f0Safe = max(f0, targetF0);
	stretchingFactor = f0Safe/targetF0;

	ap[0]   = log(0.0000000005); // セーフガード，実際はまぁ，色々とある
	axis[0] = 0.0;
	for(i = 0;i < cLen-1;i++)
	{
		ap[i+1] = log(aperiodicity[i]);
		axis[i+1] = cutOffList[cLen-i-2];
	}
	ap[cLen]   = log(aperiodicity[cLen-1]);
	axis[cLen] = (double)fs/2.0;

	for(i = 0;i <= fftl/2;i++)
	{
		w[i] = (double)(i * fs)/ (double)fftl;
	}
	interp1(axis, ap, cLen+1, w, fftl/2+1, tmpAP);
	for(i = 0;i < cLen-1;i++)
	{
		axis[i+1] *= stretchingFactor;
	}
	axis[cLen] = (double)fs/2.0 * stretchingFactor;

	interp1(axis, ap, cLen+1, w, fftl/2+1, periodicSpec);

	for(i = 0;i <= fftl/2;i++)
	{
		periodicSpec[i] = 1.0 - min(exp(tmpAP[i]*2.0), exp(periodicSpec[i]*2.0));
	}

	free(tmpAP);
	free(cutOffList);
	free(w);
	free(axis);
	free(ap);
}

//-----------------------------------------------------
// どうでも良い関数ほど下にある
// hHP:41, hLP:37の要素を持つdouble型であること．
// 一番面倒なQMFの設計．結局妥協した．
// Matlab版もfsを引数にするが，結果はfsに依存せず一定．これって何だろう？
void designQMFpairOfFilters(int fs, double *hHP, double *hLP)
{
	// hHP
	hHP[0]  =  0.00041447996898231424;
	hHP[1]  =  0.00078125051417292477;
	hHP[2]  = -0.0010917236836275842;
	hHP[3]  = -0.0019867925675967589;
	hHP[4]  =  0.0020903896961562292;
	hHP[5]  =  0.0040940570272849346;
	hHP[6]  = -0.0034025808529816698;
	hHP[7]  = -0.0074961541272056016;
	hHP[8]  =  0.0049722633399330637;
	hHP[9]  =  0.012738791249119802;
	hHP[10] = -0.0066960326895749113;
	hHP[11] = -0.020694051570247052;
	hHP[12] =  0.0084324365650413451;
	hHP[13] =  0.033074383758700532;
	hHP[14] = -0.010018936738799522;
	hHP[15] = -0.054231361405808247;
	hHP[16] =  0.011293988915051487;
	hHP[17] =  0.10020081367388213;
	hHP[18] = -0.012120546202484579;
	hHP[19] = -0.31630021039095702;
	hHP[20] =  0.51240682580627639;
	hHP[21] = -0.31630021039095702;
	hHP[22] = -0.012120546202484579;
	hHP[23] =  0.10020081367388213;
	hHP[24] =  0.011293988915051487;
	hHP[25] = -0.054231361405808247;
	hHP[26] = -0.010018936738799522;
	hHP[27] =  0.033074383758700532;
	hHP[28] =  0.0084324365650413451;
	hHP[29] = -0.020694051570247052;
	hHP[30] = -0.0066960326895749113;
	hHP[31] =  0.012738791249119802;
	hHP[32] =  0.0049722633399330637;
	hHP[33] = -0.0074961541272056016;
	hHP[34] = -0.0034025808529816698;
	hHP[35] =  0.0040940570272849346;
	hHP[36] =  0.0020903896961562292;
	hHP[37] = -0.0019867925675967589;
	hHP[38] = -0.0010917236836275842;
	hHP[39] =  0.00078125051417292477;
	hHP[40] =  0.00041447996898231424;

	// hLP
	hLP[0]  = -0.00065488170077483048;
	hLP[1]  =  0.00007561994958159384;
	hLP[2]  =  0.0020408456937895227;
	hLP[3]  = -0.00074680535322030437;
	hLP[4]  = -0.0043502235688264931;
	hLP[5]  =  0.0025966428382642732;
	hLP[6]  =  0.0076396022827566962;
	hLP[7]  = -0.0064904118901497852;
	hLP[8]  = -0.011765804538954506;
	hLP[9]  =  0.013649908479276255;
	hLP[10] =  0.01636866479016021;
	hLP[11] = -0.026075976030529347;
	hLP[12] = -0.020910294856659444;
	hLP[13] =  0.048260725032316647;
	hLP[14] =  0.024767846611048111;
	hLP[15] = -0.096178467583360641;
	hLP[16] = -0.027359756709866623;
	hLP[17] =  0.31488052161630042;
	hLP[18] =  0.52827343594055032;
	hLP[19] =  0.31488052161630042;
	hLP[20] = -0.027359756709866623;
	hLP[21] = -0.096178467583360641;
	hLP[22] =  0.024767846611048111;
	hLP[23] =  0.048260725032316647;
	hLP[24] = -0.020910294856659444;
	hLP[25] = -0.026075976030529347;
	hLP[26] =  0.01636866479016021;
	hLP[27] =  0.013649908479276255;
	hLP[28] = -0.011765804538954506;
	hLP[29] = -0.0064904118901497852;
	hLP[30] =  0.0076396022827566962;
	hLP[31] =  0.0025966428382642732;
	hLP[32] = -0.0043502235688264931;
	hLP[33] = -0.00074680535322030437;
	hLP[34] =  0.0020408456937895227;
	hLP[35] =  0.00007561994958159384;
	hLP[36] = -0.00065488170077483048;
}
