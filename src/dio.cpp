// MacRes v0.0.1   3/10/2017
#include <math.h>

#include "world.h"

// Internal functions, please don't use from outside DIO.
void rawEventByDio(double boundary_freq, double sample_rate, fft_complex *xSpec, int xLength, int fft_length, double shiftTime, double f0_floor, double f0_ceil, double *time_axis, int num_frames, 
				   double *f0Deviations, double *interpolatedF0,
				   fft_plan &forwardFFT, fft_plan &inverseFFT, double *equivalentFIR, fft_complex *eSpec);
void zeroCrossingEngine(double *waveform, int num_samples, double sample_rate,
						double *eLocations, double *iLocations, double *intervals, int *iLen);
long decimateForF0(double *waveform, int num_samples, double *y, int r);
void filterForDecimate(double *waveform, int num_samples, double *y, int r);
void nuttallWindow(int target_num_samples, double *y);
void postprocessing(double frame_period, double f0_floor, int candidates, int num_samples, int sample_rate, double **f0Map,
	               double **stability_map, double *bestF0, double *f0);
void interp1(double *t, double *y, int iLen, double *t1, int oLen, double *y1);
void histc(double *waveform, int num_samples, double *y, int target_num_samples, int *index);

// Get the number of elements in the F0 track（So we can reserve memory for them in advance）
// frame_period's units are milliseconds
int GetNumDIOSamples(int sample_rate, int num_samples, double frame_period)
{
	// Total milliseconds divided by milliseconds per frame, plus one.
	return (int)((double)num_samples / (double)sample_rate * 1000.0 / frame_period) + 1;
}

// DIO (Distributed Inline filter Operation) to filter waveform into F0.
// waveform	: Input signal
// num_samples : Length of sample
// f0	: where to put result
// time_axis: Already populated, the time axis for f0.
void dio(double *waveform, int num_samples, int sample_rate, double frame_period, 
		 double *time_axis, double *f0)
{
	int i, j; // For incrementing loops.

	// Initial conditions (can adjust these to try to make improvements)
	double f0_floor = 80;
	double f0_ceil = 640;
	double bands_per_octave = 2;
	double target_sample_rate = 4000;

	// Calculate foundational parameters.
	int decimation_ratio = (int)(sample_rate/target_sample_rate);
	double event_sample_rate = (double)sample_rate/(double)decimation_ratio;
	int num_bands = (int)(log((double)f0_ceil/(double)f0_floor)/log(2.0) * bands_per_octave);

	// List of boundary frequencies, including the floor and ceiling.
	double *boundary_freqs = (double *)malloc(sizeof(double) * (num_bands+1));
	for(i = 0;i <= num_bands;i++) {
		boundary_freqs[i] = f0_floor * pow(2.0, i / bands_per_octave);
	}

	// Calculate the fft length
	int target_num_samples = (1 + (int)(num_samples/decimation_ratio));
	int fft_length = (int)pow(2.0, 1.0 + (int)(log((double)target_num_samples + 
		(double)(4*(int)(1.0 + (double)sample_rate/boundary_freqs[0]/2.0)) ) / log(2.0)));
	double *y = (double *)malloc(sizeof(double) * fft_length);
	
	// Downsampling.
	decimateForF0(waveform, num_samples, y, decimation_ratio);

	// Get rid of the direct current component. y = y - mean(y)
	double meanY = 0.0;
	for(i = 0;i < target_num_samples;i++) {
		meanY += y[i];
	}
	meanY /= (double)target_num_samples;
	for(i = 0;i < target_num_samples;i++) {
		y[i] -= meanY;
	}

	// Zero out the remainder of y.
	for(i = target_num_samples; i < fft_length;i++) {
		y[i] = 0.0;
	}

	// Preserve the interim data.
	int		num_frames; // Number of samples in the F0 track.
	num_frames = GetNumDIOSamples(sample_rate, num_samples, frame_period);
	int lengthInMs = 1 + (int)((double)num_samples/(double)sample_rate*1000.0);
	// Gathers together all candidates in f0Map.
	double **stability_map, ** f0Map;
	stability_map = (double **)malloc(sizeof(double *) * (num_bands+1));
	f0Map		 = (double **)malloc(sizeof(double *) * (num_bands+1));
	for(i = 0;i <= num_bands;i++)
	{
		stability_map[i] = (double *)malloc(sizeof(double) * num_frames);
		f0Map[i]		= (double *)malloc(sizeof(double) * num_frames);
	}

	// Calculate the wave spectrogram in advance（Some room for optimization here）
	fft_plan			forwardFFT;				// FFT set
	fft_complex		*ySpec;	// Spectrogram
	ySpec = (fft_complex *)malloc(sizeof(fft_complex) * fft_length);
	forwardFFT = fft_plan_dft_r2c_1d(fft_length, y, ySpec, FFT_ESTIMATE);
	fft_execute(forwardFFT); // Execute FFT.

	// temporary values
	double *	interpolatedF0;
	double *	f0Deviations;
	interpolatedF0 = (double *) malloc(sizeof(double) * lengthInMs);
	f0Deviations   = (double *) malloc(sizeof(double) * lengthInMs);

	// So FFTW can reuse the plan, the memory is framed/preserved from the DIO side. (?)
	fft_destroy_plan(forwardFFT);
	double *equivalentFIR;
	equivalentFIR = (double *)malloc(sizeof(double) * fft_length);
	//	fft_plan	forwardFFT;				// FFT Set
	fft_plan	inverseFFT;
	fft_complex		*eSpec;	// Spectrogram
	eSpec = (fft_complex *)malloc(sizeof(fft_complex) * fft_length);
	forwardFFT = fft_plan_dft_r2c_1d(fft_length, equivalentFIR, eSpec, FFT_ESTIMATE);
	inverseFFT = fft_plan_dft_c2r_1d(fft_length, eSpec, equivalentFIR, FFT_ESTIMATE);

	// Event calculation (The four zero crossings are listed in detail in the thesis.)
	for(i = 0;i <= num_bands;i++)
	{
		rawEventByDio(boundary_freqs[i], event_sample_rate, ySpec, target_num_samples, fft_length, frame_period/1000.0, f0_floor, f0_ceil, time_axis, num_frames, 
			f0Deviations, interpolatedF0, forwardFFT, inverseFFT, equivalentFIR, eSpec);
		for(j = 0; j < num_frames;j++)
		{
			stability_map[i][j] = f0Deviations[j] / (interpolatedF0[j]+0.00000001);
			f0Map[i][j]        = interpolatedF0[j];
		}
	}

	free(equivalentFIR);
	free(eSpec);

	// Select the best candidate (the one most similar to the standard wave)
	double *bestF0;
	bestF0 = (double *)malloc(sizeof(double) * num_frames);
	/* Changed so that postprocessing handles the process of choosing the best candidate.
	double tmp;
	for(i = 0;i < num_frames;i++)
	{
		tmp = stability_map[0][i];
		bestF0[i] = (stability_map[0][i] < 0.002) ? f0Map[0][i] : 0.0;
		for(j = 1;j <= num_bands;j++)
		{
			if(tmp > stability_map[j][i] && stability_map[j][i] < 0.002)
			{
				tmp = stability_map[j][i];
				bestF0[i] = f0Map[j][i];
			}
		}
	}
	*/

	// Uses stability_map to correct F0. Performs other postprocessing.
	postprocessing(frame_period/1000.0, f0_floor, num_bands+1, num_samples, sample_rate, f0Map, stability_map, bestF0, f0);

	// Cleanup (freeing up memory)
	free(bestF0);
	free(interpolatedF0);
	free(f0Deviations);
	fft_destroy_plan(forwardFFT);
	fft_destroy_plan(inverseFFT);
	free(ySpec);
	for(i = 0;i <= num_bands;i++)
	{
		free(stability_map[i]);
		free(f0Map[i]);
	}
	free(stability_map);
	free(f0Map);
	free(boundary_freqs);
	free(y);
}

// Check whether this number of events happened.
// Last resort in case number of events exceeds long range. (?)
int checkEvent(int num_events)
{
	if(num_events > 0) return 1;
	return 0;
}

// Cleanup (4th step)
void postprocessing(double frame_period, double f0_floor, int candidates, int num_samples, int sample_rate, double **f0Map,
						double **stability_map, double *bestF0, double *f0)
{
	int i, j, k;
	int voiceRangeMinimum = (int)(0.5 + 1.0/frame_period/f0_floor);
	int f0Len = (int)((double)num_samples / (double)sample_rate / frame_period) + 1;
	// This is a 5ms standard, so adjust according to the frame_period.
	double allowedRange = 0.1 * frame_period/0.005;

	// Select best F0 (Moved here from the main DIO function)
	double tmp;
	double *bestF0Stab;
	bestF0Stab = (double *)malloc(sizeof(double) * f0Len);
	for(i = 0;i < f0Len;i++)
	{
		tmp = stability_map[0][i];
		bestF0[i] = 0.0;
		bestF0Stab[i] = stability_map[0][i];
		for(j = 1;j < candidates;j++)
		{
			if(tmp > stability_map[j][i])
			{
				tmp = stability_map[j][i];
				bestF0[i] = f0Map[j][i];
				bestF0Stab[i] = stability_map[j][i];
			}
		}
	}

	// Exclude F0's that are not very stable.
	int addCount = 0;
	double addValue = 0.0;
	for(i = 0;i < f0Len;i++) 
	{
		if(bestF0Stab[i] < 0.05){ addCount++; addValue += bestF0Stab[i];} 
	}
	addValue = addValue * 2.0 / addCount;
	for(i = 0;i < f0Len;i++) if(bestF0Stab[i] > addValue) bestF0[i] = 0.0;

	// It's possible to conserve more memory here, but the data size isn't too large so
	// multiple arrays are used to make code more readable.
	double *f0Base;
	f0Base = (double *)malloc(sizeof(double) * f0Len);
	double *f0Step1;
	f0Step1 = (double *)malloc(sizeof(double) * f0Len);
	double *f0Step2;
	f0Step2 = (double *)malloc(sizeof(double) * f0Len);
	double *f0Step3;
	f0Step3 = (double *)malloc(sizeof(double) * f0Len);
	double *f0Step4;
	f0Step4 = (double *)malloc(sizeof(double) * f0Len);

	// Initialization
	for(i = 0;i < voiceRangeMinimum;i++) f0Base[i] = 0;
	for(;i < f0Len-voiceRangeMinimum;i++) f0Base[i] = bestF0[i];
	for(;i < f0Len;i++) f0Base[i] = 0;
	for(i = 0;i < f0Len;i++) f0Step1[i] = 0.0;

	// Step 1 (Prevent F0 from skipping)
	for(i = voiceRangeMinimum;i < f0Len;i++)
		if(fabs((f0Base[i]-f0Base[i-1])/(0.00001+f0Base[i]) ) < allowedRange)
			f0Step1[i] = f0Base[i];

	// Step 2 (Detach unvoiced sections)
	for(i = 0;i < f0Len;i++) f0Step2[i] = f0Step1[i];
	for(i = voiceRangeMinimum;i < f0Len;i++)
	{
		for(j = 0;j < voiceRangeMinimum;j++)
		{
			if(f0Step1[i-j] == 0)
			{
				f0Step2[i] = 0.0;
				break;
			}
		}
	}

	// Also Step 2 (Detach unvoiced sections) the opposite direction
	for(i = f0Len-1-voiceRangeMinimum;i >= 0;i--)
	{
		for(j = 0;j < voiceRangeMinimum;j++)
		{
			if(f0Step1[i+j] == 0)
			{
				f0Step2[i] = 0.0;
				break;
			}
		}
	}

	// Detect the number of islands
	int *positiveIndex, *negativeIndex;
	positiveIndex = (int *)malloc(sizeof(int) * f0Len);
	negativeIndex = (int *)malloc(sizeof(int) * f0Len);
	int positiveCount, negativeCount;
	positiveCount = negativeCount = 0;
	for(i = 1;i < f0Len;i++)
	{
		if(f0Step2[i] == 0 && f0Step2[i-1] != 0)
			negativeIndex[negativeCount++] = i-1;
		else if (f0Step2[i-1] == 0 && f0Step2[i] != 0)
			positiveIndex[positiveCount++] = i;
	}

	// Step 3（forward-facing corrections）
	double refValue1, refValue2, bestError, errorValue;
	for(i = 0;i < f0Len;i++) f0Step3[i] = f0Step2[i];
	for(i = 0;i < negativeCount;i++)
	{
		for(j = negativeIndex[i];j < f0Len-1;j++)
		{
			if(f0Step3[j+1] != 0) break;
			refValue1 = f0Step3[j]*2 - f0Step3[j-1];
			refValue2 = f0Step3[j];
			//bestError = fabs(refValue - f0Map[0][j+1]);
			bestError = min(fabs(refValue1 - f0Map[0][j+1]), fabs(refValue2 - f0Map[0][j+1]));
			for(k = 1;k < candidates;k++)
			{
				//errorValue = fabs(refValue - f0Map[k][j+1]);
				errorValue = min(fabs(refValue1 - f0Map[k][j+1]), fabs(refValue2 - f0Map[k][j+1]));
				// Don't use low-stability F0's
				if(errorValue < bestError && stability_map[k][j+1] < 0.1)
				{
					bestError = errorValue;
					f0Step3[j+1] = f0Map[k][j+1];
				}
			}
			//if(bestError / (refValue+0.0001) > allowedRange)
			if(min(bestError / (refValue1+0.0001), bestError / (refValue2+0.0001)) > allowedRange)
			{
				f0Step3[j+1] = 0.0;
				break;
			}
			if(i != negativeCount && j == positiveIndex[i+1]-1)
			{
				negativeIndex[j] = j;
				break;
			}
		}
	}

	// Step 4（backward-facing corrections）
	for(i = 0;i < f0Len;i++) f0Step4[i] = f0Step3[i];
	for(i = positiveCount-1;i >= 0;i--)
	{
		for(j = positiveIndex[i]/*+1*/;j > 1;j--) //tn_fnds v0.0.4 
		{
			if(f0Step4[j-1] != 0) break;
			refValue1 = f0Step4[j]*2 - f0Step4[j-1];
			refValue2 = f0Step4[j];
			//refValue = f0Step4[j]*2 - f0Step4[j+1];
			bestError = min(fabs(refValue1 - f0Map[0][j+1]), fabs(refValue2 - f0Map[0][j+1]));
			//bestError = fabs(refValue - f0Map[0][j-1]);
			for(k = 1;k < candidates;k++)
			{
				errorValue = min(fabs(refValue1 - f0Map[k][j-1]), fabs(refValue2 - f0Map[k][j-1]));
				//errorValue = fabs(refValue - f0Map[k][j-1]);
				// Don't use low-stability F0's.
				if(min(bestError / (refValue1+0.0001), bestError / (refValue2+0.0001)) > allowedRange && stability_map[k][j-1] < 0.1)
				//if(errorValue < bestError)
				{
					bestError = errorValue;
					f0Step4[j-1] = f0Map[k][j-1];
				}
			}
			if(min(bestError / (refValue1+0.0001), bestError / (refValue2+0.0001)) > allowedRange)
			//if(bestError / (refValue+0.0001) > allowedRange)
			{
				f0Step4[j-1] = 0.0;
				break;
			}
			if(i != 0 && j == negativeIndex[i-1]+1) break;
		}
	}

	// Copy
	for(i = 0;i < f0Len;i++) f0[i] = f0Step4[i];
	/* Step 5 doesn't seem to improve things, so temporarily leaving it out.
	// Step 5（Detach isolated islands a second time）
	int voiceRangeMinimum2 = 2+(int)(voiceRangeMinimum/2);
	for(i = 0;i < f0Len;i++) f0[i] = f0Step4[i];
	for(i = voiceRangeMinimum2; i < f0Len-voiceRangeMinimum2;i++)
	{
		for(j = 0;j < voiceRangeMinimum2;j++)
		{
			if(f0Step4[i-j] == 0)
				break;
		}
		for(k = 0;k < voiceRangeMinimum2;k++)
		{
			if(f0Step4[i+k] == 0)
				break;
		}
		f0[i] = j != voiceRangeMinimum2 && k != voiceRangeMinimum2 ? 
			0 : f0Step4[i];
	}
	*/

	// For debugging, write the results of candidate correction into a file.
	/*
	FILE *file;
	file = fopen("/tmp/f0map.txt", "w");
	for(k = 1;k < candidates;k++) fprintf(file,",map%d",k);
	for(k = 1;k < candidates;k++) fprintf(file,",stab%d",k);
	fprintf(file,",BEST,base,step1,step2,step3,step4\n");
	for(i = 0; i < f0Len; i++)
	{
		fprintf(file,"%d",i);
		for(k = 1;k < candidates;k++) fprintf(file,",%f",f0Map[k][i]);
		for(k = 1;k < candidates;k++) fprintf(file,",%f",stability_map[k][i]);
		fprintf(file,",%f", bestF0[i]);
		fprintf(file,",%f", f0Base[i]);
		fprintf(file,",%f", f0Step1[i]);
		fprintf(file,",%f", f0Step2[i]);
		fprintf(file,",%f", f0Step3[i]);
		fprintf(file,",%f\n",f0Step4[i]);
	}
	fclose(file);
	*/
	
	// Free all remaining memory.
	free(bestF0Stab);
	free(f0Base);
	free(f0Step1); free(f0Step2); free(f0Step3); free(f0Step4);
}

// Internal function for calculating events (Mutates local variables, doesn't return anything.)
void rawEventByDio(double boundary_freq, double sample_rate, fft_complex *xSpec, int xLength, int fft_length, double frame_period, double f0_floor, double f0_ceil, double *time_axis, int num_frames, 
				   double *f0Deviations, double *interpolatedF0,
				   fft_plan &forwardFFT, fft_plan &inverseFFT, double *equivalentFIR, fft_complex *eSpec)
{
	int i;
	int halfAverageLength = (int)(sample_rate / boundary_freq / 2 + 0.5);
	int indexBias = halfAverageLength*2;
	for(i = halfAverageLength*2;i < fft_length;i++) equivalentFIR[i] = 0.0;
	nuttallWindow(halfAverageLength*4, equivalentFIR);

	// Calculating forwardFFT was moved to the DIO function's main body.
	fft_execute(forwardFFT); // Apply the FFT

	// Complex number multiplication.
	double tmp;
	//for(i = 0;i <= fft_length-1;i++)　
	for(i = 0;i <= fft_length>>1;i++)	// Should be fine with half of FFT lengthe.
	{
		tmp = xSpec[i][0]*eSpec[i][0] - xSpec[i][1]*eSpec[i][1];
		eSpec[i][1] = xSpec[i][0]*eSpec[i][1] + xSpec[i][1]*eSpec[i][0];
		eSpec[i][0] = tmp;
	}

	// Low-pass filter ring
	// Calculating inverseFFT was moved to the DIO function's main body.
	fft_execute(inverseFFT);
	// Remove bias (latency from the low-pass filter.
	for(i = 0;i < xLength;i++) equivalentFIR[i] = equivalentFIR[i+indexBias];

	// Four zero crossings, e:event, i:interval
	double *nELocations, *pELocations, *dnELocations, *dpELocations;
	double *nILocations, *pILocations, *dnILocations, *dpILocations;
	double *nIntervals, *pIntervals, *dnIntervals, *dpIntervals;
	int nLen, pLen, dnLen, dpLen;
	nELocations = (double *)malloc(sizeof(double) * xLength); // xLength is insurance
	pELocations = (double *)malloc(sizeof(double) * xLength);
	dnELocations = (double *)malloc(sizeof(double) * xLength);
	dpELocations = (double *)malloc(sizeof(double) * xLength);
	nILocations = (double *)malloc(sizeof(double) * xLength);
	pILocations = (double *)malloc(sizeof(double) * xLength);
	dnILocations = (double *)malloc(sizeof(double) * xLength);
	dpILocations = (double *)malloc(sizeof(double) * xLength);
	nIntervals = (double *)malloc(sizeof(double) * xLength);
	pIntervals = (double *)malloc(sizeof(double) * xLength);
	dnIntervals = (double *)malloc(sizeof(double) * xLength);
	dpIntervals = (double *)malloc(sizeof(double) * xLength);

	zeroCrossingEngine(equivalentFIR, xLength, sample_rate, 
		nELocations, nILocations, nIntervals, &nLen);

	for(i = 0;i < xLength;i++) equivalentFIR[i] = -equivalentFIR[i];
	zeroCrossingEngine(equivalentFIR, xLength, sample_rate, 
		pELocations, pILocations, pIntervals, &pLen);

	for(i = 0;i < xLength-1;i++) equivalentFIR[i] = equivalentFIR[i]-equivalentFIR[i+1];
	zeroCrossingEngine(equivalentFIR, xLength-1, sample_rate, 
		dnELocations, dnILocations, dnIntervals, &dnLen);

	for(i = 0;i < xLength-1;i++) equivalentFIR[i] = -equivalentFIR[i];
	zeroCrossingEngine(equivalentFIR, xLength-1, sample_rate, 
		dpELocations, dpILocations, dpIntervals, &dpLen);

	int usableChannel;
	usableChannel = checkEvent(nLen-2) * checkEvent(pLen-2) * 
		checkEvent(dnLen-2) * checkEvent(dpLen-2);

	double *interpolatedF0Set[4];
	if(usableChannel <= 0) 
	{ // Finish with the no-candidates.
		for(i = 0;i < num_frames;i++)
		{
			f0Deviations[i] = 100000.0;
			interpolatedF0[i] = 0.0;
		}
	}
	else
	{
		for(i = 0;i < 4;i++)
			interpolatedF0Set[i] = (double *)malloc(sizeof(double) * num_frames);
		// Four zero crossings.
		interp1(nILocations , nIntervals , nLen , time_axis, num_frames, interpolatedF0Set[0]);
		interp1(pILocations , pIntervals , pLen , time_axis, num_frames, interpolatedF0Set[1]);
		interp1(dnILocations, dnIntervals, dnLen, time_axis, num_frames, interpolatedF0Set[2]);
		interp1(dpILocations, dpIntervals, dpLen, time_axis, num_frames, interpolatedF0Set[3]);

		for(i = 0;i < num_frames;i++)
		{
			interpolatedF0[i] = (interpolatedF0Set[0][i] + interpolatedF0Set[1][i] + 
				interpolatedF0Set[2][i] + interpolatedF0Set[3][i]) / 4.0;

			f0Deviations[i]   = sqrt( ((interpolatedF0Set[0][i]-interpolatedF0[i])*(interpolatedF0Set[0][i]-interpolatedF0[i])
				+ (interpolatedF0Set[1][i]-interpolatedF0[i])*(interpolatedF0Set[1][i]-interpolatedF0[i])
				+ (interpolatedF0Set[2][i]-interpolatedF0[i])*(interpolatedF0Set[2][i]-interpolatedF0[i])
				+ (interpolatedF0Set[3][i]-interpolatedF0[i])*(interpolatedF0Set[3][i]-interpolatedF0[i])) / 3.0);

			if(interpolatedF0[i] > boundary_freq || interpolatedF0[i] < boundary_freq/2.0 
				|| interpolatedF0[i] > f0_ceil || interpolatedF0[i] < FLOOR_F0) // Don't deal with freqs below 70 Hz.
			{
				interpolatedF0[i] = 0.0;
				f0Deviations[i]   = 100000.0;
			}
		}

		for(i = 0;i < 4;i++) free(interpolatedF0Set[i]);
	}


	// Free up memory.
	free(nELocations); free(pELocations); free(dnELocations); free(dpELocations);
	free(nILocations); free(pILocations); free(dnILocations); free(dpILocations);
	free(nIntervals); free(pIntervals); free(dnIntervals); free(dpIntervals);
}

// Calculate the zero crossing.
void zeroCrossingEngine(double *waveform, int num_samples, double sample_rate,
						double *eLocations, double *iLocations, double *intervals, int *iLen)
{
	int i;
	int *negativeGoingPoints;
	negativeGoingPoints = (int *)malloc(sizeof(int) * num_samples);

	int tmp1, tmp2;
	for(i = 0;i < num_samples-1;i++) // No use in calculating the remaining points each time.
	{
		tmp1 = waveform[i] * waveform[i+1] < 0 ? 1 : 0;
		tmp2 = waveform[i+1] < waveform[i]   ? 1 : 0;
		negativeGoingPoints[i] = (i+1) * tmp1 * tmp2;
	}
	negativeGoingPoints[num_samples-1] = 0;

	// Detect the valid events.
	int *edges;
	edges = (int *)malloc(sizeof(int) * num_samples);
	int count = 0;
	for(i = 0;i < num_samples;i++)
	{
		if(negativeGoingPoints[i] > 0) edges[count++] = negativeGoingPoints[i];
	}
	// Prepare to calculate the final return value.
	double *fineEdges;
	fineEdges = (double *)malloc(sizeof(double) * count);
	for(i = 0;i < count;i++)
	{
		fineEdges[i] = (double)edges[i] - waveform[edges[i]-1]/(waveform[edges[i]]-waveform[edges[i]-1]);
	}

	*iLen = count-1;
	for(i = 0;i < *iLen;i++)
	{
		intervals[i] = sample_rate / (fineEdges[i+1] - fineEdges[i]);
		iLocations[i] = (fineEdges[i]+fineEdges[i+1])/2.0/sample_rate;
		eLocations[i] = fineEdges[i]/sample_rate;
	}
	// The 0 corner case.
	if(count != 0) eLocations[count-1] = fineEdges[count-1]/sample_rate;

	free(fineEdges);
	free(edges);
	free(negativeGoingPoints);
}

// Nuttall window.  It looks like a magic equation but this is correct.
void nuttallWindow(int target_num_samples, double *y)
{
	int i;
	double tmp;
	for(i = 0;i < target_num_samples;i++)
	{
		tmp  = ((double)(i+1) - (double)(target_num_samples+1)/2.0) / (double)(target_num_samples+1);
		y[i] = 0.355768+0.487396*cos(2*PI*tmp)+0.144232*cos(4*PI*tmp)+0.012604*cos(6*PI*tmp);
	}
}
