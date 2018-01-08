#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926535897932384

#define FFT_FORWARD 1
#define FFT_BACKWARD 2
#define FFT_ESTIMATE 3

// debug
#define MAXIMUM_FFTL 4194304

typedef double fft_complex[2];
typedef struct {
	int n;
	int sign;
	unsigned int flags;
	fft_complex *c_in;
	double *in;
	fft_complex *c_out;
	double *out;
	double *input;
	int *ip;
	double *w;
} fft_plan;

fft_plan fft_plan_dft_1d(int n, fft_complex *in, fft_complex *out, int sign, unsigned int flags);
fft_plan fft_plan_dft_c2r_1d(int n, fft_complex *in, double *out, unsigned int flags);
fft_plan fft_plan_dft_r2c_1d(int n, double *in, fft_complex *out, unsigned int flags);
void fft_execute(const fft_plan p);
void fft_destroy_plan(fft_plan p);

// hidden functions
void rdft(int n, int isgn, double *a, int *ip, double *w);
void cdft(int n, int isgn, double *a, int *ip, double *w);