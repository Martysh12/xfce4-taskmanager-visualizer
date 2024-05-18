#include "fft.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define cfromreal(re) (re)
#define cfromimag(re) (re * I)
#define addcc(a, b) ((a) + (b))
#define subcc(a, b) ((a) - (b))
#define mulcc(a, b) ((a) * (b))

static float fft_in_buf[FFT_SIZE];
static float fft_window_buf[FFT_SIZE];
static float complex fft_raw[FFT_SIZE];

float fft_out[FFT_SIZE];

void fft_stream_callback(void *bufferData, unsigned int frames) {
	// Push the new values to the buffer
	float (*fs)[2] = bufferData;

	for (size_t i = 0; i < frames; i++) {
		memmove(fft_in_buf, fft_in_buf + 1, (FFT_SIZE - 1) * sizeof(fft_in_buf[0]));
		fft_in_buf[FFT_SIZE - 1] = fs[i][0];
	}
}

static void fft(float in[], float complex out[], size_t stride, size_t n) {
	if (n == 1) {
		out[0] = cfromreal(in[0]);
		return;
	}

	fft(in, out, stride*2, n/2);
	fft(in + stride, out + n/2, stride*2, n/2);

	for (size_t k = 0; k < n/2; ++k) {
		float t = (float)k/n;
		float complex v = mulcc(cexpf(cfromimag(-2*M_PI*t)), out[k + n/2]);
		float complex e = out[k];
		out[k]       = addcc(e, v);
		out[k + n/2] = subcc(e, v);
	}
}


static inline float amp(float complex z) {
	float a = crealf(z);
	float b = cimagf(z);
	return logf(a*a + b*b);
}

void fft_analyze(void) {
	// Apply the Hann Window function on the input buffer.
	for (size_t i = 0; i < FFT_SIZE; i++) {
		float t = (float) i / (FFT_SIZE - 1);
		float hann = 0.5 - 0.5 * cosf(2 * M_PI * t);

		fft_window_buf[i] = fft_in_buf[i] * hann;
	}

	// Perform FFT
	fft(fft_window_buf, fft_raw, 1, FFT_SIZE);

	// "Squash" into the Logarithmic Scale
	float step = 1.06;
	float lowf = 1.0f;
	size_t m = 0;
	float max_amp = 1.0f;
	for (float f = lowf; (size_t) f < FFT_SIZE/2; f = ceilf(f*step)) {
		float f1 = ceilf(f*step);
		float a = 0.0f;
		for (size_t q = (size_t) f; q < FFT_SIZE/2 && q < (size_t) f1; ++q) {
			float b = amp(fft_raw[q]);
			if (b > a) a = b;
		}
		if (max_amp < a) max_amp = a;
		fft_out[m++] = a;
	}

	// Normalize
	for (size_t i = 0; i < m; i++) {
		fft_out[i] /= 15;
	}
}
