#ifndef FFT_H
#define FFT_H

#define FFT_SIZE (1<<14)

extern float fft_out[FFT_SIZE];

void fft_stream_callback(void *bufferData, unsigned int frames);

void fft_analyze(void);

#endif /* !FFT_H */
