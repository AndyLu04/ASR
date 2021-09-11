#include "fftw3.h"
#include <float.h>
fftw_complex* fft_1d(int N, double* in);

double** data_filtering(double** data, int row_size, int column_size, int sampling_rate);

double* eegfiltfft(double* data, int sampling_rate, double lowcut, double highcut, int column_size, int epochframes, int revfilt);

