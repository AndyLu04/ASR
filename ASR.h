#include "lib.h"

typedef struct ASR_PSW
{
    double sampling_rate;
    double cutoff;
    double filter_A[9];
    double filter_B[9];
    double* M;
    double channels;
    double data_length;
} ASR_PSW;

ASR_PSW create_ASR(int cutoff, int sampling_rate, int channels);

double* test_eeg_dist_revi(double* X, int X_size, double min_clean_fraction, double max_dropout_fraction, double* truncate_quant, double* clwin_step_sizes, double* shape_range);

int dcomp (const void * elem1, const void * elem2);
