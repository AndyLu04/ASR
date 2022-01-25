#include "lib.h"

typedef struct state
{
    double A[9];
    double B[9];
    double* cov;
    double* carry;
    double** iir;
    double* last_R;
    bool last_trivial;
} state;

typedef struct ASR_PSW
{
    int sampling_rate;
    double cutoff;
    double filter_A[9];
    double filter_B[9];
    double* M;
    double* T;
    double* fsm;
    int channels;
    int data_length;
    state the_state;
} ASR_PSW;

typedef struct find_clean_ASR_return_val
{
    double** X;
    int column_size;
} find_clean_ASR_return_val;

ASR_PSW create_ASR(int cutoff, int sampling_rate, int channels);

void update_ASR(ASR_PSW* the_ASR, double** data);

void subspace_ASR(ASR_PSW* the_ASR, double** data);

find_clean_ASR_return_val find_clean_ASR(ASR_PSW* the_ASR, double** data);

double* test_eeg_dist_revi(double* X, int X_size, double min_clean_fraction, double max_dropout_fraction, double* truncate_quant, double* clwin_step_sizes, double* shape_range);

double** covInASR(ASR_PSW* the_ASR, int C, int S, double** data);

double* block_geometric_median(double** X, int row, int column);

double* geometric_median(double** X, int row, int column);

double* test_asr_process(double* data, int size, double srate, ASR_PSW* the_ASR, double windowlen, double lookahead, double maxdims, int maxmem);

int dcomp (const void * elem1, const void * elem2);

int icomp (const void * elem1, const void * elem2);

int remove_duplicated(int* arr, int the_size);

void filter(int ord, double *a, double *b, int np, double *x, double *y, double* initial_state, double *iirstate);

double* reconstruct(ASR_PSW* the_ASR, double** data);

double** moving_average(int N, double** X, int x_row, int x_column, double* Zi);

void write_data_to_file(char file_name[], double* data, int row_size, int column_size);

double* inverse(double* A, int N);
