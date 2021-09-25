#include "ASR.h"

ASR_PSW create_ASR(int cutoff, int sampling_rate, int channels)
{
    ASR_PSW the_ASR;
    the_ASR.sampling_rate = sampling_rate;
    the_ASR.cutoff = cutoff;
    the_ASR.channels = channels;
    the_ASR.M = NULL;

    if((sampling_rate/2 - 1) > 80)
    {
        ;
        double filter_A_value[9] = {1, -0.991323609939396, 0.315956314546930, -0.0708347481677546, -0.0558793822071134,
                                    -0.253961902647894, 0.247305661525119, -0.0420478437473110, 0.00774557183344621};
        double filter_B_value[9] = {1.44894833258029, -2.66925147648027, 2.08139706207309, -0.973667887704926,
                                    0.105460506035271, -0.188910169231485, 0.611133163659227, -0.361648301307508, 0.183431306077666};
        for(int i=0; i<9; i++)
        {
            the_ASR.filter_A[i] = filter_A_value[i];
            the_ASR.filter_B[i] = filter_B_value[i];
        }
    }
    else
    {
        // reserved for implementing dynamic sampling rate calculate IIR
        ;
    }

    return the_ASR;
}

void update_ASR(ASR_PSW* the_ASR, double** data)
{
    if(!the_ASR->M)
    {
        subspace_ASR(the_ASR, data);
    }

}

void subspace_ASR(ASR_PSW* the_ASR, double** data)
{
    find_clean_ASR(the_ASR, data);

}

void find_clean_ASR(ASR_PSW* the_ASR, double** data)
{
    int clwin_window_len = 1;
    double window_overlap = 0.66;

    int row = the_ASR->channels;
    int length = the_ASR->data_length;
    int sampling_rate = the_ASR->sampling_rate;

    int N = clwin_window_len*the_ASR->sampling_rate;

    int wnd[N];
    for(int i=0; i<N; i++)
    {
        wnd[i] = i;
    }

    int offset_size = 1 + (the_ASR->data_length-N)/(the_ASR->sampling_rate*(1-window_overlap));
    int offsets[offset_size];

    for(int i=1, j=0; i<the_ASR->data_length-N; i+=the_ASR->sampling_rate*(1-window_overlap), j++)
    {
        offsets[j] = i;
    }

    double min_clean_fraction =0.25;
    double max_dropout_fraction = 0.1;
    double truncate_quant[2] = {0.022, 0.6};
    double clwin_step_sizes[2] = {0.01, 0.01};
    double shape_range[13] = {1.7, 1.85, 2, 2.15, 2.30, 2.45, 2.6, 2.75, 2.9, 3.05, 3.2, 3.35, 3.5};

    double X[length];
    //double** X = double_2d_array_allocate(the_ASR->channels, the_ASR->data_length);

    for(int c=the_ASR->channels-1; c>=0; c--)
    {
        for(int j=0; j<the_ASR->data_length; j++)
        {
            X[j] = data[c][j]*data[c][j];
        }

        int bsxfun_plus[sampling_rate][offset_size];
        //int** bsxfun_plus = int_2d_array_allocate(the_ASR->sampling_rate, offset_size);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<offset_size; j++)
            {
                bsxfun_plus[i][j] = offsets[j] + wnd[i];
            }
        }

        double tmp[sampling_rate][offset_size];
        //int** tmp = int_2d_array_allocate(the_ASR->sampling_rate, offset_size);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<offset_size; j++)
            {
                tmp[i][j] = X[bsxfun_plus[i][j]-1];
            }
        }

        double sum[offset_size];
        for(int i=0; i<offset_size; i++)
        {
            for(int j=0; j<N; j++)
            {
                sum[i] += tmp[j][i];
            }

            sum[i] = sqrt(sum[i]/N);
        }

        double* mu_and_sig;
        mu_and_sig = test_eeg_dist_revi(sum, offset_size, min_clean_fraction, max_dropout_fraction, truncate_quant, clwin_step_sizes, shape_range);

    }

}

double* test_eeg_dist_revi(double* X, int X_size, double min_clean_fraction, double max_dropout_fraction, double* truncate_quant, double* clwin_step_sizes, double* shape_range)
{
    double* mu_and_sig[2];
    int n = X_size;

    qsort(X, (size_t)n, sizeof(double), dcomp);

    for(int i=0; i<n; i++)
    {
        printf("%d : %f\n", i, X[i]);
    }


    return mu_and_sig;
}

int dcomp (const void * elem1, const void * elem2)
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}


