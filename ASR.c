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

double* test_eeg_dist_revi(double* X, int X_size, double min_clean_fraction, double max_dropout_fraction, double* quants, double* step_sizes, double* beta)
{
    double* mu_and_sig[2];
    int n = X_size;

    qsort(X, (size_t)n, sizeof(double), dcomp);

    double gama_inverse[13][2] = {{-1.60862458455965,0.182115306908469},
                                  {-1.50450959704516,0.180266821894673},
                                  {-1.42417727110355,0.179143454621292},
                                  {-1.36063558014906,0.178505967053939},
                                  {-1.30933712563468,0.178199325699606},
                                  {-1.26720913658822,0.178120305557329},
                                  {-1.23210755288029,0.178198673402093},
                                  {-1.20249492545494,0.178385807546646},
                                  {-1.17724258667216,0.178647578848539},
                                  {-1.15550487523500,0.178959768933642},
                                  {-1.13663673415458,0.179305049360544},
                                  {-1.12013828424610,0.179670948476032},
                                  {-1.10561666439142,0.180048458596905}};

    double rescale[13] = {0.560384511868011,
                          0.562928077207612,
                          0.564189583547756,
                          0.564583630220663,
                          0.564388419735723,
                          0.563793492242989,
                          0.562929662214962,
                          0.561888180766006,
                          0.560733236562082,
                          0.559510254710613,
                          0.558251494699939,
                          0.556979881664586,
                          0.555711663295697};

    double lower_min = quants[0];
    double max_width = quants[1] - quants[0];
    double min_width = min_clean_fraction*max_width;

    int tmp_size = (int)round(n*max_width);
    double tmp1[tmp_size];
    for(int i=0; i<tmp_size; i++)
    {
        tmp1[i] = i+1;
    }

    double tmp2[11] = {0.0220,0.0320,0.0420,0.0520,0.0620,0.0720,0.0820,0.0920,0.1020,0.1120,0.1220};
    for(int i=0; i<11; i++)
    {
        tmp2[i] = round(n * tmp2[i]);
    }

    int bsxfun_plus[tmp_size][11];
    for(int i=0; i<tmp_size; i++)
    {
        for(int j=0; j<11; j++)
        {
            bsxfun_plus[i][j] = tmp1[i] + tmp2[j];
        }
    }

    double new_X[tmp_size][11];
    for(int i=0; i<tmp_size; i++)
    {
        for(int j=0; j<11; j++)
        {
            new_X[i][j] = X[bsxfun_plus[i][j] - 1];
        }
    }

    double X1[11];
    for(int i=0; i<11; i++)
    {
        X1[i] = new_X[0][i];
    }

    for(int i=0; i<tmp_size; i++)
    {
        for(int j=0; j<11; j++)
        {
            new_X[i][j] = new_X[i][j] - X1[j];
        }
    }

    double m[44];
    double j=0.578;
    for(int i=0; i<44; i++, j-=step_sizes[1])
    {
        m[i] = round(j*n);
    }
    int new_size = remove_duplicated(m, 44); // now m is a new_size big array with no duplicated element

    qsort(m, (size_t)44, sizeof(double), dcomp);

    for(int i=0; i<new_size; i++)
    {
        nbins = round(3 * log2(1 + m[i]/2));

        double divide[11];
        for(int j=0; j<11; j++)
        {
            divide[j] = nbins/new_X[m[i]][j];
        }
        double** H = double_2d_array_allocate(m[i], 11);
        for(int j=0; j<m[i]; j++)
        {
            for(int k=0; k<11; k++)
            {
                H[j][k] = new_X[j][k] * divide[k];
            }
        }
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

int remove_duplicated(double* arr, int size)
{
    for(int i=0; i<size; i++)
    {
        for(int j=i+1; j<size; j++)
        {
            /* If any duplicate found */
            if(arr[i] == arr[j])
            {
                /* Delete the current duplicate element */
                for(int k=j; k < size - 1; k++)
                {
                    arr[k] = arr[k + 1];
                }

                /* Decrement size after removing duplicate element */
                size--;

                /* If shifting of elements occur then don't increment j */
                j--;
            }
        }
    }
    return size;
}


