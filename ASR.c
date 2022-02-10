#include "ASR.h"
#include "sqrtm.h"

ASR_PSW create_ASR(int cutoff, int sampling_rate, int channels)
{
    ASR_PSW the_ASR;
    the_ASR.sampling_rate = sampling_rate;
    the_ASR.cutoff = cutoff;
    the_ASR.channels = channels;
    the_ASR.M = NULL;
    the_ASR.T = NULL;
    the_ASR.fsm = NULL;
    the_ASR.the_state.cov = NULL;
    the_ASR.the_state.last_R = NULL;
    the_ASR.the_state.iir = (double**)malloc(channels * sizeof(double*));
    the_ASR.the_state.last_trivial = true;
    for(int i=0; i<channels; i++)
    {
        the_ASR.the_state.iir[i] = (double*)malloc(8 * sizeof(double));
    }

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

double* reconstruct(ASR_PSW* the_ASR, double** data)
{
    int availableRAM_GB = 8;
    double asr_windowlen = 0.5;
    int substract_val = round(asr_windowlen * the_ASR->sampling_rate / 2);
    int end_size = the_ASR->data_length - substract_val;
    int new_size = the_ASR->data_length + substract_val;

    double last_times_two[the_ASR->channels];
    for(int i=0; i<the_ASR->channels; i++)
    {
        last_times_two[i] = 2 * data[i][the_ASR->data_length-1];
    }

    double* sig = (double*)malloc(the_ASR->channels * new_size * sizeof(double));
    for(int i=0; i<the_ASR->channels; i++)
    {
        for(int j=0; j<the_ASR->data_length; j++)
        {
            sig[i*new_size + j] = data[i][j];
        }
        for(int j=the_ASR->data_length, k=1; j<new_size; j++, k++)
        {
            sig[i*new_size + j] = last_times_two[i] - sig[i*new_size - 1 + the_ASR->data_length - k ];
        }
    }
    double maxdims = 0.66;

    double* return_data = test_asr_process(sig, new_size, the_ASR->sampling_rate, the_ASR, asr_windowlen, asr_windowlen/2, maxdims, availableRAM_GB);

    int carry_size = (asr_windowlen/2)*the_ASR->sampling_rate;
    int new_column = new_size-carry_size;
    double* data_ASR = malloc(the_ASR->channels*new_column * sizeof(double));
    int index_dst = 0;
    int index_src = 50;
    for(int i=0; i<the_ASR->channels; i++)
    {
        memcpy(data_ASR+index_dst, return_data+index_src, new_column * sizeof(double));
        index_dst += new_column;
        index_src += new_size;
    }

    for(int i=0; i<the_ASR->channels; i++)
        free(data[i]);
    free(data);

    free(return_data);

    return data_ASR;
}

void subspace_ASR(ASR_PSW* the_ASR, double** data)
{
    find_clean_ASR_return_val  return_val = find_clean_ASR(the_ASR, data);
    double** X = return_val.X;
    int S = return_val.column_size;

    double window_len = 0.5;
    double min_clean_fraction = 0.25;
    double window_overlap = 0.66;
    double max_dropout_fraction = 0.1;


    double** uc_data = (double**)malloc(the_ASR->channels * sizeof(double*));
    for(int i=0; i<the_ASR->channels; i++)
        uc_data[i] = (double*)malloc(S * sizeof(double));
    double** X_transpose = (double**)malloc(S * sizeof(double*));
    for(int i=0; i<S; i++)
        X_transpose[i] = (double*)malloc(the_ASR->channels * sizeof(double));
    double val[S];
    double tmp[S];
    for(int i=0; i<the_ASR->channels; i++)
    {
        for(int j=0; j<S; j++)
        {
            val[j] = X[i][j];
        }
        filter(8, the_ASR->filter_A, the_ASR->filter_B, S-1, val, tmp, NULL, the_ASR->the_state.iir[i]);
        for(int k=0; k<S; k++)
        {
            uc_data[i][k] = tmp[k];
            X_transpose[k][i] = uc_data[i][k];
        }
    }
    for(int i=0; i<the_ASR->channels; i++)
        free(X[i]);
    free(X);

    the_ASR->M = covInASR(the_ASR, the_ASR->channels, S, uc_data);
    for(int i=0; i<the_ASR->channels; i++)
        free(uc_data[i]);
    free(uc_data);


    double* temp = (double*)malloc(the_ASR->channels*the_ASR->channels * sizeof(double)); // use temp because after calculate eigenvector, the array value will change
    memcpy(temp, the_ASR->M, the_ASR->channels*the_ASR->channels*sizeof(double));
    int N = window_len * the_ASR->sampling_rate;
    eigen* the_eigen = eigenvector(temp, the_ASR->channels);
    free(temp);
    double* V = the_eigen->eigen_vector; // V's row & column are inversed, -> V is V' in matlab

    double p_n[19] = {true, false, false, true, true, true, true, true, false, true, true, false, false, false, false, true, false, false, true};

    // for making P/N equal to matlab result
//    for(int i=0; i<the_ASR->channels; i++)
//    {
//        if((p_n[i] && (V[i*the_ASR->channels] < 0)) || (!p_n[i] && (V[i*the_ASR->channels] > 0)))
//        {
//            for(int j=0; j<the_ASR->channels; j++)
//            {
//                V[i*the_ASR->channels + j] = -V[i*the_ASR->channels + j];
//            }
//        }
//    }

    double** new_X = (double**)malloc(S * sizeof(double*));
    for(int i=0; i<S; i++)
        new_X[i] = (double*)malloc(the_ASR->channels * sizeof(double));

    for(int i=0; i<S; i++)
    {
        for(int j=0; j<the_ASR->channels; j++)
        {
            new_X[i][j] = 0;
            for(int k=0; k<the_ASR->channels; k++)
            {
                new_X[i][j] += X_transpose[i][k] * V[j*the_ASR->channels+k];
            }
            new_X[i][j] = fabs(new_X[i][j]);
        }
    }

    for(int i=0; i<S; i++)
    {
        free(X_transpose[i]);
    }
    free(X_transpose);

    double mu[the_ASR->channels];
    double sig[the_ASR->channels];

    for(int c=the_ASR->channels-1; c>=0; c--)
    {
        double* rms = (double*)malloc(S * sizeof(double));
        for(int i=0; i<S; i++)
        {
            rms[i] = new_X[i][c] * new_X[i][c];
        }

        int size = floor((S-N-1)/(N*(1-window_overlap)) + 1);
        int diff = N*(1-window_overlap);
        int tmp1[size];
        for(int i=0; i<size; i++)
        {
            tmp1[i] = 1 + diff*i;
        }

        double tmp2[size];
        for(int i=0; i<size; i++)
        {
            double sum = 0;
            for(int j=0; j<N; j++)
            {
                sum += rms[j + tmp1[i]-1];
            }
            tmp2[i] = sqrt(sum/N);
        }

        double* mu_and_sig = test_eeg_dist_revi(tmp2, size, min_clean_fraction, max_dropout_fraction, NULL, NULL, NULL);
        mu[c] = mu_and_sig[0];
        sig[c] = mu_and_sig[1];

        free(rms);
    }

    the_ASR->T = (double*)malloc(the_ASR->channels*the_ASR->channels * sizeof(double));

    for(int i=0; i<the_ASR->channels; i++)
    {
        mu[i] = mu[i] + the_ASR->cutoff*sig[i];
    }

    for(int i=0; i<the_ASR->channels; i++)
    {
        for(int j=0; j<the_ASR->channels; j++)
        {
            the_ASR->T[i*the_ASR->channels+j] = mu[i]*V[i*the_ASR->channels + j];
        }
    }

//    if(the_ASR->fsm == NULL)
//    {
//        double** Y_0 = (double**)malloc(the_ASR->channels * sizeof(double*));
//        for(int i=0; i<the_ASR->channels; i++)
//            Y_0[i] = (double*)malloc(S * sizeof(double));
//        for(int i=0; i<the_ASR->channels; i++)
//        {
//            for(int j=0; j<S; j++)
//            {
//                Y_0[i][j] = 0;
//                for(int k=0; k<the_ASR->channels; k++)
//                {
//                    Y_0[i][j] += V[k*the_ASR->channels + i] * uc_data[k][j];
//                }
//                printf("asd");
//            }
//        }
//        printf("asd");
//    }

    for(int i=0; i<S; i++)
        free(new_X[i]);
    free(new_X);

    free(the_eigen->eigen_value);
    free(the_eigen->eigen_vector);
}

find_clean_ASR_return_val find_clean_ASR(ASR_PSW* the_ASR, double** data)
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

    double wz[the_ASR->channels][offset_size];

    for(int c=the_ASR->channels-1; c>=0; c--)
    {
        for(int j=0; j<the_ASR->data_length; j++)
        {
            X[j] = data[c][j]*data[c][j];
        }

        int bsxfun_plus[sampling_rate][offset_size];
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<offset_size; j++)
            {
                bsxfun_plus[i][j] = offsets[j] + wnd[i];
            }
        }

        double tmp[sampling_rate][offset_size];
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
            sum[i] = 0;
        }
        for(int i=0; i<offset_size; i++)
        {
            for(int j=0; j<N; j++)
            {
                sum[i] += tmp[j][i];
            }

            sum[i] = sqrt(sum[i]/N);
        }

        double* mu_and_sig = test_eeg_dist_revi(sum, offset_size, min_clean_fraction, max_dropout_fraction, truncate_quant, clwin_step_sizes, shape_range);

        for(int i=0; i<offset_size; i++)
        {
            wz[c][i] = (sum[i] - mu_and_sig[0]) / mu_and_sig[1];
        }

    }

    int max_bad_channels = round(row*0.2);
    double zthresholds[2] = {-3.5, 5};
    double swz[the_ASR->channels][offset_size];
    double tmp2[the_ASR->channels];

    for(int i=0; i<offset_size; i++)
    {
        for(int j=0; j<the_ASR->channels; j++)
        {
            tmp2[j] = wz[j][i];
        }

        qsort(tmp2, (size_t)the_ASR->channels, sizeof(double), dcomp);
        for(int j=0; j<the_ASR->channels; j++)
        {
            swz[j][i] = tmp2[j];
        }
    }

    bool remove_mask[offset_size];
    for(int i=0; i<offset_size; i++)
    {
        remove_mask[i] = false;
    }

    for(int i=0; i<offset_size; i++)
    {
        if(swz[the_ASR->channels - max_bad_channels - 1][i] > zthresholds[1])
        {
            remove_mask[i] = true;
        }
        if(swz[1 + max_bad_channels - 1][i] < zthresholds[0])
        {
            remove_mask[i] = true;
        }
    }

    int count =0;
    for(int i=0; i<offset_size; i++)
    {
        if(remove_mask[i])
        {
            count += 1;
        }
    }
    int removed_windows[count];
    for(int i=0, j=0; i<offset_size; i++)
    {
        if(remove_mask[i])
        {
            removed_windows[j] = i;
            j += 1;
        }
    }

    int offset_val[count];
    for(int i=0; i<count; i++)
    {
        offset_val[i] = offsets[removed_windows[i]];
    }

    bool sample_mask[the_ASR->data_length];
    for(int i=0; i<the_ASR->data_length; i++)
    {
        sample_mask[i] = true;
    }
    for(int i=0; i<count; i++)
    {
        for(int j=0; j<N; j++)
        {
            sample_mask[offset_val[i] + wnd[j]-1] = false;
        }
    }

    count =0;
    for(int i=0; i<the_ASR->data_length; i++)
    {
        if(sample_mask[i])
        {
            count += 1;
        }
    }

    static double** data_clean;
    data_clean = (double **)malloc(the_ASR->channels * sizeof(double *));
    for (int i=0; i<the_ASR->channels; i++)
         data_clean[i] = (double *)malloc(count * sizeof(double));

    for(int i=0; i<the_ASR->channels; i++)
    {
        int column = 0;
        for(int j=0; j<the_ASR->data_length; j++)
        {
            if(sample_mask[j])
            {
                data_clean[i][column] = data[i][j];
                column += 1;
            }
        }
    }

    find_clean_ASR_return_val return_val;
    return_val.X = data_clean;
    return_val.column_size = count;

    return return_val;
}

double* test_eeg_dist_revi(double* origin_X, int X_size, double min_clean_fraction, double max_dropout_fraction, double* quants, double* step_sizes, double* beta)
{
    int MAX1 = 0x7fffffff;
    if(min_clean_fraction == -777)
        min_clean_fraction = 0.25;
    if(max_dropout_fraction == -777)
        max_dropout_fraction = 0.1;
    if(quants == NULL)
    {
        quants = (double*)malloc(2 * sizeof(double));
        quants[0] = 0.022;
        quants[1] = 0.6;
    }
    if(step_sizes == NULL)
    {
        step_sizes = (double*)malloc(2 * sizeof(double));
        step_sizes[0] = 0.01;
        step_sizes[1] = 0.01;
    }
    if(beta == NULL)
    {
        beta = (double*)malloc(13 * sizeof(double));
        for(int i=0; i<13; i++)
        {
            beta[i] = 1.7 + 0.15*i;
        }
    }

    double X[X_size];
    for(int i=0; i<X_size; i++)
    {
        X[i] = origin_X[i];
    }
    int n = X_size;

    qsort(X, (size_t)n, sizeof(double), dcomp);

    // b_t is the zbounds
    double b_t[13][2] = {{-1.60862458455965,0.182115306908469},
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
            new_X[i][j] = 0;
        }
    }
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

    int temp_m[44];
    double j=0.578;
    for(int i=0; i<44; i++, j-=step_sizes[1])
    {
        temp_m[i] = round(j*n);
    }
    int new_size = remove_duplicated(temp_m, 44); // now m is a new_size big array with no duplicated element

    qsort(temp_m, (size_t)new_size, sizeof(int), icomp);

    int m[new_size];
    memcpy(m, temp_m, new_size*sizeof(int));

    double opt_val = MAX1;
    double opt_beta = 0;
    double opt_bounds[2];
    double opt_lu[2];
    int the_m = 0;
    for(int i=0; i<new_size; i++)
    {
        the_m = m[i];

        int nbins = round(3 * log2(1 + m[i]/(double)2));

        double divide[11];
        for(int j=0; j<11; j++)
        {
            divide[j] = nbins/new_X[m[i]-1][j];
        }

        double H[the_m][11];
        for(int j=0; j<the_m; j++)
        {
            for(int k=0; k<11; k++)
            {
                if(new_X[j][k] != 0 && divide[k] != 0)
                {
                    H[j][k] = new_X[j][k] * divide[k];
                }
                else
                {
                    H[j][k] = 0;
                }
            }
        }

        int bounds[nbins + 1];
        double counts[nbins + 1][11]; // counts is logq
        for(int j=0; j<nbins+1; j++)
        {
            for(int k=0; k<11; k++)
            {
                counts[j][k] = 0;
            }
            bounds[j] = j;
        }
        bounds[nbins] = MAX1;

        for(int j=0; j<11; j++)
        {
            for(int k=0; k<the_m; k++)
            {
                for(int l=0; l<nbins; l++)
                {
                    if(H[k][j] >= bounds[l] && H[k][j] < bounds[l+1])
                    {
                        counts[l][j] += 1;
                        break;
                    }
                }
            }
        }

        for(int j=0; j<nbins + 1; j++)
        {
            for(int k=0; k<11; k++)
            {
                counts[j][k] = log(counts[j][k] + 0.01);
            }
        }

        int size = 1;
        for(double j=0.5; j<MAX1; j++)
        {
            if((j+1) <= nbins-0.5)
            {
                size += 1;
            }
            else
            {
                break;
            }
        }
        double tmp3[size];
        for(int j=0; j<size; j++)
        {
            tmp3[j] = (0.5 + j) / nbins;
        }

        double x_m[13][size];

        for(int j=0; j<13; j++)
        {
            double diff = b_t[j][1] - b_t[j][0];
            for(int k=0; k<size; k++)
            {
                x_m[j][k] = b_t[j][0] + tmp3[k]*diff;
            }
        }

        double p_m[13][size];
        for(int j=0; j<13 ;j++)
        {
            for(int k=0; k<size; k++)
            {
                p_m[j][k] = exp(-pow(fabs(x_m[j][k]), beta[j]))*rescale[j];
            }
        }

        for(int j=0; j<13 ;j++)
        {
            double sum = 0;
            for(int k=0; k<size; k++)
            {
                sum += p_m[j][k];
            }
            for(int k=0; k<size; k++)
            {
                p_m[j][k] = p_m[j][k]/sum;
            }
        }

        double logq_m[13][size][11];
        for(int j=0; j<13; j++)
        {
            for(int k=0; k<size; k++)
            {
                for(int l=0; l<11; l++)
                {
                    logq_m[j][k][l] = counts[k][l];
                }
            }
        }

        double kl_m[13][11] = {0};
        for(int j=0; j<13; j++)
        {
            for(int k=0; k<11; k++)
            {
                for(int l=0; l<size; l++)
                {
                    kl_m[j][k] += p_m[j][l]*(log(p_m[j][l]) - logq_m[j][l][k]);
                }
                kl_m[j][k] += log(m[i]);
            }
        }

        int idx[13];
        double min_val_t[13];
        for(int j=0; j<13; j++)
        {
            min_val_t[j] = kl_m[j][0];
            idx[j] = 0;
            for(int k=1; k<11; k++)
            {
                if(dcomp(&min_val_t[j], &kl_m[j][k]) > 0)
                {
                    min_val_t[j] = kl_m[j][k];
                    idx[j] = k;
                }
            }
        }
        int idx_b = 0;
        double min_val = min_val_t[0];
        for(int j=1; j<13; j++)
        {
            if(dcomp(&min_val, &min_val_t[j]) > 0)
            {
                min_val = min_val_t[j];
                idx_b = j;
            }
        }

        if(min_val > opt_val)
        {
            continue;
        }

        opt_val = min_val;

        opt_beta = beta[idx_b];

        opt_bounds[0] = b_t[idx_b][0];
        opt_bounds[1] = b_t[idx_b][1];
        opt_lu[0] = X1[idx[idx_b]];
        opt_lu[1] = X1[idx[idx_b]] + new_X[m[i]-1][idx[idx_b]];
    }

    double alpha = (opt_lu[1] - opt_lu[0]) / (opt_bounds[1] - opt_bounds[0]);
    static double* mu_and_sig;
    mu_and_sig = (double*)malloc(2 * sizeof(double));
    mu_and_sig[0] = opt_lu[0] - opt_bounds[0]*alpha;
    double beta_val = opt_beta;

    mu_and_sig[1] = sqrt(pow(alpha,2) * tgamma(3/beta_val) / tgamma(1/beta_val));

    return mu_and_sig;
}

double* covInASR(ASR_PSW* the_ASR, int C, int S, double** data)
{
    int blocksize = 10;
    int length = 1;
    for(int i=1; i<S; i+=10)
    {
        if((i+10) <= S)
        {
           length += 1;
        }
    }

    double ** U = (double**)malloc(length * sizeof(double*));
    for(int i=0; i<length; i++)
        U[i] = (double*)malloc(C*C * sizeof(double));
    for(int i=0; i<length; i++)
    {
        for(int j=0; j<C*C; j++)
        {
            U[i][j] = 0;
        }
    }

    for(int k=1; k<=blocksize; k++)
    {
        length = 0;
        for(int i=k; i<=S+k-1; i+=blocksize)
        {
            length += 1;
        }

        int* range = (int*)malloc(length * sizeof(int));
        for(int i=0; i<length; i++)
        {
            if(S > (k + i*10))
            {
                range[i] = k + i*10;
            }
            else
            {
                range[i] = S;
            }
        }

        double** temp = (double**)malloc(length * sizeof(double*));
        for(int i=0; i<length; i++)
            temp[i] = (double*)malloc(the_ASR->channels * sizeof(double));
        for(int i=0; i<length; i++)
        {
            for(int j=0; j<the_ASR->channels; j++)
            {
                temp[i][j] = data[j][range[i]-1];
            }
        }

        double*** temp2 = (double***)malloc(length * sizeof(double**));
        for(int i=0; i<length; i++)
        {
            temp2[i] = (double**)malloc(the_ASR->channels * sizeof(double*));
            for(int j=0; j<the_ASR->channels; j++)
            {
                temp2[i][j] = (double*)malloc(the_ASR->channels * sizeof(double));
            }
        }
        for(int d3=0; d3<the_ASR->channels; d3++)   // d1:row, d2:column, d3:hight
        {
            for(int d1=0; d1<length; d1++)
            {
                for(int d2=0; d2<the_ASR->channels; d2++)
                {
                    temp2[d1][d2][d3] = temp[d1][d3] * temp[d1][d2]; //bsxfun(@times,reshape(X(range,:),[],1,C),reshape(X(range,:),[],C,1))
                }
            }
        }

        for(int i=0; i<length; i++)
        {
            for(int j=0; j<C*C; j++)
            {
                int column = j%C;
                int hight = j/C;
                U[i][j] = U[i][j] + temp2[i][column][hight];
            }
        }

        for(int i=0; i<length; i++)
        {
            free(temp[i]);
        }
        free(temp);

        for(int i=0; i<length; i++)
        {
            for(int j=0; j<the_ASR->channels; j++)
            {
                free(temp2[i][j]);
            }
            free(temp2[i]);
        }
        free(temp2);

        free(range);
    }

    for(int i=0; i<length; i++)
    {
        for(int j=0; j<C*C; j++)
        {
            U[i][j] = U[i][j] / blocksize;
        }
    }

    double* y = block_geometric_median(U, length, C*C);

    sqrtm(y, C);

    free(U);

    return y;
}

double* block_geometric_median(double** X, int row, int column)
{
    int blocksize = 1;

    return geometric_median(X, row, column);
}

double* geometric_median(double** X, int row, int column)
{
    double tol = 0.00001;
    double* y = (double*)malloc(column * sizeof(double)); // y is the median array
    double* old_y = (double*)malloc(column * sizeof(double));
    double* temp = (double*)malloc(row * sizeof(double));
    for(int i=0; i<column; i++)
    {
        for(int j=0; j<row; j++)
        {
            temp[j] = X[j][i];
        }
        qsort(temp, (size_t)row, sizeof(double), dcomp);
        if(row%2 != 0)
        {
            y[i] = temp[row/2];
        }
        else
        {
            y[i] = (temp[row/2] + temp[row/2 - 1])/2;
        }
    }
    free(temp);

    int max_iter = 500;
    double invnorms[row];
    double sum = 0;
    double diff[column];
    double norm_diff = 0;
    double norm_y = 0;
    for(int i=0; i<max_iter; i++)
    {
        for(int j=0; j<row; j++)
        {
            invnorms[j] = 0;
            for(int k=0; k<column; k++)
            {
                invnorms[j] += pow((X[j][k] - y[k]), 2);
            }
            invnorms[j] = 1/sqrt(invnorms[j]);
        }

        sum = 0;
        for(int j=0; j<row; j++)
        {
            sum += invnorms[j];
        }

        for(int j=0; j<column; j++)
        {
            old_y[j] = y[j];
        }

        for(int j=0; j<column; j++)
        {
            y[j] = 0;
            for(int k=0; k<row; k++)
            {
                y[j] += invnorms[k] * X[k][j];
            }
            y[j] /= sum;
        }

        norm_diff = 0;
        for(int j=0; j<column; j++)
        {
            norm_diff += pow(y[j] - old_y[j], 2);
        }
        norm_diff = sqrt(norm_diff);

        norm_y = 0;
        for(int j=0; j<column; j++)
        {
            norm_y += pow(y[j], 2);
        }
        norm_y = sqrt(norm_y);

        if((norm_diff / norm_y) < tol)
            break;
    }

    free(old_y);

    return y;
}

double* test_asr_process(double* data, int size, double srate, ASR_PSW* the_ASR, double windowlen, double lookahead, double maxdims, int maxmem)
{
    maxmem = maxmem * pow(2, 30) / pow(2, 21);
    if(windowlen < 1.5*the_ASR->channels/the_ASR->sampling_rate)
    {
        windowlen = 1.5*the_ASR->channels/the_ASR->sampling_rate;
    }
    int stepsize = 32;
    bool usegpu = false;
    if(maxdims < 1)
    {
        maxdims = round(the_ASR->channels*maxdims);
    }

    int C = the_ASR->channels;
    int S = size;
    int N = round(windowlen*srate);
    int P = round(lookahead*srate);
    the_ASR->the_state.carry = (double*)malloc(the_ASR->channels * P * sizeof(double));
    for(int i=0; i<the_ASR->channels; i++)
    {
        for(int j=0; j<P; j++)
        {
            the_ASR->the_state.carry[i*P + j] = 2*data[i*S] - data[i*S + P - j];
        }
    }

    int new_column = S + P;
    double* new_data = (double*)malloc(the_ASR->channels*new_column * sizeof(double));
    for(int i=0; i<the_ASR->channels; i++)
    {
        for(int j=0; j<P; j++)
        {
            new_data[i*new_column+j] = the_ASR->the_state.carry[i*P + j];
        }
        for(int k=0; k<S; k++)
        {
            new_data[i*new_column+P+k] = data[i*size+k];
        }
    }

    double a = C*C*S*8*8 + C*C*8*S/stepsize + C*S*8*2 + S*8*5;
    double b = (double)(maxmem*1024 - (double)C*C*P*8*3/1024);
    double splits = a/b;
    splits = ceil(splits/1024);

    for(int i1=0; i1<splits; i1++)
    {
        int end = 0;
        if(S <= floor(i1*S/splits))
        {
            end = S;
        }
        else
        {
            end = floor(i1*S/splits);
        }
        int start = 1+floor((i1-1)*S/splits);

        double** filter_X = (double**)malloc(C * sizeof(double));
        for(int i=0; i<C; i++)
        {
            filter_X[i] = (double*)malloc(S * sizeof(double));
        }
        for(int i=0; i<the_ASR->channels; i++)
        {
            double* val = (double*)malloc(S * sizeof(double));
            for(int j=0; j<S; j++)
            {
                val[j] = new_data[i*new_column+j+P];
            }
            filter(8, the_ASR->filter_A, the_ASR->filter_B, S-1, val, filter_X[i], the_ASR->the_state.iir[i], the_ASR->the_state.iir[i]);
            free(val);
        }
        double** input_x = (double**)malloc(C*C * sizeof(double*));
        for(int i=0; i<C*C; i++)
        {
            input_x[i] = (double*)malloc(S * sizeof(double));
        }

        for(int i=0; i<C; i++)
        {
            for(int j=0; j<C; j++)
            {
                for(int k=0; k<S; k++)
                {
                    input_x[i*C + j][k] = filter_X[j][k] * filter_X[i][k];
                }
            }
        }
        free(filter_X);

        double** return_val = moving_average(N, input_x, C*C, S, the_ASR->the_state.cov);
        free(input_x);

        double* Xcov = return_val[0];
        the_ASR->the_state.cov = return_val[1];

        int size = floor((S+stepsize-1)/32);
        int update_at_temp[size];
        for(int i=0, j=32; i<size; i++, j+=32)
        {
            update_at_temp[i] = j>S? S-1 : j-1;
        }

        int len_update_at = 0;
        int* update_at = NULL;
        if(the_ASR->the_state.last_R == NULL)
        {
            update_at = (int*)malloc((size+1) * sizeof(int));
            len_update_at = size+1;
            int one[1] = {0};
            memcpy(update_at, one, sizeof(int));
            memcpy(update_at+1, update_at_temp, size*sizeof(int));
            the_ASR->the_state.last_R = (double*)malloc(C*C * sizeof(double));
            for(int i=0; i<C; i++)
            {
                for(int j=0; j<C; j++)
                {
                    the_ASR->the_state.last_R[i*C + j] = i==j? 1 : 0;
                }
            }
        }
        else
        {
            update_at = (int*)malloc(size * sizeof(int));
            len_update_at = size;
            memcpy(update_at, update_at_temp, size*sizeof(int));
        }

        int last_n= 0;
        double** new_Xcov = (double**)malloc(len_update_at * sizeof(double*));
        for(int j=0; j<len_update_at; j++)
        {
            new_Xcov[j] = (double*)malloc(C*C * sizeof(double));
        }
        for(int j=0, column=0; j<len_update_at; j++)
        {
            column = update_at[j];
            for(int k=0; k<C*C; k++)
            {
                new_Xcov[j][k] = Xcov[k*S + column];
            }
        }
        free(Xcov);

        for(int j=0; j<len_update_at; j++)
        {
            double* R = (double*)malloc(C*C * sizeof(double));
            eigen* the_eigen = eigenvector(new_Xcov[j], C);
            double* V = the_eigen->eigen_vector;
            double* D = the_eigen->eigen_value;
            qsort(D, (size_t)C, sizeof(double), dcomp);
            double* temp = (double*)malloc(C*C * sizeof(double));

            for(int k=0; k<C; k++)
            {
                for(int l=0; l<C; l++)
                {
                    temp[k*C + l] = 0;
                    for(int m=0; m<C; m++)
                    {
                        temp[k*C + l] += the_ASR->T[k*C + m]*V[l*C + m];
                    }
                    temp[k*C + l] = temp[k*C + l]*temp[k*C + l];
                }
            }
            bool trivial = true;
            bool keep[C];
            for(int k=0; k<C; k++)
            {
                keep[k] = true;
            }
            for(int k=0; k<C; k++)
            {
                for(int l=1; l<C; l++)
                {
                    temp[k] += temp[l*C + k];
                }
                if(!((D[k] < temp[k]) || ((k+1)<C-maxdims)))
                {
                    trivial = false;
                    keep[k] = 0;
                }
                else
                {
                    keep[k] = 1;
                }
            }
            free(temp);

            if(!trivial)
            {
                double* temp = (double*)malloc(C*C * sizeof(double));
                for(int k=0; k<C; k++)
                {
                    for(int l=0; l<C; l++)
                    {
                        temp[k*C + l] = 0;
                        for(int m=0; m<C; m++)
                        {
                            temp[k*C + l] += V[k*C + m]*the_ASR->M[m*C + l];
                        }
                    }
                }
                for(int k=0; k<C; k++)
                {
                    for(int l=0; l<C; l++)
                    {
                        temp[k*C + l] = temp[k*C + l]*keep[k];
                    }
                }

                double* pinv = inverse(temp, C);

                for(int k=0; k<C; k++)
                {
                    for(int l=0; l<C; l++)
                    {
                        temp[k*C + l] = 0;
                        for(int m=0; m<C; m++)
                        {
                            temp[k*C + l] += the_ASR->M[k*C + m]*pinv[m*C + l];
                        }
                    }
                }
                for(int k=0; k<C; k++)
                {
                    for(int l=0; l<C; l++)
                    {
                        R[k*C + l] = 0;
                        for(int m=0; m<C; m++)
                        {
                            R[k*C + l] += temp[k*C + m]*V[m*C + l];
                        }
                    }
                }
                free(temp);
            }
            else
            {
                for(int k=0; k<C; k++)
                {
                    for(int l=0; l<C; l++)
                    {
                        if(k == l)
                        {
                            R[k*C + l] = 1;
                        }
                        else
                        {
                            R[k*C + l] = 0;
                        }
                    }
                }
            }

            int n = update_at[j];
            if(!trivial || !the_ASR->the_state.last_trivial)
            {
                int subrange_start = last_n+1;
                int subrange_end = n;
                int blend_size = n-last_n;
                double* blend = (double*)malloc(blend_size * sizeof(double));
                double pi = acos(-1);
                for(int k=0; k< blend_size; k++)
                {
                    blend[k] = (1 - cos(pi * (k+1)/blend_size)) / 2;
                }

                double temp;
                double temp2;
                double answer[C][blend_size];
                for(int k=0; k<C; k++)
                {
                    for(int l=subrange_start; l<=subrange_end; l++)
                    {
                        temp = 0;
                        temp2 = 0;
                        for(int m=0; m<C; m++)
                        {
                            temp += R[k*C + m] * new_data[m*new_column + l];
                            temp2 += the_ASR->the_state.last_R[k*C + m] * new_data[m*new_column + l];
                        }

                        answer[k][l-subrange_start] = temp*blend[l - subrange_start] + temp2*(1-blend[l - subrange_start]);
                    }
                }
                for(int k=0; k<C; k++)
                {
                    for(int l=subrange_start; l<=subrange_end; l++)
                    {
                        new_data[k*new_column + l] = answer[k][l-subrange_start];
                    }
                }
                free(blend);
            }
            last_n = n;
            if(the_ASR->the_state.last_R == NULL)
            {
                the_ASR->the_state.last_R = (double*)malloc(C*C*sizeof(double));
            }
            else
            {
                free(the_ASR->the_state.last_R);
                the_ASR->the_state.last_R = (double*)malloc(C*C*sizeof(double));
            }
            memcpy(the_ASR->the_state.last_R,R, C*C*sizeof(double));
            free(R);
            the_ASR->the_state.last_trivial = trivial;

            free(V);
            free(D);
        }
        free(update_at);
    }

    for(int i=0; i<C; i++)
    {
        for(int j=0; j<P; j++)
        {
            the_ASR->the_state.carry[i*P + j] = new_data[i*new_column + new_column-2*P+j];
        }
    }

    double* outdata = (double*)malloc(C*(new_column-P) * sizeof(double));
    unsigned long long int index_dst = 0;
    unsigned long long int index_src = 0;
    for(int i=0; i<C; i++)
    {
        memcpy(outdata+index_dst, new_data+index_src, (new_column-P) * sizeof(double));
        index_dst += (new_column-P);
        index_src += new_column;
    }

    free(new_data);
    return outdata;

}

double** moving_average(int N, double** X, int x_row, int x_column, double* Zi)
{
    int zero_size = 0;
    if(Zi == NULL)
    {
        zero_size = N;
    }
    double* Y = (double*)malloc(x_row*(zero_size+x_column) * sizeof(double));
    for(int i=0; i<x_row; i++)
    {
        for(int j=0; j<zero_size; j++)
        {
            Y[i*(zero_size+x_column)+j] = 0;
        }
        for(int j=0; j<x_column; j++)
        {
            Y[i*(zero_size+x_column)+j+zero_size] = X[i][j];
        }
    }
    int M = zero_size+x_column;
    int row = x_row;
    double* new_X = (double*)malloc(row*x_column * sizeof(double));
    double* temp = (double*)malloc(row * sizeof(double));
    for(int i=0; i<row; i++)
    {
        temp[i] = 0;
    }
    for(int i1=0, i2=0; i1<x_column; i1+=1, i2+=1)
    {
        for(int j=0; j<row; j++)
        {
            temp[j] = Y[j*(zero_size+x_column) + i2]/(-N) + temp[j];
        }
        for(int j=0; j<row; j++)
        {
            new_X[j*x_column + i1] = Y[j*(zero_size+x_column) + i2+zero_size]/N + temp[j];
            temp[j] = new_X[j*x_column + i1];
        }
    }

    double* Zf = (double*)malloc(row*N * sizeof(double));
    for(int i=0; i<row; i++)
    {
        Zf[i*N] = -(new_X[i*x_column+x_column-1]*N - Y[i*(zero_size+x_column) + zero_size+x_column-N]);
    }
    for(int i=0; i<row; i++)
    {
        for(int j=1; j<N; j++)
        {
            Zf[i*N + j] = Y[i*(zero_size+x_column) + zero_size+x_column-N+j];
        }
    }

    free(Y);
    free(temp);

    double** return_val = (double**)malloc(2 * sizeof(double*));
    return_val[0] = new_X;
    return_val[1] = Zf;

    return return_val;
}

int dcomp (const void * elem1, const void * elem2)
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int icomp (const void * elem1, const void * elem2)
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int remove_duplicated(int* arr, int the_size)
{
    for(int i=0; i<the_size; i++)
    {
        for(int j=i+1; j<the_size; j++)
        {
            /* If any duplicate found */
            if(arr[i] == arr[j])
            {
                /* Delete the current duplicate element */
                for(int k=j; k < the_size - 1; k++)
                {
                    arr[k] = arr[k + 1];
                }

                /* Decrement size after removing duplicate element */
                the_size--;

                /* If shifting of elements occur then don't increment j */
                j--;
            }
        }
    }
    return the_size;
}

void filter(int ord, double *a, double *b, int np, double *x, double *y, double* initial_state, double *iirstate)
{
    int ini_len = 0;
    if(initial_state == NULL)
    {
        ini_len = ord+1;
        y[0]=b[0]*x[0];
        for (int i=1;i<ini_len;i++)
        {
            y[i]=0.0;
            for (int j=0;j<i+1;j++)
                y[i]=y[i]+b[j]*x[i-j];
            for (int j=0;j<i;j++)
                y[i]=y[i]-a[j+1]*y[i-j-1];
        }
        /* end of initial part */
    }
    else
    {
        ini_len = ord;
        for(int i=0; i<ini_len; i++)
        {
            y[i] = 0;
            for(int j=0, k=i; j<i+1; j++, k--)
            {
                y[i] += b[j]*x[k];
            }
            for (int j=0;j<i;j++)
            {
                y[i] -= a[j+1]*y[i-j-1];
            }
            y[i] += initial_state[i];
        }
    }

    for (int i=ini_len;i<np+1;i++)
    {
        y[i]=0.0;
        for (int j=0;j<ord+1;j++)
            y[i]=y[i]+b[j]*x[i-j];
        for (int j=0;j<ord;j++)
            y[i]=y[i]-a[j+1]*y[i-j-1];
    }

    for(int i=0; i<ord; i++)
    {
        iirstate[i] = 0;
        for(int j=i, k=0; k<ord-i; j++, k++)
        {
            iirstate[i] += b[j+1]*x[np-k];
        }
        for(int j=i, k=0; k<ord-i; j++, k++)
        {
            iirstate[i] = iirstate[i] - a[j+1]*y[np-k];
        }
    }

    return;
} /* end of filter */

void write_data_to_file(char file_name[], double* data, int row_size, int column_size)
{
    FILE *fpt = fopen( file_name,"w" );
    for(int i=0; i< row_size; i++)
    {
        for(int j=0; j< column_size; j++)
        {
            fprintf(fpt,"%f,", data[i*column_size + j]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
}

double* inverse(double* data, int N)
{
    int sockfd = 0;
    sockfd = socket(AF_INET , SOCK_STREAM , 0);

    if (sockfd == -1){
        printf("Fail to create a socket.");
    }

    //socket connection
    struct sockaddr_in info;
    bzero(&info,sizeof(info));
    info.sin_family = PF_INET;

    //localhost
    info.sin_addr.s_addr = inet_addr("127.0.0.1");
    info.sin_port = htons(8000);

    int err = connect(sockfd,(struct sockaddr *)&info,sizeof(info));
    if(err==-1){
        printf("Connection error");
    }

    //Send a message to server
    int size = N*N;

    send(sockfd,&N,sizeof(int),0);

    send(sockfd,data,size*sizeof(double),0);

    char server_reply[5000] = {};
    recv(sockfd,server_reply,sizeof(server_reply),0);
    double* pinv = (double*)malloc(size * sizeof(double));
    memcpy(pinv, server_reply, size * sizeof(double));

    //printf("close Socket\n");
    close(sockfd);

    return pinv;
}
