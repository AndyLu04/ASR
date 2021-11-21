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
    find_clean_ASR_return_val  return_val = find_clean_ASR(the_ASR, data);
    double** X = return_val.X;
    int S = return_val.column_size;

    double window_len = 0.5;
    double min_clean_fraction = 0.25;
    double window_overlap = 0.66;
    double max_dropout_fraction = 0.1;


    double** uc_data = (double**)malloc(the_ASR->channels * sizeof(double*));
    for(int i=0; i< the_ASR->channels; i++)
        uc_data[i] = (double*)malloc(S * sizeof(double));
    double val[S];
    double tmp[S];
    for(int i=0; i<the_ASR->channels; i++)
    {
        for(int j=0; j<S; j++)
        {
            val[j] = X[i][j];
        }
        filter(8, the_ASR->filter_A, the_ASR->filter_B, S-1, val, tmp);
        for(int k=0; k<S; k++)
        {
            uc_data[i][k] = tmp[k];
        }
    }

    double** M = covInASR(the_ASR, the_ASR->channels, S, uc_data);

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

    int m[44];
    double j=0.578;
    for(int i=0; i<44; i++, j-=step_sizes[1])
    {
        m[i] = round(j*n);
    }
    int new_size = remove_duplicated(m, 44); // now m is a new_size big array with no duplicated element

    qsort(m, (size_t)44, sizeof(int), icomp);

    double opt_val = MAX1;
    double opt_beta = 0;
    double opt_bounds[2];
    double opt_lu[2];

    for(int i=0; i<new_size; i++)
    {
        int nbins = round(3 * log2(1 + m[i]/2));

        double divide[11];
        for(int j=0; j<11; j++)
        {
            divide[j] = nbins/new_X[m[i]-1][j];
        }

        int the_m = m[i];
        double H[the_m][11];
        for(int j=0; j<the_m; j++)
        {
            for(int k=0; k<11; k++)
            {
               H[j][k] = 0;
            }
        }

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

        for(int j=0; j<the_m; j++)
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

//        printf("the_m: %d\nopt_val: %f\nopt_beta: %f\nopt_bounds: %f %f\nopt_lu: %f %f\n\n", the_m, opt_val, opt_beta, opt_bounds[0], opt_bounds[1], opt_lu[0], opt_lu[1]);
    }

    double alpha = (opt_lu[1] - opt_lu[0]) / (opt_bounds[1] - opt_bounds[0]);
    static double* mu_and_sig;
    mu_and_sig = (double*)malloc(2 * sizeof(double));
    mu_and_sig[0] = opt_lu[0] - opt_bounds[0]*alpha;
    double beta_val = opt_beta;

    mu_and_sig[1] = sqrt(pow(alpha,2) * tgamma(3/beta_val) / tgamma(1/beta_val));

    return mu_and_sig;
}

double** covInASR(ASR_PSW* the_ASR, int C, int S, double** data)
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

    int* range = (int*)malloc(length * sizeof(int));
    for(int k=0; k<blocksize; k++)
    {
        length = 0;
        for(int i=1; i<=S+k-1; i+=blocksize)
        {
            length += 1;
        }
        for(int i=0; i<length; i++)
        {
            if(S > ((k + 1) + i*10))
            {
                range[i] = (k + 1) + i*10;
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
                int column = j%19;
                int hight = j/19;
                U[i][j] = U[i][j] + temp2[i][column][hight];
            }
        }

        for(int i=0; i<length; i++)
        {
            for(int j=0; j<the_ASR->channels; j++)
            {
                free(temp2[i][j]);
            }
            free(temp2[i]);
        }
        free(temp2);
    }

    for(int i=0; i<length; i++)
    {
        for(int j=0; j<C*C; j++)
        {
            U[i][j] = U[i][j] / blocksize;
        }
    }

    double* y = block_geometric_median(U, length, C*C);


    printf("asd");
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

    return y;
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

void filter(int ord, double *a, double *b, int np, double *x, double *y)
{
    int i,j;
    y[0]=b[0]*x[0];
    for (i=1;i<ord+1;i++)
    {
        y[i]=0.0;
        for (j=0;j<i+1;j++)
            y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<i;j++)
            y[i]=y[i]-a[j+1]*y[i-j-1];
    }
    /* end of initial part */
    for (i=ord+1;i<np+1;i++)
    {
        y[i]=0.0;
        for (j=0;j<ord+1;j++)
            y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<ord;j++)
            y[i]=y[i]-a[j+1]*y[i-j-1];
    }
    return;
} /* end of filter */
