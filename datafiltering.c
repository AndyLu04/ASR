#include "datafiltering.h"

double** data_filtering(double** data, int row_size, int column_size, int sampling_rate)
{
    double** result = (double **)malloc(sizeof(double*)*row_size);
    double lowcut = 1;
    double highcut = 50;
    for(int i=0; i<row_size; i++)
    {
        result[i] = eegfiltfft(data[i], sampling_rate, lowcut, highcut, column_size, 0, 0);
    }

    return result;
}

double* eegfiltfft(double* data, int sampling_rate, double lowcut, double highcut, int column_size, int epochframes, int revfilt)
{
    int chans = 1;
    int frames = column_size;
    if(epochframes == 0)
    {
        epochframes = frames;
    }

    int epochs = (double)frames/(double)epochframes;

    double divided_frequency = (double)sampling_rate/(double)epochframes;
    // init frequency vector
    double freq_vector[epochframes];
    for(int i=0; i<epochframes; i++)
    {
        freq_vector[i] = i * divided_frequency;
    }

    // find index (idxh) closest to 1 in frequency vector
    double temp[epochframes];
    for(int i=0; i< epochframes; i++)
    {
        temp[i] = 0;
    }
    double min = 65535;
    int idxl = 0;
    int idxh = 0;
    if(lowcut != 0)
    {
        for(int i=0; i<epochframes; i++)
        {
            temp[i] = fabs(freq_vector[i] - lowcut);
            if(temp[i] < min)
            {
                min = temp[i];
                idxl = i;
            }
        }
    }
    else
    {
        idxl = 0;
    }

    // find index (idxh) closest to 50 in frequency vector
    if(highcut != 0)
    {
        min = highcut;
        for(int i=0; i<epochframes; i++)
        {
            temp[i] = fabs(freq_vector[i] - highcut);
            if(temp[i] < min)
            {
                min = temp[i];
                idxh = i;
            }
        }
    }
    else
    {
        idxh = ceil(column_size/2);
    }

    // declare & init fft variables
    double *in = data;
    fftw_complex *out;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * column_size);
    fftw_plan fft_r2c_plan = fftw_plan_dft_r2c_1d(column_size, in, out, FFTW_ESTIMATE);
    fftw_plan fft_c2r_plan = fftw_plan_dft_c2r_1d(column_size, out, in, FFTW_ESTIMATE);

    // compute fft & ifft
    for(int e=0; e<epochs; e++)
    {
        for(int c=0; c<chans; c++)
        {
            fftw_execute(fft_r2c_plan);

            if(revfilt)
            {
                for(int i=idxl+1; i<idxh; i++)
                {
                    out[i][0] = 0;
                    out[i][1] = 0;
                }
                for(int i=column_size/2; i<column_size; i++)
                {
                    out[i][0] = 0;
                    out[i][1] = 0;
                }
            }
            else
            {
                for(int i=0; i<idxl; i++)
                {
                    out[i][0] = 0;
                    out[i][1] = 0;
                }
                for(int i=column_size-idxl; i<column_size; i++)
                {
                    out[i][0] = 0;
                    out[i][1] = 0;
                }
                for(int i=idxh; i< column_size; i++)
                {
                    out[i][0] = 0;
                    out[i][1] = 0;
                }
            }

            fftw_execute(fft_c2r_plan);

            for(int i=0; i< column_size; i++)
            {
                in[i] = (double)in[i]/(double)column_size;
            }
        }
    }

    fftw_destroy_plan(fft_r2c_plan);
    fftw_destroy_plan(fft_c2r_plan);
    fftw_free(out);

    return in;
}
