#include "datafiltering.h"
#include "ASR.h"
#include "2d_array.h"

double** read_data(char* file_name, int row_size, int column_size)
{
    FILE *fp;
    fp = fopen(file_name, "r");

    static double** data;
    data = (double **)malloc(row_size * sizeof(double *));
    for (int i=0; i<row_size; i++)
         data[i] = (double *)malloc(column_size * sizeof(double));

    unsigned long int line_size = 1000000;
    char line[line_size];
    char *result = NULL;
    int row = 0;
    int column = 0;
    double value = 0;

    while(fgets(line, line_size, fp) != NULL) {

        result = strtok(line, ",");

        while( result != NULL ) {
            value = atof(result);
            data[row][column] = value;
            result = strtok(NULL, ",");
            column += 1;
        }
        row += 1;
        column = 0;
    }

    fclose (fp);
    return data;
}

void dft_r2c_1d(int N)
{
    int digs = DECIMAL_DIG;
    double *in;
    fftw_complex *out;
    in = (double*) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    //double** fft_input = read_data("con_emg_cell{2,1}.csv", 19, 5601);

    double** fft_input = read_data("fft_input.csv", 1, 256);
    for(int i=0; i< N; i++)
    {
        in[i] = fft_input[0][i];
    }


    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE );
    fftw_execute(p);
    FILE *fpt = fopen( "fft_r2c_result.csv","w" );
    for(int i=0; i< N; i++)
    {
        //printf("%d : %.*e %.*e\n", i, digs, out[i][0], digs, out[i][1]);
        fprintf(fpt,"%.*e+%.*ei,", digs, out[i][0], digs, out[i][1]);
    }
    fclose(fpt);


    fftw_plan p2 = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE );
    fftw_execute(p2);


    fftw_destroy_plan(p);
    fftw_destroy_plan(p2);

    FILE *fpt2 = fopen( "fft_c2r_result.csv","w" );
    for(int i=0; i< N; i++)
    {
        printf("%d : %.*e\n", i, digs, in[i]/N);
        fprintf(fpt,"%.*e,", digs, in[i]/N);
    }
    fclose(fpt2);

    fftw_free(in);
    fftw_free(out);
}

void write_to_file(char file_name[], double** data, int row_size, int column_size)
{
    FILE *fpt = fopen( file_name,"w" );
    for(int i=0; i< row_size; i++)
    {
        for(int j=0; j< column_size; j++)
        {
            fprintf(fpt,"%f,", data[i][j]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
}

int main()
{
    int sampling_rate = 200;
    int channels = 19;
    //double** test_unclean = read_data("con_emg_cell{2,1}.csv", 19, 5601);
    //double** test_clean = read_data("pure_data_cell{2,1}.csv", 19, 5601);
    //printf("%lf\n%lf\n", *test_unclean[0][0], *test_clean[0][0]);

    //test_clean = data_filtering(test_clean, 19, 5601, sampling_rate);
    //write_to_file("test_clean.csv", test_clean, 19, 5601);

    //test_unclean = data_filtering(test_unclean, 19, 5601, sampling_rate);
    //write_to_file("test_unclean.csv", test_unclean, 19, 5601);

    double** clean_data = read_data("data_clean.csv", 19, 48000);
    double** unclean_data = read_data("data_unclean.csv", 19, 48000);


    // create ASR
    ASR_PSW my_ASR = create_ASR(20, sampling_rate, channels);
    my_ASR.data_length = 48000;

    // update_ASR(&my_ASR, unclean_data);

    int seg_size = my_ASR.data_length/(20*my_ASR.sampling_rate);
    int seg[seg_size];
    for(int i=0, j=0; i<seg_size; i++, j+=(20*my_ASR.sampling_rate))
        seg[i] = j;

    int index_src;
    for(int i=0; i<seg_size-1; i++)
    {
        my_ASR.data_length = 20*my_ASR.sampling_rate+1;
        index_src = seg[i];
        double** temp_data = (double**)malloc(my_ASR.channels * sizeof(double*));
        for(int j=0; j<my_ASR.channels; j++)
        {
            temp_data[j] = (double*)malloc(my_ASR.data_length * sizeof(double));
            memcpy(temp_data[j], unclean_data[j]+index_src, my_ASR.data_length * sizeof(double));
        }
        subspace_ASR(&my_ASR, temp_data);
    }

    my_ASR.data_length = 48000;
    double* data_processed = reconstruct(&my_ASR, unclean_data);

    write_data_to_file("data_processed.csv", data_processed, my_ASR.channels, my_ASR.data_length);
    //dft_r2c_1d(256);

    return 0;
}
