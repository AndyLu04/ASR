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
//    int sampling_rate = 200;
//    int channels = 19;
//    int data_len = 48000;
//    int cutoff = 20;
//    double filter_A_value[9] = {1, -0.991323609939396, 0.315956314546930, -0.0708347481677546, -0.0558793822071134,
//                                    -0.253961902647894, 0.247305661525119, -0.0420478437473110, 0.00774557183344621};
//    double filter_B_value[9] = {1.44894833258029, -2.66925147648027, 2.08139706207309, -0.973667887704926,
//                                    0.105460506035271, -0.188910169231485, 0.611133163659227, -0.361648301307508, 0.183431306077666};
//    double** unclean_data = read_data("data_unclean.csv", 19, 48000);

    int sampling_rate = 500;
    int channels = 8;
    int data_len = 40000;
    int cutoff = 5;
    double filter_A_value[9] = {1, -4.689332908445262, 10.5989986701082, -14.9691518101371, 14.332035839974143,
                                -9.49243170691796, 4.24258996189885, -1.17156009751804, 0.153804842771786};
    double filter_B_value[9] = {2.31335200869725, -11.9471223009152, 29.1067166493392, -43.7550171007296, 44.3385767452337,
                                -30.9965523846531, 14.6209883020843, -4.27434124003581, 0.598255358378799};
    double** unclean_data = read_data("2_filtered_40000.csv", channels, data_len);

    //double** test_unclean = read_data("con_emg_cell{2,1}.csv", 19, 5601);
    //double** test_clean = read_data("pure_data_cell{2,1}.csv", 19, 5601);

    //test_clean = data_filtering(test_clean, 19, 5601, sampling_rate);
    //write_to_file("test_clean.csv", test_clean, 19, 5601);

    //test_unclean = data_filtering(test_unclean, 19, 5601, sampling_rate);
    //write_to_file("test_unclean.csv", test_unclean, 19, 5601);

    //double** clean_data = read_data("data_clean.csv", 19, 48000);




    printf("start create ASR\n");
    // create ASR
    ASR_PSW my_ASR = create_ASR(cutoff, sampling_rate, channels, filter_A_value, filter_B_value);
    my_ASR.data_length = data_len;

    // update_ASR(&my_ASR, unclean_data);
    int window_len = 20;
    int seg_size = my_ASR.data_length/(window_len*my_ASR.sampling_rate);
    int seg[seg_size];
    for(int i=0, j=0; i<seg_size; i++, j+=(window_len*my_ASR.sampling_rate))
        seg[i] = j;

    int index_src;
    for(int i=0; i<seg_size-1; i++)
    {
        my_ASR.data_length = window_len*my_ASR.sampling_rate+1;
        index_src = seg[i];
        double** temp_data = (double**)malloc(my_ASR.channels * sizeof(double*));
        for(int j=0; j<my_ASR.channels; j++)
        {
            temp_data[j] = (double*)malloc(my_ASR.data_length * sizeof(double));
            memcpy(temp_data[j], unclean_data[j]+index_src, my_ASR.data_length * sizeof(double));
        }
        subspace_ASR(&my_ASR, temp_data);
        printf("subspace %d\n", i);
    }

    my_ASR.data_length = data_len;
    printf("start reconstruct\n");
    double* data_processed = reconstruct(&my_ASR, unclean_data);
    printf("finished reconstruct\n");

    write_data_to_file("2_data_processed.csv", data_processed, my_ASR.channels, my_ASR.data_length);
    //dft_r2c_1d(256);

    return 0;
}
