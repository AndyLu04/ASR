void offline_ASR(int sampling_rate, int channels, int data_len);
double** read_data(char* file_name, int row_size, int column_size);
void write_to_file(char file_name[], double** data, int row_size, int column_size);
void dft_r2c_1d(int N);
