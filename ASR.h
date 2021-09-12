typedef struct ASR_PSW
{
    double sampling_rate;
    double cutoff;
    double filter_A[9];
    double filter_B[9];
} ASR_PSW;

ASR_PSW create_ASR(int cutoff, int sampling_rate);



