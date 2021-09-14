#include "ASR.h"

ASR_PSW create_ASR(int cutoff, int sampling_rate)
{
    ASR_PSW the_ASR;
    the_ASR.sampling_rate = sampling_rate;
    the_ASR.cutoff = cutoff;
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
{    if(!*the_ASR->M)
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
    int window_overlap = 0.66;


}
