double* AB_filter(int sampling_rate, char type)
{
    double* filter = (double*)malloc(9 * sizeof(double));
    if(type == 'A')
    {
        switch(sampling_rate)
        {
            case 200:
                filter[0] = 1;
                filter[1] = -0.991323609939396;
                filter[2] = 0.315956314546930;
                filter[3] = -0.0708347481677546;
                filter[4] = -0.0558793822071134;
                filter[5] = -0.253961902647894;
                filter[6] = 0.247305661525119;
                filter[7] = -0.0420478437473110;
                filter[8] = 0.00774557183344621;
                break;
            case 500:
                filter[0] = 1;
                filter[1] = -4.689332908445262;
                filter[2] = 10.5989986701082;
                filter[3] = -14.9691518101371;
                filter[4] = 14.332035839974143;
                filter[5] = -9.49243170691796;
                filter[6] = 4.24258996189885;
                filter[7] = -1.17156009751804;
                filter[8] = 0.153804842771786;
                break;
        }
    }
    else if(type == 'B')
    {
        switch(sampling_rate)
        {
            case 200:
                filter[0] = 1.44894833258029;
                filter[1] = -2.66925147648027;
                filter[2] = 2.08139706207309;
                filter[3] = -0.973667887704926;
                filter[4] = 0.105460506035271;
                filter[5] = -0.188910169231485;
                filter[6] = 0.611133163659227;
                filter[7] = -0.361648301307508;
                filter[8] = 0.183431306077666;
                break;
            case 500:
                filter[0] = 2.31335200869725;
                filter[1] = -11.9471223009152;
                filter[2] = 29.1067166493392;
                filter[3] = -43.7550171007296;
                filter[4] = 44.3385767452337;
                filter[5] = -30.9965523846531;
                filter[6] = 14.6209883020843;
                filter[7] = -4.27434124003581;
                filter[8] = 0.598255358378799;
                break;
        }
    }

    return filter;
}
