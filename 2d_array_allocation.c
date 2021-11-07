#include "2d_array.h"

int** int_2d_array_allocate(int r, int c)
{
    static int** ptr;
    ptr = (int**)malloc(r * sizeof(int));
    for(int i=0; i<c; i++)
    {
        ptr[i] = (int*)malloc(c * sizeof(int));
    }

    return ptr;
}

double** double_2d_array_allocate(int r, int c)
{
    static double** ptr;
    ptr = (double**)malloc(r * sizeof(double*));
    for(int i=0; i<c; i++)
    {
        ptr[i] = (double*)malloc(c * sizeof(double));
    }

    return ptr;
}

