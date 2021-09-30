int** int_2d_array_allocate(int r, int c)
{
    int* ptr = malloc(r * sizeof(int));
    for(int i=0; i<c; i++)
    {
        ptr[i] = malloc(c * sizeof(int));
    }

    return ptr;
}

double** double_2d_array_allocate(int r, int c)
{
    static double** ptr;
    ptr = (double*)malloc(r * sizeof(double*));
    for(int i=0; i<c; i++)
    {
        ptr[i] = (double*)malloc(c * sizeof(double));
    }

    return ptr;
}

