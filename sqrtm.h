#include "lib.h"

int dsyevr(char JOBZ, char RANGE, char UPLO, int N,
       double *A, int LDA, double VL, double VU,
       int IL, int IU, double ABSTOL, int *M,
       double *W, double *Z, int LDZ, int *ISUPPZ,
       double *WORK, int LWORK, int *IWORK, int LIWORK);

void set_entry(double *A, int i, int j, double val, int N);

double get_entry(const double *A, int i, int j, int N);

static double dlamch(char CMACH);

void sqrtm(double* A, int N);

double* eigenvector(double* A, int N);
