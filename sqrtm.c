#include "sqrtm.h"

int dsyevr(char JOBZ, char RANGE, char UPLO, int N,
       double *A, int LDA, double VL, double VU,
       int IL, int IU, double ABSTOL, int *M,
       double *W, double *Z, int LDZ, int *ISUPPZ,
       double *WORK, int LWORK, int *IWORK, int LIWORK)
{
  extern void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                      double *A, int *LDAp, double *VLp, double *VUp,
                      int *ILp, int *IUp, double *ABSTOLp, int *Mp,
                      double *W, double *Z, int *LDZp, int *ISUPPZ,
                      double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                      int *INFOp);
  int INFO;
  dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
          &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
          WORK, &LWORK, IWORK, &LIWORK, &INFO);
  return INFO;
}

void set_entry(double *A, int i, int j, double val, int N)
{
    A[j*N+i] = val;
}

double get_entry(const double *A, int i, int j, int N)
{
  return A[j*N+i];
}

static double dlamch(char CMACH)
{
  extern double dlamch_(char *CMACHp);
  return dlamch_(&CMACH);
}

void sqrtm(double* A, int N)
{
    double *B, *W, *Z, *WORK;
    int *ISUPPZ, *IWORK;
    int  i, j;
    int  M;

    /* allocate space for the output parameters and workspace arrays */
    W = malloc(N*sizeof(double));
    Z = malloc(N*N*sizeof(double));
    ISUPPZ = malloc(2*N*sizeof(int));
    WORK = malloc(26*N*sizeof(double));
    IWORK = malloc(10*N*sizeof(int));

    /* get the eigenvalues and eigenvectors */
    dsyevr('V', 'A', 'L', N, A, N, 0, 0, 0, 0, dlamch('S'), &M,
            W, Z, N, ISUPPZ, WORK, 26*N, IWORK, 10*N);

    /* allocate and initialise a new matrix B=Z*D */
    B = malloc(N*N*sizeof(double));
    for (j=0; j<N; ++j) {
        double  lambda=sqrt(W[j]);
        for (i=0; i<N; ++i) {
            set_entry(B, i, j, get_entry(Z,i,j, N)*lambda, N);
        }
    }

    /* calculate the square root A=B*Z^T */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, N, N,
                  1, B, N, Z, N, 0, A, N);

    /* emit the result */
//    for (i=0; i<N; ++i) {
//        for (j=0; j<N; ++j) {
//            double  x = get_entry(A, i, j, N);
//            printf("%6.2f", x);
//        }
//        putchar('\n');
//    }
//    putchar('\n');

    /* check the result by calculating A*A */
//    memcpy(B, A, N*N*sizeof(double));
//    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N,
//                  1, A, N, B, N, 0, Z, N);
//
//    for (i=0; i<N; ++i) {
//        for (j=0; j<N; ++j) {
//            double  x = get_entry(Z, i, j, N);
//            printf("%6.2f", x);
//        }
//        putchar('\n');
//    }

    free(W);
    free(Z);
    free(ISUPPZ);
    free(WORK);
    free(IWORK);
    free(B);
}

double* eigenvector(double* A, int N)
{
    double *W, *Z, *WORK;
    int *ISUPPZ, *IWORK;
    int  M;

    W = malloc(N*sizeof(double));
    Z = malloc(N*N*sizeof(double));
    ISUPPZ = malloc(2*N*sizeof(int));
    WORK = malloc(26*N*sizeof(double));
    IWORK = malloc(10*N*sizeof(int));

    dsyevr('V', 'A', 'L', N, A, N, 0, 0, 0, 0, dlamch('S'), &M,
            W, Z, N, ISUPPZ, WORK, 26*N, IWORK, 10*N);

    free(W);
    free(Z);
    free(ISUPPZ);
    free(WORK);
    free(IWORK);

    return Z;
}

