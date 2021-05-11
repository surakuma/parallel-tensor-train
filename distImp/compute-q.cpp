#include <iostream>
#include <cmath>
#include <vector>
//#include "lapacke.h"
using namespace std;
extern "C" {
  extern int dgeqrf_(int const & M,
          int const & N,
          double*  A,
          int const & LDA,
          double*  TAU,
          double* WORK,
          int const& LWORK,
          int & info);


  extern int dorgqr_(int const& M, //The number of rows of the matrix Q. M >= 0
          int const& N, //The number of columns of the matrix Q. M >= N >= 0
          int const& K, //The number of elementary reflectors whose product defines the
                        //matrix Q. N >= K >= 0
          double* A,    // The matrix returned by DGEQRF in the first k columns
          int const& LDA, //The first dimension of the array A. LDA >= max(1,M)
          double * TAU, // the scalar factor of the elementary reflector H(i)
          double *WORK, // double precision array with dimension (MAX(1,LWORK))
          int const& LWORK, //LWORK >= max(1,N), for optimum performance LWORK >= N*NB
          int & info);
}

//inplace computation when m >= n
double* computeQ(double *A, int m, int n)
{
    int lda = m;
    vector<double> tau(min(m,n));
    int lwork = n*n;
    vector<double> work(lwork);
    int info;

    dgeqrf_(m, n, A, m, tau.data(), work.data(), lwork, info);
    
    dorgqr_(m, n, n, A, m, tau.data(), work.data(), lwork, info);

    return A;
}
