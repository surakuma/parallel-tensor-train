#ifndef LAPACK_HEADER_EXTERN_H_
#define LAPACK_HEADER_EXTERN_H_

extern "C" {
    extern int dgemm_(const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K,
            const double* ALPHA, const double* A, const int* LDA, const double* B,
            const int* LDB, const double* BETA, double* C, const int* LDC);
}


extern "C" {
    extern int dlaqps_(int const &     M,     // [in] The number of rows of the matrix A. M >= 0.
            int const &     N,     // [in] Number of columns 
            int const &          OFFSET, // [in] Number of columns factorized in previous steps
            int const &          NB,     // [in] number of cols to factorize
            int const &          KB,     // [out] number of columns actually factorized
            double*  A,      // [in,out] dimension(LDA,N) input matrix
            int const &          LDA,    // [in] leading dimension of A
            int*     JPVT,   // [in,out] JPVT(I) = K <==> Column K of A to position I in AP.
            double*  TAU,    // [out] dimension KB scala factor of elementary reflectors
            double*  VN1,    // [in,out] dim N partial column norms
            double*  VN2,    // [in,out] dim N exact column norm
            double*  AUXV,   // [in,out] dim NB auxiliary vector
            double*  F,      // [in,out] dim (LDF,NB) F.T = L Y.T A
            int const &          LDF     // [in] leading dimension of F
            );

    extern int dgeqp3_(int const & M,
            int const & N,
            double*  A,
            int const & LDA,
            int*  JPVT,
            double*  TAU,
            double* WORK,
            int const& LWORK,
            int & info);

    extern int dgeqpf_(int const & M,
            int const & N,
            double*  A,
            int const & LDA,
            int*  JPVT,
            double*  TAU,
            double* WORK,
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


    extern int dgesvd_(char&   JOBU, // [in] 
            char&   JOBVT, // [in]
            int& M, // [in] number of rows
            int& N, // [in] number of columns
            double *A, // [in,out] A is DOUBLE PRECISION array
            int& LDA, // [in] leading dimension of the array A
            double* S, // [out] singular values of A, sorted so that S(i) >= S(i+1)
            double *U, // [out] left singular vectors
            int& LDU, // [in] leading dimension of the array U
            double* VT, // [out] right singular vectors
            int& LDVT, // [in] leading dimension of the array VT
            double* work,
            int& lwork, // [in] To be safe >= 4* min(M,N) * min(M, N) + 7* min(M, N)
            int& info);


}
#endif //LAPACK_HEADER_EXTERN_H_

