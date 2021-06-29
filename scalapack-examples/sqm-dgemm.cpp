#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include "mpi.h"


using namespace std;

extern "C" {
    /*
    void pdgemm_(
            char* TRANSA, // = 'N', A is used in the computation
            char* TRANSB, // = 'N', B is used in the computation
            int* M, // numbers of rows in submatrix C
            int* N,  // number of columns in submatrix C
            int* K, // If transa = 'N', it is the number of columns in submatrix A
            double* ALPHA, // scalar
            double* A, //
            int* IA, // first row of the submatrix A
            int* JA, // first column of the submatrix A
            int* DESCA, // array descriptor for global matrix A
            double* B, 
            int* IB, 
            int* JB, 
            int* DESCB,
            double* BETA,
            double* C, 
            int* IC, 
            int* JC, 
            int* DESCC
            );

    */
    void pdgemm_(const char& transpose_a, const char& transpose_b, const int& ac_rows,
            const int& bc_columns, const int& ab_columns_rows, const double& alpha,
            const double* a_matrix, const int& ia, const int& ja, int* a_descriptor,
            const double * b_matrix, const int& ib, const int& jb, int* b_descriptor,
            const double& beta,
            double* c_matrix, const int& ic, const int& jc, int* c_descriptor);



    //void pdgemm_(char*, char*, int*, int*, int*, double*, double*, int*, int*, int*, double*, int*, int*, int*, double*, double*, int*, int*, int*);

    void blacs_get_(int*, int*, int*);
    void blacs_gridmap_(int* ctx,
            int* usermap,
            int* ldumap,
            int* nprow,
            int* npcol);
    //blacs_gridmap_(&ctx, usermap.data(), &nprow, &nprow, &npcol);
    void blacs_gridinfo_(int* ctx,
            int* nprow,
            int* npcol,
            int* myrow,
            int* mycol);


    int numroc_(int* nglobal,
            int* block_size,
            int* processindex,
            int* firstprocessrowcol,
            int* nproworcol);

    void descinit_(int* desc,
            int* m,
            int* n,
            int* mb,
            int* nb,
            int* irsrc,
            int* icsrc,
            int* ctx,
            int* lddA,
            int* info);

    void blacs_gridexit_(int* ctx);

    /*
    void pdgemm_(const char* TRANSA, 
            const char* TRANSB, 
            const int* ncrows, 
            const int* nccols,
            const int* nacols,
            const double* ALPHA, 
            const double* A, 
            const int* ia,
            const int* ja,
            int* descA,
            const double* B, 
            const int* ib,
            const int* jb,
            int* descB,
            const double *BETA,
            double* C, 
            const int* ic,
            const int* jc,
            int* descC);
    */
}
extern "C" void blacs_gridinit_(int*, char*, int*, int*);
int main(int argc, char **argv)
{
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    int nprow = 2;
    int npcol = 2;

    assert(nprow*npcol == nprocs);

    int m = 4;
    int n = 4;

    int nb = 2;
    int mb = 2;



    int zero = 0;
    int ctx = -1;

    //retrieve a defualt system context
    blacs_get_(&zero, &zero, &ctx);

    vector<int> usermap(nprow * npcol);

    //col major order of mapping
    for(int it=0; it < nprow * npcol; it++)
        usermap[it] = it;

    /*
    usermap[0] = 0;
    usermap[1] = 3;
    usermap[2] = 2;
    usermap[3] = 1;
    */
    //blacs_gridmap_(&ctx, usermap.data(), &nprow, &nprow, &npcol);
    char layout='R'; // Block cyclic, Row major processor mapping
    blacs_gridinit_(&ctx, &layout, &nprow, &npcol); // Context -> Initialize the grid

    int myrow, mycol;

    //obtain identity of myrow and mycol in the grid
    blacs_gridinfo_(&ctx, &nprow, &npcol, &myrow, &mycol);

    std::cout << "my rank = " << myrank << ", gridid=(" << myrow << ", " << mycol << ")" << endl; ;

    int firstprocessrowcol = 0;
    //number of rows for the present process
    int mlocal    = numroc_(&n, &nb, &myrow, &firstprocessrowcol, &nprow); 
    //number of cols for the present process
    int nlocal   = numroc_( &n, &nb, &mycol, &firstprocessrowcol, &npcol);


    std::cout << "Grid id=(" << myrow << ", " << mycol << "), localdatasize= (" << mlocal <<", " << nlocal << ")" << endl;


    double *A = new double [mlocal*nlocal];
    double *B = new double [mlocal*nlocal];
    double *C = new double [mlocal*nlocal];

    if(A == nullptr) std::cout << "Memory allocation error for A on proc (" << myrow << ", " << mycol << ")" << endl;
    if(B == nullptr) std::cout << "Memory allocation error for B on proc (" << myrow << ", " << mycol << ")" << endl;
    if(C == nullptr) std::cout << "Memory allocation error for C on proc (" << myrow << ", " << mycol << ")" << endl;


    for(size_t k = 0; k < mlocal * nlocal; k++) {
        A[k] = myrank + 1.0 + k;
        B[k] = myrank*myrank +1.0 +k;
        B[k] = myrank +1.0 +k;
        C[k] = 0.0;
    }
    int descA[9];
    int descB[9];
    int descC[9];

    int info;
    int ldA = mlocal>1?mlocal:1;
    int ldB = ldA;
    int ldC = ldA;

    descinit_(descA, &n, &n, &nb, &nb, &firstprocessrowcol, &firstprocessrowcol, &ctx, &ldA, &info);
    assert(info == 0);
    descinit_(descB, &n, &n, &nb, &nb, &firstprocessrowcol, &firstprocessrowcol, &ctx, &ldB, &info);
    assert(info == 0);
    descinit_(descC, &n, &n, &nb, &nb, &firstprocessrowcol, &firstprocessrowcol, &ctx, &ldC, &info);
    assert(info == 0);



    //calling pdgemm function here
    double alpha = 1.0;
    double beta = 0;


    
    int one=1;
    char notrans='N';



    //pdgemm_(&notrans, &notrans, &n, &n, &n, &alpha, A, &one, &one, descA,
    //        B, &one, &one, descB, &beta, C, &one, &one, descC);

    pdgemm_('N', 'N', n, n, n, 1, A, 1, 1, descA, B, 1, 1, descB, 0, C, 1, 1, descC);
    printf("My rank = %d, output values = %lf, %lf, %lf, %lf\n", myrank, C[0], C[1], C[2], C[3]);
    
    delete [] A;
    delete [] B;
    delete [] C;
    blacs_gridexit_(&ctx);
    MPI_Finalize();
    return 0;
}
