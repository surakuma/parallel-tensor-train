#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include "mpi.h"


using namespace std;

extern "C" {

    //declaration from http://www.netlib.org/scalapack/explore-html/d6/da2/pdgemm___8c_source.html
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

    //declaration based on David's code
    /*
    void pdgemm_(const char& transpose_a, const char& transpose_b, const int& ac_rows,
            const int& bc_columns, const int& ab_columns_rows, const double& alpha,
            const double* a_matrix, const int& ia, const int& ja, int* a_descriptor,
            const double * b_matrix, const int& ib, const int& jb, int* b_descriptor,
            const double& beta,
            double* c_matrix, const int& ic, const int& jc, int* c_descriptor);
    */

    //declaration from http://www.netlib.org/scalapack/explore-html/d9/dc9/blacs__get___8c_source.html
    void blacs_get_(int *ctx, int *what, int *val);


    /*
     * http://www.netlib.org/scalapack/explore-html/de/dfb/blacs__map___8c_source.html
     */

    void blacs_gridmap_(int* ctx,
            int* usermap,
            int* ldumap,
            int* nprow,
            int* npcol);


    /*
     * http://www.netlib.org/scalapack/explore-html/d9/de2/blacs__info___8c_source.html
     */
    void blacs_gridinfo_(int* ctx,
            int* nprow,
            int* npcol,
            int* myrow,
            int* mycol);

    /*
     *
     */

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
    /*
     * http://www.netlib.org/scalapack/explore-html/db/dcd/blacs__grid___8c_source.html
     */

    void blacs_gridexit_(int* ctx);
    /*
     * http://www.netlib.org/scalapack/explore-html/d5/db4/blacs__init___8c_source.html
     */
    void blacs_gridinit_(int *Ctx, char* order, int *nprow, int *npcol);

}
int main(int argc, char **argv)
{
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    int nprow = 2;
    int npcol = 2;

    assert(nprow*npcol == nprocs);

    int m = 3;
    int k = 5;
    int n = 2;

    int mb = 2;
    int kb = 3;
    int nb = 1;



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
    int mlocalrow    = numroc_(&m, &mb, &myrow, &firstprocessrowcol, &nprow); 
    int klocalrow    = numroc_(&k, &kb, &myrow, &firstprocessrowcol, &nprow); 
    //int nlocalrow    = numroc_(&n, &nb, &myrow, &firstprocessrowcol, &nprow); 
    //number of cols for the present process
    //int mlocalcol    = numroc_(&m, &mb, &mycol, &firstprocessrowcol, &npcol); 
    int klocalcol    = numroc_(&k, &kb, &mycol, &firstprocessrowcol, &npcol); 
    int nlocalcol    = numroc_(&n, &nb, &mycol, &firstprocessrowcol, &npcol); 


    std::cout << "Grid id=(" << myrow << ", " << mycol << "), A localdatasize= (" << mlocalrow <<", " << klocalcol << ")" <<
        ", B localdatasize= (" << klocalrow <<", " << nlocalcol << ")" << 
        ", C localdatasize= (" << mlocalrow <<", " << nlocalcol << ")" << endl;


    double *A = new double [mlocalrow * klocalcol];
    double *B = new double [klocalrow * nlocalcol];
    double *C = new double [mlocalrow * nlocalcol];

    if(A == nullptr) std::cout << "Memory allocation error for A on proc (" << myrow << ", " << mycol << ")" << endl;
    if(B == nullptr) std::cout << "Memory allocation error for B on proc (" << myrow << ", " << mycol << ")" << endl;
    if(C == nullptr) std::cout << "Memory allocation error for C on proc (" << myrow << ", " << mycol << ")" << endl;


    /*
    for(size_t k = 0; k < mlocal * nlocal; k++) {
        A[k] = myrank + 1.0 + k;
        B[k] = myrank*myrank +1.0 +k;
        B[k] = myrank +1.0 +k;
        C[k] = 0.0;
    }
    */

    for(size_t it = 0; it < mlocalrow * klocalcol; it++) {
        A[it] = myrank + 1.0 + it;
    }

    for(size_t it = 0; it < klocalrow * nlocalcol; it++) {
        B[it] = myrank + 1.0 + it;
    }

    for(size_t it = 0; it < mlocalrow * nlocalcol; it++) {
        C[it] = 0.0;
    }

    int descA[9];
    int descB[9];
    int descC[9];

    int info;
    int ldA = mlocalrow>1?mlocalrow:1;
    int ldB = klocalrow>1?klocalrow:1;
    int ldC = mlocalrow>1?mlocalrow:1;

    descinit_(descA, &m, &k, &mb, &kb, &firstprocessrowcol, &firstprocessrowcol, &ctx, &ldA, &info);
    assert(info == 0);
    descinit_(descB, &k, &n, &kb, &nb, &firstprocessrowcol, &firstprocessrowcol, &ctx, &ldB, &info);
    assert(info == 0);
    descinit_(descC, &m, &n, &mb, &nb, &firstprocessrowcol, &firstprocessrowcol, &ctx, &ldC, &info);
    assert(info == 0);



    //calling pdgemm function here
    double alpha = 1.0;
    double beta = 0;


    
    int one=1;
    char notrans='N';



    pdgemm_(&notrans, &notrans, &m, &n, &k, &alpha, A, &one, &one, descA,
            B, &one, &one, descB, &beta, C, &one, &one, descC);

    //pdgemm_('N', 'N', n, n, n, 1, A, 1, 1, descA, B, 1, 1, descB, 0, C, 1, 1, descC);
    //printf("My rank = %d, output values = %lf, %lf, %lf, %lf, %lf\n", myrank, C[0], C[1], C[2], C[3], C[4]);
    for(size_t it=0; it < mlocalrow * nlocalcol; it++)
    {
        printf("My rank = %d, C[%d] = %lf\n", myrank, it, C[it]);
    }
    
    delete [] A;
    delete [] B;
    delete [] C;
    blacs_gridexit_(&ctx);
    MPI_Finalize();
    return 0;
}
