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


            void pdgesvd_(const char& jobu, const char& jobvt, const int& matrix_height,
            const int& matrix_width,
            double* matrix, const int& ia, const int& ja, int* matrix_descriptor,
            double* singular_values,
            double* u_matrix, const int& iu, const int& ju, int* u_descriptor,
            double* vt_matrix, const int& ivt, const int& jvt, int* vt_descriptor,
            double* work, int const& lwork, int const& info);

    */

    //declaration from http://www.netlib.org/scalapack/explore-html/d2/de6/pdgesvd_8f_source.html

    void pdgesvd_(
            char* JOBU,
            char* JOBVT,
            int* M,
            int* N,
            double* A,
            int* IA,
            int* JA,
            int* DESCA,
            double* S,
            double* U,
            int* IU,
            int* JU,
            int* DESCU,
            double* VT,
            int* IVT,
            int* JVT,
            int* DESCVT,
            double *work,
            int* lwork,
            int* info); 

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
    void blacs_gridinit_(int *ctx, char* order, int *nprow, int *npcol);
    
    /*
     * http://www.netlib.org/scalapack/explore-html/de/df8/_s_l__init_8f_source.html
     * Row major ordering of the processes
     */
    void sl_init_(int *ctx, int *nprow, int *npcol);

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


    int m = 4;
    int n = 6;

    int mb = (m+nprow - 1)/nprow;
    int nb = (n+npcol -1)/npcol;

    double globalA[4*6] = {
        0.170340,
        0.429255,
        0.528681,
        0.792433,
        0.479592,
        0.079933,
        0.556482,
        0.668215,
        0.872455,
        0.790935,
        0.623004,
        0.462168,
        0.689222,
        0.371529,
        0.153106,
        0.181449,
        0.783312,
        0.698733,
        0.984562,
        0.203782,
        0.856288,
        0.689554,
        0.844158,
        0.166957
    };
    int zero = 0;
    int ctx = -1;

    //retrieve a defualt system context
    blacs_get_(&zero, &zero, &ctx);


    //Row major mapping of processors
    sl_init_(&ctx, &nprow, &npcol);


    int myrow, mycol;

    //obtain identity of myrow and mycol in the grid
    blacs_gridinfo_(&ctx, &nprow, &npcol, &myrow, &mycol);

    std::cout << "my rank = " << myrank << ", gridid=(" << myrow << ", " << mycol << ")" << endl; ;

    int firstprocessrowcol = 0;
    //number of rows for the present process
    int mlocalrow    = numroc_(&m, &mb, &myrow, &firstprocessrowcol, &nprow); 
    int nlocalcol    = numroc_(&n, &nb, &mycol, &firstprocessrowcol, &npcol); 


    int startlocalrow = myrow * mlocalrow;
    if(mlocalrow != mb)
        startlocalrow +=  m % nprow;

    int endlocalrow = startlocalrow + mlocalrow -1;


    int startlocalcol = mycol * nlocalcol;
    if(nlocalcol != nb)
        startlocalcol += n % nlocalcol;

    int endlocalcol = startlocalcol + nlocalcol -1;

    std::cout << "Grid id=(" << myrow << ", " << mycol << "), localdatasize= (" << mlocalrow <<", " << nlocalcol << 
        "), row=(" << startlocalrow << ", " << endlocalrow << "), col=(" << startlocalcol << ", " << endlocalcol << ")" << endl;
    
    double *A = new double [mlocalrow * nlocalcol];
    if(A == nullptr) std::cout << "Memory allocation error for A on proc (" << myrow << ", " << mycol << ")" << endl;

    int nlocalelement = 0;
    for(int jt = startlocalcol; jt <=endlocalcol; jt++)
        for(int it = startlocalrow; it <=endlocalrow; it++)
            A[nlocalelement++] = globalA[jt*m + it];

    assert(nlocalelement == mlocalrow * nlocalcol);

    std::cout << "My rank = " << myrank << " elements = " << A[0] << " " << A[1] <<  " " <<A[2] << " " << A[3] << endl;

    delete [] A;
    blacs_gridexit_(&ctx);
    MPI_Finalize();
 
    return 0;
}





