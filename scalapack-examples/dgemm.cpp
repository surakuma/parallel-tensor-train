#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <algorithm>
#include <cassert>
#include "mpi.h"

/**
 * Compile with something like
 *
 * [*] With Intel MKL
 *
 * mpiicpc dgemm.cpp -o dgemm  -I${MKLROOT}/include ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
 *
 * [*] In general
 *
 * mpicxx dgemm.cpp -o dgemm \
 *     -L/.../scalapack/2.1.0/lib \
 *     -lscalapack
 */

extern "C" void blacs_get_(int*, int*, int*);
extern "C" void blacs_pinfo_(int*, int*);
extern "C" void blacs_gridinit_(int*, char*, int*, int*);
extern "C" void blacs_gridinfo_(int*, int*, int*, int*, int*);
extern "C" void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern "C" void pdgemm_(char*, char*, int*, int*, int*, double*, double*, int*, int*, int*, double*, int*, int*, int*, double*, double*, int*, int*, int*);
extern "C" void blacs_gridexit_(int*);
extern "C" int numroc_(int*, int*, int*, int*, int*);

int main(int argc, char **argv) {
    int izero=0;
    int ione=1;
    int myrank_mpi, nprocs_mpi;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

    int n = 1000;       // (Global) Matrix size
    int nprow = 2;   // Number of row procs
    int npcol = 2;   // Number of column procs
    int nb = 256;      // (Global) Block size
    char notrans = 'N';
    char layout='R'; // Block cyclic, Row major processor mapping

    printf("Usage: ./test matrix_size block_size nprocs_row nprocs_col\n");

    if(argc > 1) {
        n = atoi(argv[1]);
    }
    if(argc > 2) {
        nb = atoi(argv[2]);
    }
    if(argc > 3) {
        nprow = atoi(argv[3]);
    }
    if(argc > 4) {
        npcol = atoi(argv[4]);
    }

    assert(nprow * npcol == nprocs_mpi);

    // Initialize BLACS
    int iam, nprocs;
    int zero = 0;
    int ictxt, myrow, mycol;
    blacs_pinfo_(&iam, &nprocs) ; // BLACS rank and world size
    blacs_get_(&zero, &zero, &ictxt ); // -> Create context
    blacs_gridinit_(&ictxt, &layout, &nprow, &npcol ); // Context -> Initialize the grid
    blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    // Compute the size of the local matrices
    int mpA    = numroc_( &n, &nb, &myrow, &izero, &nprow ); // My proc -> row of local A
    int nqA    = numroc_( &n, &nb, &mycol, &izero, &npcol ); // My proc -> col of local A

    printf("Hi. Proc %d/%d for MPI, proc %d/%d for BLACS in position (%d,%d)/(%d,%d) with local matrix %dx%d, global matrix %d, block size %d\n",myrank_mpi,nprocs_mpi,iam,nprocs,myrow,mycol,nprow,npcol,mpA,nqA,n,nb);

    // Allocate and fill the matrices A, B and C
    // A[I,J] = (I == J ? 5*n : I+J)
    double *A;
    double *B;
    double *C;
    A = (double *)calloc(mpA*nqA,sizeof(double)) ;
    B = (double *)calloc(mpA*nqA,sizeof(double)) ;
    C = (double *)calloc(mpA*nqA,sizeof(double)) ;
    if (A==NULL){ printf("Error of memory allocation A on proc %dx%d\n",myrow,mycol); exit(0); }
    if (B==NULL){ printf("Error of memory allocation B on proc %dx%d\n",myrow,mycol); exit(0); }
    if (C==NULL){ printf("Error of memory allocation C on proc %dx%d\n",myrow,mycol); exit(0); }
    for(size_t k = 0; k < mpA * nqA; k++) {
	A[k] = 1.0;
        B[k] = 1.0;
        C[k] = 0.0;
    }

    // Create descriptor
    int descA[9];
    int descB[9];
    int descC[9];
    int info;
    int lddA = mpA > 1 ? mpA : 1;
    int lddB = lddA;
    int lddC = lddA;
    descinit_( descA,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &lddA, &info);
    if(info != 0) {
        printf("Error in descinit, info = %d\n", info);
    }
    descinit_( descB,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &lddB, &info);
    if(info != 0) {
        printf("Error in descinit, info = %d\n", info);
    }
    descinit_( descC,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &lddC, &info);
    if(info != 0) {
        printf("Error in descinit, info = %d\n", info);
    }

    // Run dgemm and time
    double alpha = 1.0;
    double beta = 0.0;
    double MPIt1 = MPI_Wtime();
    printf("[%dx%d] Starting dgemm\n", myrow, mycol);
    pdgemm_(&notrans, &notrans, &n, &n, &n,
            &alpha,
            A, &ione, &ione, descA,
	    B, &ione, &ione, descB,
            &beta,
            C, &ione, &ione, descC);
    double MPIt2 = MPI_Wtime();
    printf("[%dx%d] Done, time %e s.\n", myrow, mycol, MPIt2 - MPIt1);

    size_t error = 0;
    for(size_t k = 0; k < mpA * nqA; k++) error += (C[k] != n);
    if(error) printf("\n\n====> ERROR (%zd)\n\n",error);

    // Exit and finalize
    free(A); free(B); free(C);
    blacs_gridexit_(&ictxt);
    MPI_Finalize();
    return 0;
}
