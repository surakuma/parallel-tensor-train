#include <iostream>
#include <cmath>
#include <cassert>
#include <cstring>
#include "mpi.h"

#include <random>
#include <vector>
#include <algorithm>
#include <vector>


using namespace std;

extern "C" {


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



extern int dgemm_(const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K,
                  const double* ALPHA, const double* A, const int* LDA, const double* B,
          const int* LDB, const double* BETA, double* C, const int* LDC);
}


template <class T>
void rearrange_data(T* data, int nrows_per_block, int ncols_per_block, int ncols)
{
    T* tmp = new T [nrows_per_block * ncols_per_block * ncols];
    memcpy(tmp, data, nrows_per_block * ncols_per_block * ncols * sizeof(T));

    for(int icol_per_block =0; icol_per_block <ncols_per_block; icol_per_block++)
    {
        for(int icol=0; icol<ncols; icol++)
        {
            //TODO: check this computation
            int prev_offset = icol * nrows_per_block * ncols_per_block + icol_per_block * nrows_per_block ;
            int curr_offset = icol_per_block * nrows_per_block * ncols + icol * nrows_per_block; 
            memcpy(&data[curr_offset], &tmp[prev_offset], ncols_per_block*sizeof(T));
        }

    }

    delete [] tmp;
}
//Perform A^T(1:n1, 1:m1)*B(1:m2, 1:n2)
double* multiplyPortionofATwithB(double *A, const int& m1, const int& n1, const int& lda, double *B, const int& m2, const int& n2)
{
    assert(m1 == m2);
    double *output = new double[n1*n2];
    
    int ldb = m2;
    int ldc = n1;

    char transA = 'T';
    char transB = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    dgemm_(&transA, &transB, &n1, &n2, &m1, &alpha, A, &lda, B, &ldb, &beta, output, &ldc);
    //dgemm_(&trans, &trans, &number_of_my_rows, &required_rank, &number_of_my_columns, &alpha, our_data, &lda, rand_matrix, &ldb, &beta, output_after_rand_matrix, &ldc);
    return output;
}


double* reduceAlongColsAndGather(MPI_Comm original_comm, int nproc_row, int nproc_col, double *local_data, int nrows, int ncols)
{
    int rank, nprocess;
    MPI_Comm_rank(original_comm, &rank);
    MPI_Comm_size(original_comm, &nprocess);

    assert(nprocess == nproc_row * nproc_col);


    int rank_row = rank % nproc_row;
    int rank_col = rank / nproc_row;

    MPI_Comm col_comm;
    MPI_Comm_split(original_comm, rank_col, rank, &col_comm);

    double *col_reduced_data = new double [nrows*ncols];
    MPI_Reduce(local_data, col_reduced_data, nrows*ncols, MPI_DOUBLE, MPI_SUM, 0, col_comm);

    MPI_Barrier(col_comm);
    MPI_Comm_free(&col_comm);
    
    double *output = nullptr;
    if(rank == 0)
    {
        output = new double [nrows*ncols*nproc_col];
        std::vector<int> recv_counts;
        std::vector<int> displacements;
        //int recv_counts = new int [nprocess];

        int offset = 0;
        for(int iprocess =0; iprocess <nprocess; iprocess++)
        {
            int val =0;
            if(iprocess % nproc_row == 0)
                val = nrows*ncols;
            recv_counts.push_back(val);
            displacements.push_back(offset);
            offset += val;
        }

        MPI_Gatherv(col_reduced_data, nrows*ncols, MPI_DOUBLE, output, recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, original_comm);

    }
    else
    {
        MPI_Gatherv(col_reduced_data, nrows*ncols, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0, original_comm);
    }


    delete [] col_reduced_data;

    return output;

}

void computesvd(double* &U, double* &S, double* &VT, double* input, int m, int n)
{
    int minmn = min(m,n);

	char schar='S';
    U = new double [m*minmn];
    S = new double [minmn];
    VT = new double [minmn * n];

    int lda = m;
    int ldu = m;
    int ldvt = minmn;
    int lwork = 10*minmn + max(m,n);
    int info;
    vector<double> work(lwork);

    assert(input != nullptr);

    dgesvd_(schar, schar, m, n, input, lda, S, U, ldu, VT, ldvt, work.data(), lwork, info);
}


double* multiply(double* input1, const int &m1, const int&n1, double* input2, const int &m2, const int &n2)
{
    assert(n1 == m2);

    double *output = new double[m1*n2];

    int lda = m1;
    int ldb = m2;
    int ldc = m1;
    char trans = 'N';

    double alpha = 1.0;
    double beta = 0.0;

    dgemm_(&trans, &trans, &m1, &n2, &n1, &alpha, input1, &lda, input2, &ldb, &beta, output, &ldc);
    return output;
}


//inplace computation when m >= n
double* computeQ(double *A, int m, int n)
{
    double *copyA = new double [m*n];
    int lda = m;
    vector<double> tau(min(m,n));
    int lwork = n*n;
    vector<double> work(lwork);
    int info;

    dgeqrf_(m, n, copyA, m, tau.data(), work.data(), lwork, info);
    
    dorgqr_(m, n, n, copyA, m, tau.data(), work.data(), lwork, info);

    return copyA;
}

int main(int argc, char* argv[])
{

    int rank, nprocess;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int nproc_row = sqrt(nprocess);
    int nproc_col = nprocess/nproc_row;

    int rank_row = rank % nproc_row;
    int rank_col = rank / nproc_row;


    //assert(nprocess == 16);
    assert(nproc_row == nproc_col);

    int number_of_my_columns = 4;
    int number_of_my_rows = 4;
    int number_of_elements = number_of_my_rows * number_of_my_columns;

    double *our_data = new double [number_of_elements];

    if(rank == 0)
    {

        int global_nrow = nproc_row * number_of_my_rows;
        int global_ncol = nproc_col *number_of_my_columns;

        double input[global_nrow * global_ncol];

        for(int i=0; i<global_ncol; i++)
            for(int j=0; j<global_nrow; j++)
                input[i*global_nrow + j] = i*global_nrow +j+1;

        for(int j=0; j<nproc_col; j++)
            for(int i=0; i<nproc_row; i++)
            {
                int recv_rank = j*nproc_row + i;


                int start_row = i * global_nrow/nproc_row;
                int end_row = (i+1) * global_nrow/nproc_row - 1;

                int start_col = j * global_ncol/nproc_col;
                int end_col = (j+1) * global_ncol/nproc_col -1;


                double temp[number_of_elements];

                int nelements = 0;
                for(int it=start_col; it<=end_col; it++)
                    for(int jt=start_row; jt<=end_row; jt++)
                        temp[nelements++] = input[it*global_nrow + jt];


                assert(nelements == number_of_elements);

//                    if(recv_rank ==2)
//                    {
//                        for(int ielement=0; ielement<nelements; ielement++)
//                            std::cout << temp[ielement];
//                    
//                        std::cout << "\n";
//                    }

                if(recv_rank == rank) 
                {
                    memcpy(our_data, temp, number_of_elements*sizeof(double));
                    continue;
                }
                
                int tag = recv_rank;
                MPI_Send(temp, number_of_elements, MPI_DOUBLE, recv_rank, tag, MPI_COMM_WORLD);
            }
    }
    else
    {
        MPI_Status status;
        int tag = rank;
        MPI_Recv(our_data, number_of_elements, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

//        if(rank ==2)
//        {
//            for(int ielement=0; ielement<nelements; ielement++)
//                std::cout << temp[ielement];
//
//            std::cout << "\n";
//        }
    }


    MPI_Barrier(MPI_COMM_WORLD);

    //end of data distribution

    int required_rank = 4;
    int number_of_elements_in_rmatrix = number_of_my_columns * required_rank;

    double *rand_matrix = new double [number_of_elements_in_rmatrix];
    double *output_after_rand_matrix = new double [number_of_my_rows * required_rank];




    int seed = rank_col;

    //std::random_device rd;
    //std::default_random_engine generator (seed);
    //std::mt19937 e2(rd());

    // construct a trivial (32-bit) random generator engine from column rank
    std::mt19937 generator (seed);
    std::normal_distribution<double> dist(0.0, 1.0);

    for(int i=0; i<number_of_elements_in_rmatrix; i++)
        rand_matrix[i] = dist(generator);


    //calling GEMM kernel for each processor

    char trans = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    int lda = number_of_my_rows;
    int ldb = number_of_my_columns;
    int ldc = number_of_my_rows;
    
    dgemm_(&trans, &trans, &number_of_my_rows, &required_rank, &number_of_my_columns, &alpha, our_data, &lda, rand_matrix, &ldb, &beta, output_after_rand_matrix, &ldc);


    //perform MPI reduce operation
    
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank_row, rank, &row_comm);


    double *row_reduced_data = new double [number_of_my_rows * required_rank];
    MPI_Reduce(output_after_rand_matrix, row_reduced_data, number_of_my_rows * required_rank, MPI_DOUBLE, MPI_SUM, 0, row_comm);

    MPI_Barrier(row_comm);
    MPI_Comm_free(&row_comm);
    //perform gather operation and rearrange the data (also look at MPI_Gatherv)

    double *qfactor;
    if(rank == 0)
    {
        double *gathered_data = new double [number_of_my_rows * required_rank * nproc_row];
        std::vector<int> recv_counts;
        std::vector<int> displacements;
        //int recv_counts = new int [nprocess];

        int offset = 0;
        for(int iprocess =0; iprocess <nprocess; iprocess++)
        {
            int val =0;
            if(iprocess / nproc_row == 0)
                val =number_of_my_rows * required_rank;
            recv_counts.push_back(val);
            displacements.push_back(offset);
            offset += val;
        }

        MPI_Gatherv(row_reduced_data, number_of_my_rows * required_rank , MPI_DOUBLE, gathered_data, recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //rearrange the combined data in row major order
        rearrange_data(gathered_data, number_of_my_rows, required_rank, nproc_row);
        //Apply QR (Q would be in input variable)
        qfactor = computeQ(gathered_data, number_of_my_rows*nproc_row, required_rank);

        delete [] gathered_data;
    }
    else
    {
        MPI_Gatherv(row_reduced_data, number_of_my_rows * required_rank , MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    delete [] row_reduced_data;
    delete [] output_after_rand_matrix;
    delete [] rand_matrix;
    //std::cout << "rank = " << rank << " (rank_row, rank_col) = (" << rank_row << ", " << rank_col << ")" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank != 0)
        qfactor = new double [number_of_my_rows * required_rank * nproc_row];

    MPI_Bcast(qfactor, number_of_my_rows * required_rank * nproc_row, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //perform Q^T A
    int qta_lda = number_of_my_rows * nproc_row;
    double *qta_local = multiplyPortionofATwithB(&qfactor[number_of_my_rows * rank_row], number_of_my_rows, required_rank, qta_lda, our_data, number_of_my_rows, number_of_my_columns); 


    double *qta_combined_on_root = reduceAlongColsAndGather(MPI_COMM_WORLD, nproc_row, nproc_col, qta_local, required_rank, number_of_my_columns);


    if(rank == 0)
    {
        //process qfactor and qta_combined_on_root
        //svd of qta_combined_on_root

        int nrows = required_rank;
        int ncols = nproc_col *number_of_my_columns;
        double *U=nullptr;
        double *S=nullptr;
        double *VT=nullptr;

        computesvd(U, S, VT, qta_combined_on_root, nrows, ncols);

        double *combinedU = multiply(qfactor, number_of_my_rows*nproc_row, required_rank, U, nrows, nrows);

        delete [] U;

        U = combinedU;

        delete [] qta_combined_on_root;

        std::cout << "Singular values are " << std::endl;
        for(int i=0; i<required_rank; i++)
            std::cout << S[i] << " ";
        
        std::cout << endl;
    }
    else
        assert(qta_combined_on_root == nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
/*

    //delete qmatrix from all processors
    delete [] qfactor;
    //delete local data
    delete [] our_data;

    //perform svd of Q^T A
    //Multiply results with Q


//    //multily every block with a random matrix of size
//
//    MPI_Comm row_comm;
//    MPI_Comm_split(MPI_COMM_WORLD, rank_row, rank, &row_comm);
//
//    double *row_reduced_data = new double [4*4];
//    MPI_Reduce(our_data, row_reduced_data, 16, MPI_DOUBLE, MPI_SUM, 0, row_comm);
//
//
//    if(rank<4)
//    {
//        std::cout << "rank = " << rank << "\n";
//        for(int irow=0; irow<4; irow++)
//        {
//            for(int icol=0; icol<4; icol++)
//                //std::cout << our_data[icol*4 + irow] << " ";
//                std::cout << row_reduced_data[icol*4 + irow] << " ";
//
//            std::cout << "\n";
//        }
//    }
//    
//    MPI_Barrier(row_comm);
//    delete [] row_reduced_data;
//
//    MPI_Comm_free(&row_comm);

  */  
    MPI_Finalize(); 
    return 0;
}
