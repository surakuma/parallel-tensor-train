#include <iostream>
#include <cmath>
#include <cassert>
#include <cstring>
#include "mpi.h"

#include <random>
#include <vector>



using namespace std;



extern "C" {
extern int dgemm_(const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K,
                  const double* ALPHA, const double* A, const int* LDA, const double* B,
          const int* LDB, const double* BETA, double* C, const int* LDC);
}

rearrange_data(gathered_data, number_of_my_rows, required_rank, nproc_col);

template <class T>
void rearrange_data(T* data, int nrows_per_block, int ncols_per_block, int ncols)
{
    T* tmp = new T [nrows_per_block * ncols_per_block * ncols];
    memcpy(tmp, data, nrows_per_block * ncols_per_block * ncols 8 sizeof(T));

    for(int icol_per_block =0; icol_per_block <ncols_per_block; icol_per_block++)
    {
        for(int icol=0; icol<ncols; icol++)
        {
            //TODO: check this computation
            int prev_offset = icol * nrows_per_block * ncols_per_block + icol_per_block * nrows_per_block ;
            int curr_offset = icol_per_block * nrows_per_block * ncols + icol * nrows_per_block; 
            memcpy(&gathered_data[curr_offset], &tmp[prev_offset], ncols_per_block*sizeof(T));
        }

    }

    delete [] tmp;
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

    assert(nprocess == 16);


    double *our_data = new double [4*4];
    if(rank == 0)
    {
        double input[16*16];
        for(int i=0; i<16; i++)
            for(int j=0; j<16; j++)
                input[i*16 + j] = i*16 +j+1;

        for(int j=0; j<nproc_col; j++)
            for(int i=0; i<nproc_row; i++)
            {
                int recv_rank = j*nproc_row + i;


                int start_row = i * 16/nproc_row;
                int end_row = (i+1) * 16/nproc_row - 1;

                int start_col = j * 16/nproc_col;
                int end_col = (j+1) * 16/nproc_col;


                double temp[4*4];

                    int nelements = 0;
                    for(int it=start_col; it<=end_col; it++)
                        for(int jt=start_row; jt<=end_row; jt++)
                            temp[nelements++] = input[it*16+jt];

//                    if(recv_rank ==2)
//                    {
//                        for(int ielement=0; ielement<nelements; ielement++)
//                            std::cout << temp[ielement];
//                    
//                        std::cout << "\n";
//                    }

                if(recv_rank == rank) 
                {
                    memcpy(our_data, temp, 4*4*sizeof(double));
                    continue;
                }
                
                int tag = recv_rank;
                MPI_Send(temp, 4*4, MPI_DOUBLE, recv_rank, tag, MPI_COMM_WORLD);
            }
    }
    else
    {
        MPI_Status status;
        int tag = rank;
        MPI_Recv(our_data, 4*4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

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

    int number_of_my_columns = 4;
    int number_of_my_rows = 4;
    int number_of_elements = 4*4;
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

    if(rank == 0)
    {
        double *gathered_data = new double [number_of_my_rows * required_rank * nproc_col];
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

        MPI_Gatherv(row_reduced_data, number_of_my_rows * required_rank , MPI_DOUBLE, gathered_data, recv_counts.data(), displacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //rearrange the combined data in row major order
        rearrange_data(gathered_data, number_of_my_rows, required_rank, nproc_col);
        //Apply QR

        delete [] gathered_data;
    }
    else
    {
        MPI_Gatherv(row_reduced_data, number_of_my_rows * required_rank , MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    //perform svd of Q^T A
    //Multiply results with Q


    delete [] row_reduced_data;
    delete [] output_after_rand_matrix;
    delete [] rand_matrix;


    //multily every block with a random matrix of size

    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank_row, rank, &row_comm);

    double *row_reduced_data = new double [4*4];
    MPI_Reduce(our_data, row_reduced_data, 16, MPI_DOUBLE, MPI_SUM, 0, row_comm);


    if(rank<4)
    {
        std::cout << "rank = " << rank << "\n";
        for(int irow=0; irow<4; irow++)
        {
            for(int icol=0; icol<4; icol++)
                //std::cout << our_data[icol*4 + irow] << " ";
                std::cout << row_reduced_data[icol*4 + irow] << " ";

            std::cout << "\n";
        }
    }
    
    MPI_Barrier(row_comm);
    delete [] row_reduced_data;

    MPI_Comm_free(&row_comm);

    delete [] our_data;
    
    MPI_Finalize(); 
    return 0;
}
