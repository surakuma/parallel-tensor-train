#include <iostream>
#include <cmath>
#include <cassert>
#include <cstring>
#include "mpi.h"

using namespace std;
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
