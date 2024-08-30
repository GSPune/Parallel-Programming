#include<mpi.h>
#include<iostream>
#include<vector>
using namespace std;
#include <cmath>
//Distributed System
int main(int argc,char** argv)
{ 
	int rank,p;
	int n;
    n = 4;
	int A[n][n] = {{1,5,0,4},{2,6,5,0},{0,7,3,2},{4,0,1,1}};
	int x[n] = {1,2,4,6};
    
    MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int rem = n % p;
    int sendCounts[p], displs[p];
    displs[0] = 0;
    for(int i = 0; i < p; ++i) sendCounts[i] = ((n/p)+((i < rem)?1:0))*n;
    for(int j = 1; j < p; ++j) displs[j] = sendCounts[j-1] + displs[j-1];
    int recvbuf[sendCounts[rank]];
    //Distributing rows of A

    MPI_Scatterv(&A,sendCounts,displs,MPI_INT,&recvbuf,sendCounts[rank],MPI_INT,0,MPI_COMM_WORLD);
    //send sCnt elements of x to each process
    int disp[p];disp[0] = 0;
    for(int j = 1; j < p; ++j) disp[j] = sendCounts[j-1]/n + displs[j-1];
	// int rcv[sendCounts[rank]];
    // MPI_Scatterv(&x,sendCounts,disp,MPI_INT,&rcv,sendCounts[rank],MPI_INT,0,MPI_COMM_WORLD);
	// // MPI_Scatter(&x,sCnt,MPI_INT,rcv,sCnt,MPI_INT,0,MPI_COMM_WORLD);
    // //Gather all xi's on each process
    // int xGathered[n];
	// MPI_Allgatherv(rcv,sendCounts[rank],MPI_INT,xGathered,sendCounts,disp,MPI_INT,MPI_COMM_WORLD);

    // if(rank == 0){
    //     //for(int i = 0; i < p; ++i) cout << sendCounts[i] << ", ";
    //     for(int i = 0; i < p; ++i) cout << displs[i] << ", ";
    // }
    cout << "All Elements received on proc " << rank << endl;
    for(int i = 0; i < sendCounts[rank]; ++i) cout << recvbuf[i] << ", ";
    cout << endl;

	MPI_Finalize();
	return 0;
}