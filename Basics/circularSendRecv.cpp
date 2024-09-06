#include<mpi.h>
#include<bits/stdc++.h>
using namespace std;

int main(int argc,char** argv){
    int rank,size,i=0;
    // vector<int> A = {1,2,3,4,5}, B = {7,8,9,10,11};
    int A[5] = {1,2,3,4,5}, B[5] = {7,8,9,10,11};

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank % 2 == 0){
        MPI_Send(&A,5,MPI_INT,(rank+1)%size,1,MPI_COMM_WORLD);
        MPI_Recv(&B,5,MPI_INT,(rank-1+size)%size,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    else{
        MPI_Recv(&B,5,MPI_INT,(rank-1+size)%size,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(&A,5,MPI_INT,(rank+1)%size,1,MPI_COMM_WORLD);
    }
    cout << "Process Rank " << rank << endl;
    for(auto e:B) cout << e << " ";
    cout << endl;
    MPI_Finalize();
	return 0;
}