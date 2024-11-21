#include<mpi.h>
#include<iostream>
using namespace std;

int main(int argc,char** argv)
{
	// cout << "Outside MPI Env 1 and rank: unknown?" << endl;
	int rank,size;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm Cart_topo;
    int d = 2, dims[2] = {2,3}, period[2] = {1,1}, reorder = 1; 
    //dim 2*3 dim0 & dim 1
    //wraparound both axis {1,1} -> {true,true}
    //reorder -> renumber the ranks?
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,period,reorder,&Cart_topo);

    int new_rank;
    //Gets the rank of the calling MPI process in the communicator specified.
    MPI_Comm_rank(Cart_topo, &new_rank);

    // Get my coordinates in the new communicator
    int my_coords[2];
    MPI_Cart_coords(Cart_topo, new_rank, 2, my_coords); // Rank to coordinates

	// cout << "Hello World" << endl;
	cout << "I am being printed on processor " << rank << ", with coordinates: (" << my_coords[0] <<", " << my_coords[1] << ") out of " << size << " total processors" << endl;

    int left,right,shift = 1,change = 1;

    //MPI_Cart_shift( MPI_Comm comm , int direction , int disp/stepSize , int* rank_source , int* rank_dest);
    MPI_Cart_shift(Cart_topo,1 ,1,&left ,&right);
    cout << "\nLeft :: " << left << " and right :: " << right << endl;
	MPI_Finalize();
	// cout << "Outside MPI Env 2 and rank: unknown?" << endl;
	return 0;
}