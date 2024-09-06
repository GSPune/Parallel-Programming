#include<mpi.h>
#include<iostream>
#include<vector>
using namespace std;
int total = 0;

int main(int argc,char** argv)
{
	int *rank = new int,size;
	// int total = 0;
	vector<int> A = {1,2,3,4,5,6,7,8,9,10,11,12};
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
	int s = 0;
	
	if(*rank == 0){
		for(int i = 1; i < 4; ++i){
			s += A[i-1];
		}
	}
		
	else if(*rank == 1){
		for(int i = 4; i < 7; ++i){
			s += A[i-1];
		}
	}
	else if(*rank == 2){
		for(int i = 7; i < 10; ++i){
			s += A[i-1];
		}
	}
	else if(*rank == 3){
		for(int i = 10; i < 13; ++i){
			s += A[i-1];
		}
	}
	MPI_Reduce(&s,&total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	if(*rank == 0){
		cout << "TOTAL Sum = " << total << endl;
	}
	MPI_Finalize();
	delete rank;
	return 0;
}