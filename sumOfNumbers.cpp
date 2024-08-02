#include<mpi.h>
#include<iostream>
#include<vector>
using namespace std;
int total = 0;

int main(int argc,char** argv)
{
	int *rank = new int,size,i=0;
	vector<int> A = {1,2,3,4,5,6,7,8,9,10,11,12};
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
	int s = 0;
	if(*rank == 0){
		for(int i = 1; i < 4; ++i){
			s += A[i-1];
		}
		//cout << "Rank0 Sum = " << s << endl;
		for(int j = 1; j < size; ++j)
		{ 
			int temp = 0;
			MPI_Recv(&temp,1,MPI_INT,j,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			total += temp;
		}
		total += s;
		cout << "TOTAL Sum = " << total << endl; 
	}
		
	else if(*rank == 1){
		for(int i = 4; i < 7; ++i){
			s += A[i-1];
		}
		//cout << "Rank1 Sum = " << s << endl;
		MPI_Send(&s,1,MPI_INT,0,1,MPI_COMM_WORLD);
	}
	else if(*rank == 2){
		for(int i = 7; i < 10; ++i){
			s += A[i-1];
		}
		//cout << "Rank2 Sum = " << s << endl;
		MPI_Send(&s,1,MPI_INT,0,1,MPI_COMM_WORLD);
	}
	else if(*rank == 3){
		for(int i = 10; i < 13; ++i){
			s += A[i-1];
		}
		//cout << "Rank3 Sum = " << s << endl;
		MPI_Send(&s,1,MPI_INT,0,1,MPI_COMM_WORLD);
	}
	MPI_Finalize();
	delete rank;
	return 0;
}