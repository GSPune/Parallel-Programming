#include<mpi.h>
#include<iostream>
#include<vector>
using namespace std;
int total = 0;
//Here size of array must be completely divisible by no. of processors
int main(int argc,char** argv)
{ 
	int *rank = new int,size,i=0;
	//vector<int> A = {1,2,3,4,5,6,7,8,9,10,11,12};
	int A[] = {1,2,3,4,5,6,7,8,9,10,11,12};
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
	int len = 12;
	int ofs = 0;
	int s = 0,offset = 12/size;
	int B[offset];
	
	if(*rank == 0){
		for(int i = 1; i < size; i++){
			ofs = i * offset;
			MPI_Send(&A[ofs],offset,MPI_INT,i,1,MPI_COMM_WORLD);
		}
		for(int i = 0; i < offset; ++i){
			s += A[i];
		}
		total += s;
		for(int j = 1; j < size; ++j)
		{ 
			int temp = 0;
			MPI_Recv(&temp,1,MPI_INT,j,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			total += temp;
		}
		cout << "TOTAL Sum = " << total << endl; 
	}
	else{
		MPI_Recv(&B,offset,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		for(int k = 0; k < offset; k++){
			s += B[k];
		}
		cout << "Rank Sum = " << s << endl;
		MPI_Send(&s,1,MPI_INT,0,1,MPI_COMM_WORLD);
	}
	MPI_Finalize();
	delete rank;
	return 0;
}