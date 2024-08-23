#include<mpi.h>
#include<iostream>
#include<vector>
using namespace std;
#include <cmath>
//Distributed System
int main(int argc,char** argv)
{ 
	int rank,size;
	int rows,y;
	int A[4][4] = {{1,5,0,4},{2,6,5,0},{0,7,3,2},{4,0,1,1}};
	//considering square matrices ...row major?
	//vector<vector<int>> A{{1,5,0,4,5},{2,6,5,0},{0,7,3,2},{4,0,1,1}};
	rows = 4;//A.size();
	int x[4] = {1,2,4,6};
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//cout << "\n" ;
	//if(rank == 0) cout << A.size() << endl;
	int sendCountRows = rows/size;
	int *recvbuf = new int[sendCountRows*rows];
	//send (rows/size) rows OF A to each process 
	MPI_Scatter(&A,sendCountRows*rows,MPI_INT,recvbuf,sendCountRows*rows,MPI_INT,0,MPI_COMM_WORLD);
	
	int sCnt = (rows/size);
	int *rcv = new int[sCnt];
	//send sCnt elements of x to each process
	MPI_Scatter(&x,sCnt,MPI_INT,rcv,sCnt,MPI_INT,0,MPI_COMM_WORLD); 
	//gather all the component of X on each processor
	//int gCnt = (size)*sCnt;
	int *xGathered = new int[rows];
	MPI_Allgather(rcv,sCnt,MPI_INT,xGathered,sCnt,MPI_INT,MPI_COMM_WORLD);
	for (int i = 0; i < size; ++i){
		if(rank == i){
			y = 0;
			//cout << "All Elements received on proc " << rank << endl;
			//for(int i = 0; i < sendCountRows*rows; ++i) cout << recvbuf[i] << ", ";
			//cout << "and..";
			//for(int i = 0; i < sCnt; ++i) cout << rcv[i] << ", ";
			//cout << endl;
			//cout << "All Elements gathered on proc " << rank << endl;
			//for(int i = 0; i < rows; ++i) cout << xGathered[i] << ", ";
			for(int i = 0; i < rows ; ++i){ 
				y += recvbuf[i]*xGathered[i];
				
			}
			cout << y << endl;
		}
	}
	int *bGathered = new int[size];
	//MPI_Gather collects data from all processes in a given communicator and concatenates them in the given buffer on the specified process. 
	//The concatenation order follows that of the ranks.
	MPI_Gather(&y,1,MPI_INT,bGathered,1,MPI_INT,0,MPI_COMM_WORLD);
	if(rank == 0){
		for(int i = 0; i < rows ; ++i){
			//cout << "y" << i << " is " << bGathered[i] << endl;
		}
	}
	MPI_Finalize();
	return 0;
}