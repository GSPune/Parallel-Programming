#include <bits/stdc++.h>
#include<mpi.h>
using namespace std;

void printMatrix(double A[][5], int n){
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n + 1; ++j){
			cout << A[i][j] << "    ";
		}
		cout << endl;
	}
	cout << endl;
}

void Op(double *B,double *t,int size, int pivot){
	//GE requires a constant uniform scaling across the pivot row during row transformations!
	//The Gaussian elimination process requires proportional operations across the entire row to maintain the integrity of the system.
	double factor = B[pivot]/t[pivot];
	for(int i = 0; i < size; ++i){
		//B[i] -= (double)(B[i]*t[i])/t[pivot];
		 B[i] -= factor * t[i];
	}
}

int main(int argc,char** argv){
	int rank,p;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p); // p is the total number of processors used
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int n = 4;
	double A[4][5];
	double B[5],temp[5]; // for each processor
	double reducedA[4][5];
	if (rank == 0){
		double copyOfA[4][5] = {{2,3,-1,1,8},{3,2,1,-2,3},{1,1,-2,3,1},{4,0,-1,3,7}};
		memcpy(&A,&copyOfA,sizeof(A));
		//printMatrix(A,4);
	}
	MPI_Scatter(&A,5,MPI_DOUBLE,&B,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//A will receive the first row in B i.e R0
	if (rank == 0) memcpy(&temp,&B,sizeof(B));
	
	// Perform required number of operations
	for (int k = 1; k < 4; k ++){
		
		//Send necesary Row to all other processes
		MPI_Bcast(&temp, 5, MPI_DOUBLE, k-1, MPI_COMM_WORLD);
		if (rank > k-1){
			Op(B,temp,5,k-1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		//MPI_Barrier blocks all MPI processes in the given communicator until they all call this routine.
		//copy the pivot row for next iteration
		if (rank == k && rank != n-1) memcpy(&temp,&B,sizeof(B));
	}
	MPI_Gather(&B,5,MPI_DOUBLE,&reducedA[rank],5,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if (rank==0){
		printMatrix(reducedA,4);
		//Backsubstitution *
	}
	MPI_Finalize();
	return 0;
}

/*---------------------Debug stuff----------------------------
if (k == 2){
			cout << k << " Set of OPs::" << endl;
			cout << "Rank " << rank << "::  ";
			for (int t = 0; t < n+1; ++t){
				cout << B[t] << " ";
			}
			cout << endl;
		}
*/
