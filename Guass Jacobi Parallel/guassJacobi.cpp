#include<mpi.h>
#include <bits/stdc++.h>
using namespace std;
#include <cmath>
//Distributed System
int n = 4;
int main(int argc,char** argv)
{ 
	int rank,p;
	int A[n][n],b[n];
	double x[n];

    MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (rank == 0){
		ifstream fin;
		fin.open("./input.txt");
		//fin >> n;
		// int A[n][n] = {{1,5,0,4},{2,6,5,0},{0,7,3,2},{4,0,1,1}};
		
		for(int i = 0 ; i < n ; ++i){
			for (int j = 0 ; j < n; ++j){
			    fin >> A[i][j];
			}
		}
		//int b[n];
		for (int k = 0; k < n; ++k){
			fin >> b[k];x[k] = 0;
		}
		fin.close();
		
		for(int i = 0 ; i < n ; ++i){
			for (int j = 0 ; j < n; ++j){
			    cout << A[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Program Results passed onto to output file!" << endl;
		
	}
	
    int rem = n % p, k = n/p;
    int sendCountRows = k*n;
	int recvbuf[sendCountRows];
	//send scR rows OF A to each process 
	//Distributing rows of A
	MPI_Scatter(&A,sendCountRows,MPI_INT,recvbuf,sendCountRows,MPI_INT,0,MPI_COMM_WORLD);
	
	//Now to distribute elements of b
	int sendC = k;
	int recv[sendC];
	MPI_Scatter(&b,sendC,MPI_INT,recv,sendC,MPI_INT,0,MPI_COMM_WORLD);
	
	//cout << "Process " << rank << endl;
 	//or(int i = 0; i < sendC; ++i){
 	//	cout << recv[i] << " ";
 	//}
 	//cout << endl;
 	MPI_Bcast(&x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
 	for (int it = 0; it < 2; it++){
 		double xn[k];
 		for(int i = 0; i < k ; ++i){
 			cout << recv[i] << endl;
 			xn[i] = recv[i];
 			for (int j = 0; j < n; ++j){
 				if((i*n+rank*k+i) != j){
 					xn[i] -= (recvbuf[j+n*i] * x[j]);
 				}
 			}
 			xn[i] /= (double)recvbuf[i*n+rank*k+i];
 			//cout << x[i] << " " << endl;
 			//cout << xNGathered[i] << " ";
 		}
 		
 		//Gather the Xnew s from the other processors
 		double xNGathered[n];
		MPI_Allgather(&xn,k,MPI_DOUBLE,xNGathered,k,MPI_DOUBLE,MPI_COMM_WORLD);
		
		
	 	
 		for (int p = 0; p < n; p++){
 			x[p] = xNGathered[p];
 		}
 		
 		if(rank == 1){
		for(int i = 0; i < n; ++i){
	 		//cout << x[i] << endl;
	 	}}
 	}
 	
	if(rank == 0){
	 	for(int i = 0; i < n; ++i){
	 		cout << x[i] << " ";
	 	}
	 	cout << endl;
	}
	MPI_Finalize();
	return 0;
}
