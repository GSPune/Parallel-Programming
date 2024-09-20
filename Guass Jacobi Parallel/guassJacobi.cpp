#include<mpi.h>
#include<fstream>
#include <bits/stdc++.h>
using namespace std;
#include <cmath>
//Distributed System

int main(int argc,char** argv)
{ 	
	std::ifstream fin;
	fin.open("./input.txt");
	int n; fin >> n;
	// fin.close();
	int rank,p;
	int A[n][n],b[n];
	double x[n]; // decalred on each process but initialized on rank 0

    MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (rank == 0){
		// ifstream fin;
		// fin.open("./input.txt");
		
		for(int i = 0 ; i < n ; ++i){
			for (int j = 0 ; j < n; ++j){
			    fin >> A[i][j];
			}
		}
		//int b[n];
		for (int k = 0; k < n; ++k){
			fin >> b[k];x[k] = 0;
		}
		
		for(int i = 0 ; i < n ; ++i){
			for (int j = 0 ; j < n; ++j){
			    cout << A[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Program Results passed onto to output file!" << endl;
		
	}
	fin.close();
	
    int rem = n % p, k = n/p;
    int sendCountRows = k*n;
	int recvbuf[sendCountRows];
	//send {k*n} rows OF A to each process 
	//Distributing rows of A using MPI_Scatter 
	MPI_Scatter(&A,sendCountRows,MPI_INT,recvbuf,sendCountRows,MPI_INT,0,MPI_COMM_WORLD);
	
	//Now to distribute {k} elements of b to each process
	int sendC = k;
	int recv[sendC];
	MPI_Scatter(&b,sendC,MPI_INT,recv,sendC,MPI_INT,0,MPI_COMM_WORLD);
	
	//cout << "Process " << rank << endl;
 	//or(int i = 0; i < sendC; ++i){
 	//	cout << recv[i] << " ";
 	//}
 	//cout << endl;

	//Broadcasting X_old = {0,0,0,0} to all processes from source '0'
 	MPI_Bcast(&x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	// cout << "Test :: " << endl;
 	for (int it = 0; it < 50; it++){
 		double x_new[k];
 		for(int i = 0; i < k ; ++i){
 			// cout << recv[i] << endl;
 			x_new[i] = recv[i];
 			for (int j = 0; j < n; ++j){
 				if((i*n+rank*k+i)%n != j){
					// cout << recvbuf[j+n*i+rank*k] << endl;
 					x_new[i] -= (recvbuf[j+n*i] * x[j]);
 				}
 			}
 			x_new[i] /= (double)recvbuf[i*n+rank*k+i];
 			//cout << x[i] << " " << endl;
 			//cout << xNGathered[i] << " ";
 		}
 		
 		//Gather the X_new s from the other processors
 		double xNGathered[n];
		MPI_Allgather(&x_new,k,MPI_DOUBLE,xNGathered,k,MPI_DOUBLE,MPI_COMM_WORLD);
		
		//update the old x
 		for (int p = 0; p < n; p++){
 			x[p] = xNGathered[p];
 		}
 		
 		// if(rank == 1){
		// for(int i = 0; i < n; ++i){
	 	// 	//cout << x[i] << endl;
	 	// }}
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
