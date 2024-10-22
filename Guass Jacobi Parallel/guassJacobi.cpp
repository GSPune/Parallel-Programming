#include<mpi.h>
#include<fstream>
#include <bits/stdc++.h>
using namespace std;
#include <cmath>
//Distributed System

bool isDiagonallyDominant(int** mat,int rows,int cols){
    // |a_ii| >= Summation from j = 1 to n of |a_ij| s.t i != j
        for (int r = 0; r < rows; r++){
            double sum = 0.0;
            for (int c = 0; c < cols - 1; (c == r-1) ? c += 2 : c++){
                if (r == c && c == 0)
                    continue;
                sum += fabs(mat[r][c]);
            }
            if (fabs(mat[r][r]) <= sum)
                return false;
            sum = 0.0;
        }
		cout << "True" << endl;
        return true;
}

int main(int argc,char** argv)
{ 	
	std::ifstream fin;
	fin.open("./input.txt");
	int n; fin >> n;
	// fin.close();
	int rank,p;
	int b[n],N[n*n];
	int** A = new int* [n];
	for (int i = 0; i < n; ++i) {
        A[i] = new int[n];
    }
	double x[n]; // declared on each process but initialized on rank 0

    MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (n % p != 0){
			cout << "The program is meant to be run with n \% p == 0 " << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	if (rank == 0){
		// ifstream fin;
		// fin.open("./input.txt");
		
		for(int i = 0 ; i < n ; ++i){
			for (int j = 0 ; j < n; ++j){
			    fin >> A[i][j];
			}
		}

		if (!isDiagonallyDominant(A,n,n)){
			 cout << "Matrix is not Diagonally Dominant!" << endl;
			 MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		//int b[n];
		for (int k = 0; k < n; ++k){
			fin >> b[k];x[k] = 0;
		}
		
		// for(int i = 0 ; i < n ; ++i){
		// 	cout << (A+i) << " \n\n";
		// 	for (int j = 0 ; j < n; ++j){
		// 		cout <<i <<", " << j << ": " << (*(A+i)+j) << " ";
		// 	    // cout << A[i][j] << " ";
		// 	}
		// 	cout << endl << endl;
		// }

        //creating contiguous array..
        for (int i =0;i < n ; ++i){
            for (int j = 0 ; j < n; ++j){
                N[i*n+j] = A[i][j];
            }
        }
	
	}
	fin.close();
	
    int rem = n % p, k = n/p;
    int sendCountRows = k*n;
	// cout << "Count of rows is " << sendCountRows << endl;
	int recvbuf[sendCountRows];
	//send {k*n} rows OF A to each process 
	//Distributing rows of A using MPI_Scatter 
	MPI_Scatter(N,sendCountRows,MPI_INT,recvbuf,sendCountRows,MPI_INT,0,MPI_COMM_WORLD);
	
	cout << "\nProcess " << rank << " received: " << endl;
 	for(int i = 0; i < sendCountRows; ++i){
 		cout << recvbuf[i] << " ";
 	}
 	cout << endl;

	//Now to distribute {k} elements of b to each process
	int sendC = k;
	int recv[sendC];
	MPI_Scatter(&b,sendC,MPI_INT,recv,sendC,MPI_INT,0,MPI_COMM_WORLD);
	
	// cout << "\nProcess " << rank << "received: " << endl;
 	// for(int i = 0; i < sendC; ++i){
 	// 	cout << recv[i] << " ";
 	// }
 	// cout << endl;

	//Broadcasting X_old = {0,0,0,0} to all processes from source '0'
 	MPI_Bcast(&x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	// cout << "\nTest :: " << endl;
 	for (int it = 0; it < 100; it++){
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
 			// cout << x[i] << " " << endl;
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
		ofstream fout;
        fout.open("./output.txt");
		fout << "Solutions are :: " << endl; 
	 	for(int i = 0; i < n; ++i){
	 		fout << x[i] << " ";
	 	}
	 	fout << endl;
		fout.close();
		cout << "Program Results passed onto to output file!" << endl;
	}
	MPI_Finalize();
	for (int i = 0; i < n; ++i) {
        delete[] A[i];  // Delete each row
    }
    delete[] A;  // Delete the array of pointers
	return 0;
}
