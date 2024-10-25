#include <bits/stdc++.h>
#include<mpi.h>
using namespace std;

void printMatrix(const double* A, int n){
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n + 1; ++j){
			cout << A[i*(n+1)+j] << "    ";
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
    int n;
	ifstream fin;
    fin.open("./input.txt");
    fin >> n;
	double A[n][n+1];
	double reducedA[n][n+1];
	if (rank == 0){
         for (int i = 0; i < n; ++i){
            for(int j = 0; j < n+1 ; ++j){
                fin >> A[i][j];
            }
        }
        int m = n+1;
        double flatA[n*m];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                flatA[i * m + j] = A[i][j];
            }
        }
		printMatrix(flatA,n);

        if (p != 4){
			cout << "Error: Number of processors should be 4." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
	}
    int rows = n/p;
    double B[(n+1)*rows],temp[(n+1)*rows]; // for each processor
	MPI_Scatter(&A,(n+1)*rows,MPI_DOUBLE,&B,(n+1)*rows,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//A will receive the first row in B i.e R0
	if (rank == 0) memcpy(&temp,&B,sizeof(B));
	
	// Perform required number of operations
	for (int k = 1; k < n; k ++){
		
		//Send necesary Row to all other processes
		MPI_Bcast(&temp,(n+1)*rows, MPI_DOUBLE, k-1, MPI_COMM_WORLD);
		if (rank > k-1){
			Op(B,temp,(n+1),k-1); //DRY run
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		//MPI_Barrier blocks all MPI processes in the given communicator until they all call this routine.
		//copy the pivot row for next iteration
		if (rank == k && rank != n-1) memcpy(&temp,&B,sizeof(B));
	}
 	MPI_Gather(&B,(n+1)*rows,MPI_DOUBLE,&reducedA[rank],(n+1)*rows,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (rank==0){
		// printMatrix(reducedA,4);
		//Backsubstitution *
		double ans[n];
		for (int i = n-1; i >= 0; --i){
			double tmp = reducedA[i][n];
			for(int j = i+1; j < n; ++j){
				tmp -= reducedA[i][j]*ans[j];
			}
			ans[i] = tmp/reducedA[i][i];
		}

		ofstream fout;
        fout.open("./output.txt");
		for (int i = 0; i < n; ++i){
			fout << "\nx"<<i+1 << ": " << ans[i] << endl;  
		}
		cout << "Solution passed to the output file!" << endl;
		fout << endl;
		fout.close();
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
