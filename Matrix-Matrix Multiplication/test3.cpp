#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

int main(int argc,char** argv){
	int rank,p;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int n; //size of matrix
    ifstream fin;
    fin.open("./input.txt");
    fin >> n;
    int A[n][n],B[n][n];
    int submSize = int(n/sqrt(p));
    int recvbufA[submSize*submSize];
    int recvbufB[submSize*submSize];
    // cout << n/sqrt(p) << endl;
    if (rank == 0){
        for (int i = 0; i < n; ++i){
            for(int j = 0; j < n ; ++j){
                fin >> A[i][j];
            }
        }

        for (int i = 0; i < n; ++i){
            for(int j = 0; j < n ; ++j){
                fin >> B[i][j];
            }
        }

        fin.close();

        if (p != 4){
			cout << "Error: Number of processors should be 4." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}

        for (int i = 0; i < n; ++i){
            for(int j = 0; j < n ; ++j){
                cout << B[i][j] << " ";
            }
            cout << endl;
        }
    }

    int sendCounts[p], displs[p];
    displs[0] = 0;
    for (int i = 0; i < p; i++) sendCounts[i] = p;
    for (int i = 0; i < p; i++) displs[i] = ((i < p/2) ? p * i : p * (i+p+2)); // Refine the formula

    // if (rank == 0) for (int i = 0; i < p; i++) { cout << displs[i] << " \n";}
    for (int k = 0; k < p; k++){
        MPI_Scatterv(&A, sendCounts, displs ,MPI_INT,&recvbufA[k*p],p,MPI_INT,0,MPI_COMM_WORLD);
        for (int i = 0; i < p; i++) {displs[i] += n;}
    }

    for (int i = 0; i < p; i++) displs[i] = ((i < p/2) ? p * i : p * (i+p+2));
    for (int k = 0; k < p; k++){
        MPI_Scatterv(&B, sendCounts, displs ,MPI_INT,&recvbufB[k*p],p,MPI_INT,0,MPI_COMM_WORLD);
        for (int i = 0; i < p; i++) {displs[i] += n;}
    }

    cout << "\nProcess " << rank << " received: " << endl;
	for (int i = 0; i < submSize; ++i){
        for (int j = 0; j< submSize; ++j){
            cout << recvbufB[i*submSize+j] << " ";
        }
        cout << endl;
    }

    
	MPI_Finalize();
	return 0;
}