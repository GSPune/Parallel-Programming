#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

void display(const int* A, int s){
    cout << endl;
	for (int i = 0; i < s; ++i){
        for (int j = 0; j< s; ++j){
            cout << A[i*s+j] << "  ";
        }
        cout << endl;
    }
    cout << endl;
}

void matrix_multiply(const int* A, const int* B, int* C, int subm) {
    for (int i = 0; i < subm; ++i) {
        for (int j = 0; j < subm; ++j) {
            C[i * subm + j] = 0;
            for (int k = 0; k < subm; ++k) {
                C[i * subm + j] += A[i * subm + k] * B[k * subm + j];
            }
        }
    }
}

void matrix_add(const int* A, const int* B, int* C, int subm) {
    for (int i = 0; i < subm; ++i) {
        for (int j = 0; j < subm; ++j) {
            // Calculate the index in the 1D array representation
            int index = i * subm + j;
            C[index] = A[index] + B[index];
        }
    }
}


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
    int Result[n*n];
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

    //horizontal data sharing..
    int recvbufTA[submSize*submSize]; //buffer for new submatrix that comes horizontally for each processor
    int recvbufTB[submSize*submSize];

    if(rank == 0){ //processor 0 sends to proc. 1 & vice versa
        MPI_Sendrecv(&recvbufA, submSize*submSize ,MPI_INT ,1 ,0 ,&recvbufTA ,submSize*submSize ,MPI_INT ,1 ,0 ,MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
    }
    else if (rank == 1){
        MPI_Sendrecv(&recvbufA, submSize*submSize ,MPI_INT ,0 ,0 , &recvbufTA ,submSize*submSize ,MPI_INT ,0 ,0 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(rank == 2){
        MPI_Sendrecv(&recvbufA, submSize*submSize , MPI_INT , 3 , 1 , &recvbufTA , submSize*submSize , MPI_INT , 3 , 1 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (rank == 3){
        MPI_Sendrecv(&recvbufA, submSize*submSize , MPI_INT , 2 , 1 , &recvbufTA , submSize*submSize , MPI_INT , 2 , 1 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

        
    //vertical data sharing.. of B's distributed blocks
    if(rank == 0){
        MPI_Sendrecv(&recvbufB, submSize*submSize ,MPI_INT ,2 ,0 ,&recvbufTB ,submSize*submSize ,MPI_INT ,2 ,0 ,MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
    }
    else if (rank == 2){
        MPI_Sendrecv(&recvbufB, submSize*submSize ,MPI_INT ,0 ,0 ,&recvbufTB ,submSize*submSize ,MPI_INT ,0 ,0 ,MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
    }

    if(rank == 1){
         MPI_Sendrecv(&recvbufB, submSize*submSize ,MPI_INT ,3 ,0 ,&recvbufTB ,submSize*submSize ,MPI_INT ,3 ,0 ,MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
    }
    else if (rank == 3){
         MPI_Sendrecv(&recvbufB, submSize*submSize ,MPI_INT ,1 ,0 ,&recvbufTB ,submSize*submSize ,MPI_INT ,1 ,0 ,MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
    }
    
    //Matrix Multiplication
    int C[submSize*submSize]; // resultant submatrix for each processor
    int intdC1[submSize*submSize],intdC2[submSize*submSize]; //Intermediate arrays
    if (rank == 0){ //C0
        matrix_multiply(recvbufA,recvbufB,intdC1,submSize);
        // display(intdC1,submSize);
        matrix_multiply(recvbufTA,recvbufTB,intdC2,submSize);
        matrix_add(intdC1,intdC2,C,submSize);
    }
    else if(rank == 1){ //C1
        matrix_multiply(recvbufTA,recvbufB,intdC1,submSize);
        matrix_multiply(recvbufA,recvbufTB,intdC2,submSize);
        matrix_add(intdC1,intdC2,C,submSize);
    }
    else if(rank == 2){
        matrix_multiply(recvbufA,recvbufTB,intdC1,submSize);
        matrix_multiply(recvbufTA,recvbufB,intdC2,submSize);
        matrix_add(intdC1,intdC2,C,submSize);
    }
    else if(rank == 3){
        matrix_multiply(recvbufA,recvbufB,intdC1,submSize);
        matrix_multiply(recvbufTA,recvbufTB,intdC2,submSize);
        matrix_add(intdC1,intdC2,C,submSize);
    }
    
    MPI_Gather(&C,submSize*submSize,MPI_INT,&Result,submSize*submSize,MPI_INT,0,MPI_COMM_WORLD);
    if(rank == 0){
        int formattedMatr[n][n];
        int count = 0;
        for (int i = 0; i < submSize; ++i){
            count = 0;
            for(int j = 0; j < submSize; j++){
                for(int k = 0; k < submSize; ++k){
                    int a = (i == 0 || i == 1) ? j:j+submSize;
                    int b = (i == 0 || i == 2) ? k:k+submSize ;    
                    formattedMatr[a][b] = Result[i*(submSize*submSize) + (count++)];
                }
            }
        }
        ofstream fout;
        fout.open("./output.txt");
        for(int j = 0; j < n; j++){
                for(int k = 0; k < n; ++k){
                    fout << formattedMatr[j][k] << "  ";
                }
                fout << endl;
        }
        fout << endl;
        cout << "\nPassed the resultant matrix onto output file" << endl;
		fout.close();
    }
	MPI_Finalize();
	return 0;
}