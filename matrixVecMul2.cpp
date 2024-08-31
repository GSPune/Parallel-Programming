#include<mpi.h>
#include <bits/stdc++.h>
using namespace std;
#include <cmath>
//Distributed System
int main(int argc,char** argv)
{ 
	int rank,p;
	int n;
    ifstream fin;
    fin.open("inputMat.txt");
    fin >> n;
	// int A[n][n] = {{1,5,0,4},{2,6,5,0},{0,7,3,2},{4,0,1,1}};
    int A[n][n];
    for(int i = 0 ; i < n ; ++i){
        for (int j = 0 ; j < n; ++j){
            fin >> A[i][j];
        }
    }
	int x[n];
    for (int k = 0; k < n; ++k) fin >> x[k];
    fin.close();
    MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int rem = n % p;
    int sendCounts[p], displs[p];
    displs[0] = 0;
    for(int i = 0; i < p; ++i) sendCounts[i] = ((n/p)+((i < rem)?1:0))*n;
    for(int j = 1; j < p; ++j) displs[j] = sendCounts[j-1] + displs[j-1];
    int recvbuf[sendCounts[rank]];
    //Distributing rows of A

    MPI_Scatterv(&A,sendCounts,displs,MPI_INT,&recvbuf,sendCounts[rank],MPI_INT,0,MPI_COMM_WORLD);
    //send sCnt elements of x to each process
    int scnt2[p],disp[p];disp[0] = 0;
    for(int i = 0; i < p; ++i) scnt2[i] = ((n/p)+((i < rem)?1:0));
    for(int j = 1; j < p; ++j) disp[j] = scnt2[j-1] + disp[j-1];
	int rcv[scnt2[rank]];
    MPI_Scatterv(&x,scnt2,disp,MPI_INT,&rcv,sendCounts[rank],MPI_INT,0,MPI_COMM_WORLD);
	// MPI_Scatter(&x,sCnt,MPI_INT,rcv,sCnt,MPI_INT,0,MPI_COMM_WORLD);
    //Gather all xi's on each process
    int xGathered[n];
	MPI_Allgatherv(rcv,scnt2[rank],MPI_INT,xGathered,scnt2,disp,MPI_INT,MPI_COMM_WORLD);

    int y[scnt2[rank]];
    for(int j = 0; j < scnt2[rank]; ++j) y[j] = 0;

    for(int k = 0; k < sendCounts[rank] ; ++k){ 
        y[(k/n)] += recvbuf[k]*xGathered[k%n];
    }

    int bGathered[n];
    MPI_Gatherv(&y,scnt2[rank],MPI_INT,bGathered,scnt2,disp,MPI_INT,0,MPI_COMM_WORLD);
	if(rank == 0){
        ofstream fout;
        fout.open("./output.txt");
        fout << "The computed components of the y resultant matrix are given below:-" << endl;
		for(int i = 0; i < n ; ++i){
			fout << "y" << i << " is " << bGathered[i] << endl;
		}
        fout.close();
        cout << "Program Results passed onto to output file!" << endl;
	}
	MPI_Finalize();
	return 0;
}

//---Debug Process----
 // if(rank == 0){
//     //for(int i = 0; i < p; ++i) cout << scnt2[i] << ", ";
//     for(int i = 0; i < p; ++i) cout << disp[i] << ", ";
// }
// cout << "All Elements gathered on proc " << rank << endl;
// for(int i = 0; i < n; ++i) cout << xGathered[i] << ", ";
// cout << endl;