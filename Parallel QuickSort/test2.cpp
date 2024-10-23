#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

void display(vector<int>& A){
	for(auto e:A){
		cout << e << " ";
	}
	cout << endl;
}

int partition(vector<int>& A,int s,int e,int p){
	//1st element is chosen as pivot
	int i = s + 1, pivot = p, j;
	// cout << "Pivot is : " << pivot << endl;
	for (j = s + 1; j < e; j++){
	 	if(A[j] < pivot){
	 		swap(A[j],A[i]);
	 		i++;
	 		//cout << i << endl;
	 	}
	 	//display(A,6);
	}
	swap(A[s],A[i-1]);
	//cout << i << endl;
	return (i-1);
	// if (pivot == A[s]) ? return (i-1) : ;
}

int main(int argc,char** argv){
	int rank,p;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int n, rem, pvt;
	vector<int> A;
	//receive buffer for each process
	vector<int> recvbuf;
	vector<int> sendCounts(p,0), displs(p,0);
	ifstream fin;
	fin.open("input.txt");
	fin >> n;
	if (rank == 0){
		A.resize(n,0);
		for (int i = 0; i < n; ++i) fin >> A[i];
		display(A);
		fin.close();

		if (n % p != 0){
			cout << "Error: Number of processors does not completely divide Size of array ." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}

		if(p != 4){
			cout << "The program is meant to be run with 4 processors " << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
	}
	
	rem = n%p;
	displs[0] = 0;
	for(int i = 0; i < p; ++i) sendCounts[i] = ((n/p)+((i < rem)?1:0));
	for(int j = 1; j < p; ++j) displs[j] = sendCounts[j-1] + displs[j-1];
	recvbuf.resize(sendCounts[rank]);

	//Scatter the array to different processes
	MPI_Scatterv(A.data(),sendCounts.data(),displs.data(),MPI_INT,recvbuf.data(),sendCounts[rank],MPI_INT,0,MPI_COMM_WORLD);

	// cout << "\nProcess " << rank << " received: " << endl;
	// for (auto e: recvbuf) cout << e <<" ,";
	// cout << endl;

	if (rank == 0) {
		pvt = recvbuf[0];
	}
	MPI_Bcast(&pvt , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
	//partition index for each element
	int prt;
	prt = partition(recvbuf,0,sendCounts[rank],pvt);

	cout << "\nProcess " << rank << " received: " << endl;
	cout << prt <<" ";
	cout << endl;

	int peer = (rank == 0) ? 2 : 0;
	MPI_Sendrecv(recvbuf.data()+prt,sendCounts[rank] , MPI_Datatype sendtype , int dest , int sendtag , void* recvbuf , MPI_Count recvcount , MPI_Datatype recvtype , int source , int recvtag , MPI_Comm comm , MPI_Status* status);	
	// quickSort(A,0,n,A[0]);
	// display(A,n);
	MPI_Finalize();
	return 0;
}