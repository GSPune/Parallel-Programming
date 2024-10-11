#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

void display(vector<int>& A,int size){
	for(auto e:A){
		cout << e << " ";
	}
	cout << endl;
}

int partition(vector<int>& A,int s,int e){
	//1st element is chosen as pivot
	int i = s + 1, pivot = A[s], j;
	cout << "Pivot is : " << pivot << endl;
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
}

void quickSort(vector<int>& A,int s,int e){
	//cout << "In here.. " << s << ".." << e << endl;
	if(s < e){
		int k = partition(A,s,e);
		quickSort(A,s,k-1);
		quickSort(A,k+1,e);
	}
}

int main(int argc,char** argv){
	int rank,p;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int n = 6;
	// int A[n] = {7,12,5,6,2,8};
	vector<int> A = {7,12,5,6,2,8};
	display(A,6);
	quickSort(A,0,n);
	display(A,n);
	MPI_Finalize();
	return 0;
}