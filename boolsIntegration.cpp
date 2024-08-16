#include<mpi.h>
#include<iostream>
#include<vector>
#include <cmath>
using namespace std;


double f(double x){
	return sin(x);
}

double fivep_integration(int x_k,double h,double A[]){
	return ((2*h)/(45))*(double)(7*f(A[4*x_k]) + 32*f(A[4*x_k+1]) + 12*f(A[4*x_k+2]) + 32*f(A[4*x_k+3]) + 7*f(A[4*x_k+4]));
}
//(k) is assumed to be divisible by no of processors i.e. size here!

double integral = 0;

int main(int argc,char** argv)
{ 
	int *rank = new int,size;
	int a = 0,b = 1;
	int n = 120, k = n/4;
	double h = (b-a)/(double)n;
	double sum = 0;
	//if (*rank == 0) cout << "h::" << (2*h)/45 << endl;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);

	double A[n+1];//A[0]= 0;
	for (int i = 0; i < n+1 ; ++i){
		A[i] = (double)(a + i*h);
		//if (*rank == 0) cout << A[i] << endl ;
	}

	int iters = k/size;
	//if (*rank == 0) cout << "iters is " << iters << endl;
	//for loop to distribute to each processor
	for (int i = 0; i < size; ++i){
		if(*rank == i){
			int st = (*rank) * iters;
			sum = 0;
			//loop to actually do the work!
			for(int k = st; k < st + iters; k++){
				
				sum += fivep_integration(k,h,A);
				//cout << fivep_integration(k,h,A) << endl;
			}
			cout << "rank and sum are " << *rank <<","<< sum << endl;
			
		}
	}
	if (*rank == 0){
		integral += sum;
		double temp;
		for(int i = 1; i < size; ++i){
			MPI_Recv(&temp,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			integral += temp;
		}
		cout << "final integration value :- " << integral << endl;
	}
	
	else{
		MPI_Send(&sum,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
	delete rank;
	return 0;
}