#include<mpi.h>
#include<bits/stdc++.h>
using namespace std;


double f(double x){
	return 0.5*sin(x);
}

double fivep_integration(int x_k,double h,double A[]){
	return ((2*h)/(45))*(double)(7*f(A[4*x_k]) + 32*f(A[4*x_k+1]) + 12*f(A[4*x_k+2]) + 32*f(A[4*x_k+3]) + 7*f(A[4*x_k+4]));
}
//(k) is assumed to be divisible by no of processors (p) i.e. size here!
//Modification added to make it work when p is not a factor of k

double integral = 0;

int main(int argc,char** argv)
{ 
	int *rank = new int,size;
	int a = 0,b = 1;
	int n = 200, k = n/4;
	double h = (b-a)/(double)n;
	double sum = 0;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);

	double A[n+1];
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
			}			
		}
	}
	if (*rank == 0){
		integral += sum;
		double temp;
		for(int i = 1; i < size; ++i){
			MPI_Recv(&temp,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			integral += temp;
		}
        //case of remaining elements
        if(k % size != 0){
            sum = 0;
            for(int j = k-(k%size); j < k; ++j){
                sum += fivep_integration(j,h,A);
            }
            integral += sum;
        }
        ofstream fout;
        fout.open("./boolsOutput.txt");
		fout << "Final Approximate Integration Value :- " << integral << endl; 
        fout.close();
        cout << "Program Results passed onto to output file!" << endl;
	}
	
	else{
		MPI_Send(&sum,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
	delete rank;
	return 0;
}