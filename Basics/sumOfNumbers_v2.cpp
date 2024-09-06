#include<mpi.h>
#include <bits/stdc++.h>
using namespace std;
int total = 0;

int main(int argc,char** argv)
{ 
	int *rank = new int,size,i=0;
	ifstream fin;
    fin.open("input.txt");
    int len = 12;
    fin >> len;
    int A[len];
    for (int p = 0; p < len; ++p) fin >> A[p];
    fin.close();
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
	
	int ofs = 0;
	int s = 0,offset = len/size;
	int B[offset];
	
	if(*rank == 0){
		for(int i = 1; i < size; i++){
			ofs = i * offset;
			MPI_Send(&A[ofs],offset,MPI_INT,i,1,MPI_COMM_WORLD);
		}

		//handling remaining parts
		if(len % size != 0){
			int count = 1;
			while(count < (len%size)+1){
				MPI_Send(&A[size*offset+(count-1)],1,MPI_INT,count,2,MPI_COMM_WORLD);
				count++;
			}
		}

		//P0's part..
		for(int i = 0; i < offset; ++i){
			s += A[i];
		}
		total += s;

		//sum computation
		for(int j = 1; j < size; ++j)
		{ 
			int temp = 0;
			MPI_Recv(&temp,1,MPI_INT,j,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			total += temp;
		}
		//receive remaining elemants as well in case any
		if(len % size != 0){
			for(int k = 1; k < (len%size)+1; k++){
				int temp = 0;
				MPI_Recv(&temp,1,MPI_INT,k,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				total += temp;
			}
		}
        ofstream fout;
        fout.open("./output.txt");
		fout << "TOTAL Sum = " << total << endl; 
        fout.close();
        cout << "Program Results passed onto to output file!" << endl;
	}
	else{
		MPI_Recv(&B,offset,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		for(int k = 0; k < offset; k++){
			s += B[k];
		}
		//cout << "Rank Sum = " << s << endl;
		MPI_Send(&s,1,MPI_INT,0,1,MPI_COMM_WORLD);

		
		//handle remaining elements...stores in var t
		if(*rank >= 1 && *rank < (len%size)+1){
			int t = 0;
			MPI_Recv(&t,1,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			MPI_Send(&t,1,MPI_INT,0,2,MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
	delete rank;
	return 0;
}