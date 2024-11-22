#include <bits/stdc++.h>
using namespace std;
#include <mpi.h>

int main(int argc,char *argv[]){
    MPI_Init(&argc,&argv);
  	int d=2 , dims[2] = {2,2}, period[2] = {1,1}, reorder = 1, rank;
    MPI_Comm cartTopo;

    MPI_Cart_create(MPI_COMM_WORLD, d, dims, period, reorder, &cartTopo);
    MPI_Comm_rank(cartTopo,&rank);

    int n; //size of matrix
    ifstream fin;
    fin.open("./input.txt");
    fin >> n;

    if (n != 4){
        cout << "The program requires that n should be 4." << endl;
        MPI_Abort(cartTopo,EXIT_FAILURE);
	}

    int A[n][n],B[n][n],C[n][n],A_loc[n/2][n/2],B_loc[n/2][n/2],C_loc[n/2][n/2];

    if(rank == 0){
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

        if ((dims[0] * dims[1]) != 4){
			cout << "Error: Number of processors should be 4." << endl;
			MPI_Abort(cartTopo,EXIT_FAILURE);
		}
    }

    int sendCounts1[n],displs1[dims[0]*dims[1]];
    for (int i = 0; i < n; i++) sendCounts1[i] = (dims[0]*dims[1])/2; 
    for (int i = 0; i < n; i++) displs1[i] = ((i < n/2) ? sendCounts1[i] * i : sendCounts1[i] * (i+2)); // Refine the formula
    // int blockSize = n / dims[0]; // Block size for each dimension
    // for (int i = 0; i < dims[0]; ++i) {
    //     for (int j = 0; j < dims[1]; ++j) {
    //         int rank = i * dims[1] + j;?
    //         sendCounts1[rank] = blockSize * blockSize; // Elements in each block
    //         displs1[rank] = (i * blockSize * n) + (j * blockSize); // 2D displacement
    //     }
    // }

    MPI_Scatterv(&A, sendCounts1, displs1, MPI_INT, &A_loc[0], n/2, MPI_INT, 0, cartTopo);
	MPI_Scatterv(&B, sendCounts1, displs1, MPI_INT, &B_loc[0], n/2, MPI_INT, 0, cartTopo);

    int sendCounts2[n],displs2[dims[0]*dims[1]];
    for (int i = 0; i < n; i++) sendCounts2[i] = (dims[0]*dims[1])/2; 
    for (int i = 0; i < n; i++) displs2[i] = ((i < n/2) ? sendCounts1[i] * (i+2) : sendCounts1[i] * (i+4)); // Refine the formula

    MPI_Scatterv(&A, sendCounts2, displs2, MPI_INT, &A_loc[1], n/2, MPI_INT, 0, cartTopo);
	MPI_Scatterv(&B, sendCounts2, displs2, MPI_INT, &B_loc[1], n/2, MPI_INT, 0, cartTopo);
    memset(C_loc, 0, sizeof(C_loc));


    for(int step = 0; step < (n/2) ;step++)
	{   
		if(step == 0){
			if (rank == 2 || rank == 3)
			{
				MPI_Request send_request, recv_request;
				MPI_Status status;		
				int left,right;
				int shift=-1; //move left
				int change=1;
				MPI_Cart_shift(cartTopo,change,shift,&left,&right);
				MPI_Isend(&A_loc,4,MPI_INT,left,0,cartTopo,&send_request);
				MPI_Irecv(&A_loc,4,MPI_INT,right,0,cartTopo,&recv_request);
                /*MPI_Isend sends the local A_local to the left neighbor, 
                and MPI_Irecv receives A_local from the right neighbor.*/
				MPI_Wait(&send_request, &status);
				MPI_Wait(&recv_request, &status);
			}	
			if (rank == 1 || rank == 3)
			{
				MPI_Request send_request1, recv_request1;
				MPI_Status status1;		
				int up,down;
				int shift=1;
				int change=0; // dim 0
				MPI_Cart_shift(cartTopo,change,shift,&up,&down);
				MPI_Isend(&B_loc,4,MPI_INT,up,1,cartTopo,&send_request1);
				MPI_Irecv(&B_loc,4,MPI_INT,down,1,cartTopo,&recv_request1);
				MPI_Wait(&send_request1, &status1);
				MPI_Wait(&recv_request1, &status1);
			}
		}
	
        if(step != 0)
        {
            MPI_Request send_request2, recv_request2;
            MPI_Status status2;		
            int left1,right1;
            int shift1=-1;
            int change1=1;
            MPI_Cart_shift(cartTopo,change1,shift1,&left1,&right1);
            MPI_Isend(&A_loc,4,MPI_INT,left1,0,cartTopo,&send_request2);
            MPI_Irecv(&A_loc,4,MPI_INT,right1,0,cartTopo,&recv_request2);
            MPI_Wait(&send_request2, &status2);
            MPI_Wait(&recv_request2, &status2);

            MPI_Request send_request3, recv_request3;
            MPI_Status status3;		
            int up1,down1;
            int shift2=1;
            int change2=0;
            MPI_Cart_shift(cartTopo,change2,shift2,&up1,&down1);
            MPI_Isend(&B_loc,4,MPI_INT,up1,1,cartTopo,&send_request3);
            MPI_Irecv(&B_loc,4,MPI_INT,down1,1,cartTopo,&recv_request3);
            MPI_Wait(&send_request3, &status3);
            MPI_Wait(&recv_request3, &status3);
         }

        for(int i=0;i<(n/2);i++)
        { 
            for (int j=0;j<(n/2);j++)
            { 
                for (int k=0;k<2;k++)
                {
                    C_loc[i][j] =C_loc[i][j]+ A_loc[i][k]*B_loc[k][j];
                }
            }
        }
	}

    MPI_Gatherv(&C_loc[0],2, MPI_INT,&C,sendCounts1, displs1, MPI_INT, 0, cartTopo);
	MPI_Gatherv(&C_loc[1],2, MPI_INT,&C,sendCounts2, displs2, MPI_INT, 0, cartTopo);

    // Print the final result matrix on rank 0
    if (rank == 0) {
        ofstream fout;
        fout.open("./output.txt");
        fout << "Final Resultant Matrix C:" << endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                fout << C[i][j] << " ";
            }
            fout << endl;
        }
        cout << "\nPassed the resultant matrix onto output file" << endl;
		fout.close();
    }

    MPI_Finalize();
}