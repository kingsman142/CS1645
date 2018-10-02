#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

double eval_func1(double);

int main(int argc, char* argv[]){
	MPI_Status status;
	MPI_Request req, req_r;
	int count, rank;
	double integral = 0.0;
	double* arr;

	MPI_Init(&argc, &argv); //Initiate computation
	MPI_Comm_size(MPI_COMM_WORLD, &count); //Find # of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Find my id

	double sumTerm = 0.0; // Term of the trapezoid rule sum for each processor
	arr = new double[count];

	if(rank > 0){
		double step = 1.0 / (double) (count-1); // Using count-1 processors to evaluate integral
		double leftVal = eval_func1((rank-1) * step); // Evaluate the function at the given value
		double rightVal = eval_func1((rank) * step); // Evaluate the function at the given value
		sumTerm = ((leftVal + rightVal) / 2.0) * step;
	}

	MPI_Gather(&sumTerm, 1, MPI_DOUBLE,
		arr, 1, MPI_DOUBLE,
		0, MPI_COMM_WORLD);

	if(rank == 0){
		for(int i = 0; i < count; i++){
			integral += arr[i];
		}
		cout<<"Final value of the integral is "<<integral<<"!"<<endl;
	}

	delete[] arr;

	MPI_Finalize();

	return 0;
}

double eval_func1(double x){
	return 4.0 / (1.0 + pow(x, 2));
}
