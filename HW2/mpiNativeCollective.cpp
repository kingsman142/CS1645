#include <iostream>
#include <cmath>
#include <mpi.h>
#include <ctime>
#include <iomanip>

using namespace std;
#define ROOT 0

double eval_func1(double);
double trapezoid_rule(double, double, double);
double calc_chunk(int, double, double);

int main(int argc, char* argv[]){
	MPI_Status status;
	MPI_Request req, req_r;
	int count, rank, n;
	double integral = 0.0;

	MPI_Init(&argc, &argv); //Initiate computation
	MPI_Comm_size(MPI_COMM_WORLD, &count); //Find # of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Find my id

	// ask the user for the number of parts of the integral
	if(rank == 0){
		cout<<"Enter value of n: ";
		cin>>n;
	}
	MPI_Barrier(MPI_COMM_WORLD); // tell all processors to wait here until we grab n
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // send n from processor 0 to all other processors and continue

	if(n % count != 0){ // check to make sure the number of parts can be evenly spread out among the processors
		if(rank == ROOT) cout<<"n/p is not evenly divisible! Exiting..."<<endl;
		MPI_Finalize();
		return 0;
	}
	
	clock_t t1 = clock();
	
	double parts = n / count;
	double step = 1.0 / (double) n;
	double sumTerm = calc_chunk(rank, parts, step);

	MPI_Reduce(&sumTerm, &integral, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD); // sum up the results of all the processors at the ROOT
	if(rank == ROOT){
		cout<< setprecision(9) <<"Final value of the integral is "<<integral<<"!"<<endl;

		clock_t t2 = clock();
		double diff = ((double) t2 - (double) t1) / CLOCKS_PER_SEC;
		cout<<"Time taken: "<<diff<<endl;
	}

	MPI_Finalize();

	return 0;
}

// Function provided in the provide; makes it easy to swap out in the future
double eval_func1(double x){
	return 4.0 / (1.0 + pow(x, 2));
}

// Use the trapezoid rule to estimate the value of the integral
double trapezoid_rule(double leftVal, double rightVal, double step){
	return ((leftVal + rightVal) / 2.0) * step; // Term of the trapezoid rule sum for each processor
}

// Calculate the chunk of the integral assigned to each processor
double calc_chunk(int chunkNum, double parts, double step){
	double x = chunkNum*(parts*step); // current x value we're evaluating the 'part' of the 'chunk' at
	double totalSum = 0.0;
	for(int i = 0; i < parts; i++){
		double leftVal = eval_func1(x); // f(x_{i-1})
		double rightVal = eval_func1(x + step); // f(x_{i})
		totalSum += trapezoid_rule(leftVal, rightVal, step); // trapezoidal rule
		x += step; // move x to the next 'part' of the 'chunk'
	}

	return totalSum;
}
