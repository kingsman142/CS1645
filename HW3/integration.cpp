#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <pthread.h>

using namespace std;

double eval_func1(double);
void* calc_chunk(void* data);

pthread_mutex_t guard = PTHREAD_MUTEX_INITIALIZER; // global mutex guard
double step; // step size when computing the integral; as n increases, step decreases and the integral's precision increases
int parts; // number of subsections for each thread's 'chunk'

struct chunk_data{
	int chunkNum;
	double* totalSum;
};

int main(int argc, char* argv[]){
	double integral = 0.0;
	int numThreads, n;

	cout<<"Enter number of threads: ";
	cin >> numThreads;
	cout<<"Enter value of n: ";
	cin >> n;

	if(n % numThreads != 0){ // can't split n up evenly by the number of threads
		cout<<"N is not evenly divisible by P! Exiting..."<<endl;
		return 1;
	}

	pthread_t threads[numThreads];
	struct chunk_data datas[numThreads];

	double sums[numThreads]; // stores the sum of each thread's integral chunk
	for(int i = 0; i < numThreads; i++) sums[i] = 0.0; // initialize to 0.0 just to be safe

	// set important parameters to compute the integral
	step = 1.0 / (double) n; // step value of the integral, this is delta x for the trapezoidal rule
	parts = n / numThreads; // for each chunk (1 processor = 1 chunk), how many sub-parts of the integrals to calculate?

	// begin the timer for the multi-threaded portion of the code
	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);

	// create all threads, telling each one to run calc_chunk(...)
	for(int i = 0; i < numThreads; i++){
		datas[i].chunkNum = i;
		datas[i].totalSum = &sums[i];
		pthread_create(&threads[i], NULL, calc_chunk, &datas[i]);
	}

	// join all threads back again
	for(int i = 0; i < numThreads; i++) pthread_join(threads[i], NULL);
	
	// compute the total computation time
	clock_gettime(CLOCK_MONOTONIC, &end);
	float timeDiff = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_sec)/1000000000.0;

	// sum up their respective integral chunks
	for(int i = 0; i < numThreads; i++) integral += sums[i];
	cout<<setprecision(9)<<"Final value of the integral is "<<integral<<"!"<<endl;
	cout<<setprecision(9)<<"Execution time: "<<timeDiff<<endl;

	return 0;
}

// Integral function we want to evaluate; increases modularity of the program
double eval_func1(double x){
	return 4.0 / (1.0 + pow(x, 2));
}

// Calculate the chunk of the integral assigned to each processor
void* calc_chunk(void* data){
	struct chunk_data chunkData = *(struct chunk_data*) data;
	int chunkNum = chunkData.chunkNum;
	double* totalSum = chunkData.totalSum;

	double x = chunkNum*parts*step; // current x value we're evaluating the 'part' of the 'chunk' at
	for(int i = 0; i < parts; i++){
		double leftVal = eval_func1(x); // f(x_{i-1})
		double rightVal = eval_func1(x + step); // f(x_{i})
		*(totalSum) += (leftVal + rightVal)/2.0 * step; // trapezoidal rule
		x += step; // move x to the next 'part' of the 'chunk'
	}
}
