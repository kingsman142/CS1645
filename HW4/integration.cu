#include <iostream>
#include <cmath>
#include <iomanip>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

using namespace std;

#define ROOT 0

__device__ double eval_func1(double);
__global__ void calc_chunk(double, double*, int);

int main(int argc, char* argv[]){
	double integral = 0.0;
	int n, numBlocks, numThreads; // n = number of points, numBlocks = blocks on the device, numThreads = threads per block
	double step; // step size when computing the integral; as n increases, step decreases and the integral's precision increases
	int parts; // number of parts of the integral to compute for each thread
	double* deviceGlobalSum;

	cout<<"Enter number of points: ";
	cin >> n;
	cout<<"Enter number of blocks (up to 28): "; // GTX1080 on Pitt's CRC only allows up to 28 processors
	cin >> numBlocks;
	cout<<"Enter number of threads per block (up to 1024): "; // GTX1080 on Pitt's CRC only allows up to 1024 threads/block
	cin >> numThreads; 

	if(n % (numBlocks*numThreads) != 0){ // can't split n up evenly by the number of threads
		cout<<"N is not evenly divisible by # of Blocks! Exiting..."<<endl;
		return 1;
	} else if(numBlocks > 28 || numThreads > 1024){
		cout<<"Too many requested blocks or threads! Exiting..."<<endl;
		return 1;
	}

	step = 1.0 / (double) n; // step value of the integral, this is delta x for the trapezoidal rule
	parts = (n / numBlocks) / numThreads; // number of parts per thread

	// begin the timer for the multi-threaded portion of the code
	cudaEvent_t start, end;
	cudaEventCreate(&start);
	cudaEventCreate(&end);
	cudaEventRecord(start);

	// copy the global sum, or future integral value, from the host to the device
	cudaMalloc((void**) &deviceGlobalSum, sizeof(double));
	cudaMemcpy(deviceGlobalSum, &integral, sizeof(double), cudaMemcpyHostToDevice);

	// call a function with specified number of blocks (i.e. address spaces, similar to processors) and number of threads (i.e. threads for each block)
	calc_chunk<<<numBlocks, numThreads>>>(step, deviceGlobalSum, parts);

	// barrier for all devices to complete
	cudaDeviceSynchronize();

	// copy the global sum of the device to the integral variable on the host
	cudaMemcpy(&integral, deviceGlobalSum, sizeof(double), cudaMemcpyDeviceToHost);

	// compute the total computation time
	cudaEventRecord(end);
	float timeDiff;
	cudaEventElapsedTime(&timeDiff, start, end);
	cudaEventSynchronize(end);
	cudaEventDestroy(start);
	cudaEventDestroy(end);

	// free memory malloc'd on the device
	cudaFree(deviceGlobalSum);

	// print out the time and integral value
	cout<<setprecision(9)<<"Final value of the integral is "<<integral<<"!"<<endl;
	cout<<setprecision(9)<<"Execution time: "<<timeDiff<<endl;

	return 0;
}

// Integral function we want to evaluate; increases modularity of the program
__device__ double eval_func1(double x){
	return 4.0 / (1.0 + pow(x, 2));
}

// Calculate the chunk of the integral assigned to each thread (which belongs to a block)
__global__ void calc_chunk(double step, double* globalSum, int parts){
	__shared__ double blockSum; // sum of all the parts in this block's chunk
	if(threadIdx.x == 0) blockSum = 0.0; // let thread 0 initilize the shared variable since you can't declare and initialize on the same line
	__syncthreads(); // barrier so no thread writes to the shared variable before it's initialized

	int chunkNum = blockIdx.x * blockDim.x + threadIdx.x; // the chunk of the integral this thread is working on
	double threadSum = 0.0;
	double x = chunkNum * step * parts; // current x value we're evaluating the 'part' of the 'chunk' at
	for(int i = 0; i < parts; i++){
		double leftVal = eval_func1(x); // f(x_{i-1})
		double rightVal = eval_func1(x + step); // f(x_{i})
		threadSum += ((leftVal + rightVal)/2.0) * step; // trapezoidal rule
		x += step;
	}
	atomicAdd(&blockSum, threadSum); // all thread write to the shared variable of the block

	__syncthreads();
	if(threadIdx.x == ROOT){ 
		atomicAdd(globalSum, blockSum); // all blocks write to the 'shared' variable, global sum, of the device
	}
}
