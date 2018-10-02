#include <iostream>
#include <mpi.h>
#include <cmath>
#include <iomanip>
#include <ctime>

using namespace std;

#define ROOT 0

double evalFunc1(double);
double calcChunk(int, double, double);
double getChildrenSums(MPI_Status, int, int);
int findParentId(int, int);

int main(int argc, char* argv[]){
	int id, np;
	unsigned int n;
	MPI_Status status;
	MPI_Request req, req_r;
	clock_t t1, t2;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	// ask the user for the number of parts of the integral
	if(id == 0){
		cout<<"Enter value of n: ";
		cin>>n;
	}
	MPI_Barrier(MPI_COMM_WORLD); // tell all processors to wait here until we grab n
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // send n from processor 0 to all other processors and continue

	if(n % np != 0){ // check to make sure the number of parts can be evenly spread out among the processors
		if(id == ROOT) cout<<"n/p is not evenly divisible! Exiting..."<<endl;
		MPI_Finalize();
		return 0;
	}
	
	t1 = clock();
	double step = 1.0 / (double) n; // step value of the integral, this is delta x for the trapezoidal rule
	int parts = n / np; // for each chunk (1 processor = 1 chunk), how many sub-parts of the integrals to calculate?

	if(id == ROOT){ // root node
		double totalSum = calcChunk(id, parts, step) + getChildrenSums(status, id, np);
		cout<<setprecision(9)<<"Integral value is "<<totalSum<<endl;
	} else{ // non-root node
		int parentId = findParentId(id, np);
		double totalSum = calcChunk(id, parts, step) + getChildrenSums(status, id, np);
		MPI_Send(&totalSum, 1, MPI_DOUBLE, parentId, 0, MPI_COMM_WORLD);
	}

	t2 = clock();
	if(id == ROOT){
		float diff ((float)t2-(float)t1);
		float sec = diff / CLOCKS_PER_SEC;
		cout<<"Execution time: "<<sec<<endl;
	}
	MPI_Finalize();

	return 0;
}

// This function is the function we are computing the integral over
double evalFunc1(double x){
	return 4.0 / (1.0 + pow(x, 2)); // NOTE: This function is a sample function given in the project; you may swap out to whichever function you want it will work in any situation
}

// Calculate the chunk of the integral assigned to each processor
double calcChunk(int chunkNum, double parts, double step){
	double x = chunkNum*(parts*step); // current x value we're evaluating the 'part' of the 'chunk' at
	double totalSum = 0.0;
	for(int i = 0; i < parts; i++){
		double leftVal = evalFunc1(x); // f(x_{i-1})
		double rightVal = evalFunc1(x + step); // f(x_{i})
		totalSum += (leftVal + rightVal)/2.0 * step; // trapezoidal rule
		x += step; // move x to the next 'part' of the 'chunk'
	}

	return totalSum;
}

// For each processor, find the child processor ids and receive the 'chunk' sums from them, to pass it up the binary tree collective
double getChildrenSums(MPI_Status status, int id, int np){
	double totalSum = 0.0;
	int subTreeLen = np / 2; // i.e. for 8 processors, this starts at 4
	int currParent = 0;
	while(subTreeLen >= 1 && id != ROOT){ // iterate until we find this processor's id in the "binary tree" collective
		if(id > currParent + subTreeLen) currParent = currParent + subTreeLen; // keep on going down the tree
		else if(id == currParent + subTreeLen){ // we found this node! now, we know the width of the sub-tree with all child nodes
			currParent = id;
			subTreeLen /= 2;
			break;
		}
		subTreeLen /= 2;
	}

	double sum;
	while(subTreeLen >= 1){ // for processor 8, the children are 8 + 4, 8 + 2, and 8 + 1; notice a pattern? the 4, 2, and 1 are the subTreeLen
		MPI_Recv(&sum, 1, MPI_DOUBLE, id + subTreeLen, 0, MPI_COMM_WORLD, &status); // grab the 'chunk' sum from the child processor
		subTreeLen /= 2;
		totalSum += sum;
	}

	return totalSum;
}

// For any given node, except for ROOT, find the id of the parent processor to send the sum to
int findParentId(int id, int np){
	int subTreeLen = np / 2;
	int currParent = 0;
	while(subTreeLen >= 1){ // Go down the tree until we find the parent's subtree with this processor as its immediate child
		if(id > currParent + subTreeLen) currParent = currParent + subTreeLen;
		else if(id == currParent + subTreeLen) return currParent;
		subTreeLen /= 2;
	}

	return currParent;
}
