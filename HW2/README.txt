Background
==========
This project was written for my High Performance Computing class.  For this project, we are computing the integral over the function 4 / (1 + x^2) on the interval from 0 to 1.
The final value of the integral should be pi, which is actually what it ends up being.
Additionally, customCollective.cpp is ran.  In this program, our own binary tree collective algorithm is ran to spread work out evenly across processors.  For example, with 8 processors, each processor has a parent node and children.  The parents will collect integral values from their child nodes, sum them up, add that sum to their integral chunk value, and push their data up to their parent until they all gather in the root, or processor 0.  Once all integral chunk values have arrived at the root, they are finally totalled to the actual integral value.

This program showcases some of the key strong points and benefits of parallel processing, even though this is on a small scale.  

Tools/Requirements
==================

OpenMPI 3.1.2 ( https://www.open-mpi.org/software/ompi/v3.1/ )

Setup/Run instructions
==================
Begin on a machine with a Unix distribution
Load the OpenMPI module ( `module load openmpi` from bash )
Run `mpic++ customCollective.cpp` to compile the program with OpenMPI
Run `mpiexec -np NUM_PROC a.out`
	- NUM_PROC is the number of processors to run the program with - replace this with a number; ensure this is a power of 2
		- Furthermore, NUM_PROC cannot be more than the number of processors in your machine.  For example, most machines have 4 cores
		- As NUM_PROC increases, the time spent per processor on the program decreases, but the total sys time increases since there is more work required to send to parents and gather from children
	- When prompted, enter the value of N that you desire.  As N increases, the accuracy of the integral's computed value increases.  With N = 2^20 (1,048,576), the precision is around 7 digits of pi.
Once all user input has been processed, the integral is calculated and output to the user

Notes
=====
In the above setup, we compiled customCollective.cpp and ran the program.  This program contains a custom implemented, created by me, to implement a binary tree distributed system program.  However, if we replace customCollective.cpp with mpiNativeCollective.cpp, we can see the ground truth running times of the binary tree collective that is implemented by OpenMPI.  So, you are able to compare my implementation to theirs, which turns out to be slightly slower, but very much comparable.