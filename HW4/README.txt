Background
==========
This project was written for a High Performance Computing class.  For this project, we are computing the integral over the function 4 / (1 + x^2) on the interval from 0 to 1.
The final value of the integral should be pi, which is actually what it ends up being.  The power of blocks and threads, through the use of the CUDA framework, is displayed with the use of GPUs.

Tools/Requirements
==================

C
CUDA 9.0

Setup/Run instructions
==================
Begin on a machine with a Unix distribution
Run `module load cuda` to load CUDA
Run `nvcc -arch=sm_60 integration.cu` to compile the program
Run `./a.out`
	- When prompted, enter the number of desired threads.  As the number of blocks and threads per block increases, computation time decreases.
	- When prompted, enter the value of N that you desire.  As N increases, the accuracy of the integral's computed value increases.  With N = 2^20 (1,048,576), the precision is around 7 digits of pi.
Once all user input has been processed, the integral is calculated and output to the user

NOTE: If you want to perform GPU profiling on the program, run `nvprof ./a.out` instead of just `./a.out` to execute the program.