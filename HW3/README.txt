Background
==========
This project was written for a High Performance Computing class.  For this project, we are computing the integral over the function 4 / (1 + x^2) on the interval from 0 to 1.
The final value of the integral should be pi, which is actually what it ends up being.  The power of multithreading, through the use of the C/C++ library pthreads, is displayed.

Tools/Requirements
==================

C/C++
gcc/g++/other compiler

Setup/Run instructions (for C++)
==================
Begin on a machine with a Unix distribution
Run `g++ integration.cpp -lpthread` to compile the program
Run `./a.out`
	- When prompted, enter the number of desired threads.  As the number of threads increases, computation time decreases.  A significant increase is not easily shown until the number of threads approaches the thousands.
	- When prompted, enter the value of N that you desire.  As N increases, the accuracy of the integral's computed value increases.  With N = 2^20 (1,048,576), the precision is around 7 digits of pi.
Once all user input has been processed, the integral is calculated and output to the user