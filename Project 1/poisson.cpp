#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

#define vvd vector< vector<double> >

vector<double> subtract_vectors(vector<double>, vector<double>);
double calc_norm(vector<double>, vector<double>, int);
void solve_jacobi();
vvd solve_gauss_seidel(int, vvd, int, int, double, double);
vvd initialize_arr(int, int);
vvd set_boundary_values(vvd, int, int, double, double, vvd l(vvd, int, double), vvd t(vvd, int, double), vvd r(vvd, int, int, double), vvd b(vvd, int, int, double));
vvd verif_left(vvd, int, double);
vvd verif_top(vvd, int, double);
vvd verif_right(vvd, int, int, double);
vvd verif_bottom(vvd, int, int, double);
void print_solution(vvd);

int main(int argc, char* argv[]){
	int x, y;

	cout<<"Enter x dimension: ";
	cin >> x;
	cout<<"Enter y dimension: ";
	cin >> y;

	double xConstraint = 2.0;
	double yConstraint = 1.0;
	double dx = (double) xConstraint / x;
	double dy = (double) yConstraint / y;

	vvd poissonMatrix = initialize_arr(y, x);
	poissonMatrix = set_boundary_values(poissonMatrix, y, x, dx, dy, verif_left, verif_top, verif_right, verif_bottom);
	print_solution(poissonMatrix);
	cout<<endl;
	poissonMatrix = solve_gauss_seidel(1000, poissonMatrix, y, x, dx, dy);
	print_solution(poissonMatrix);

	return 0;
}

void print_solution(vvd a){
	for(int i = 0; i < a.size(); i++){
		for(int j = 0; j < a[0].size(); j++){
			cout<<a[i][j]<<" ";
		}
		cout<<endl;
	}
}

// TODO: Implement
void solve_jacobi(){
	return;
}

vvd solve_gauss_seidel(int iters, vvd a, int m, int n, double dx, double dy){
	for(int iter = 0; iter < iters; iter++){
		vvd output(a);
		for(int i = 1; i < m-1; i++){
			for(int j = 1; j < n-1; j++){
				double below = a[i+1][j];
				double here = a[i][j];
				double above = a[i-1][j];
				double right = a[i][j+1];
				double left = a[i][j-1];
				if(iter == 0) cout<<below<<" "<<here<<" "<<above<<" "<<right<<" "<<left<<endl;
				//output[i][j] = (below - 2*here - above) / (dy*dy) + (right - 2*here + left) / (dx*dx);
				output[i][j] = -(1.0/4.0)*(-(i*dx)*pow(dx, 2) + a[i+1][j] + a[i-1][j] + a[i][j-1] + a[i][j+1]);
				if(iter == 0) cout<<output[i][j]<<endl;
			}
		}
		a = output;
	}

	return a;
}

vvd set_boundary_values(vvd a, int m, int n, double dx, double dy, vvd left(vvd, int, double), vvd top(vvd, int, double), vvd right(vvd, int, int, double), vvd bottom(vvd, int, int, double)){
	a = left(a, m, dy);
	a = top(a, n, dx);
	a = right(a, m, n, dy);
	a = bottom(a, m, n, dx);
	return a;
}

vvd initialize_arr(int m, int n){
	vvd output(m, vector<double>(n, 0));
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			output[i][j] = 0;
		}
	}

	return output;
}

// Subtract two vectors from each other
vector<double> subtract_vectors(vector<double> t1, vector<double> t2){
	int t1Size = t1.size();
	int t2Size = t2.size();
	if(t1Size != t2Size){
		return t1;
	}

	vector<double> output(t1Size);
	for(int i = 0; i < t1Size; i++){
		int diff = t2[i] - t1[i];
		output[i] = diff;
	}

	return output;
}

// Calculate the p-norm between two vectors t1 and t2
double calc_norm(vector<double> vec, double p){
	double sum = 0;
	for(int i = 0; i < vec.size(); i++){
		sum += (p == 1) ? abs(vec[i]) : pow(vec[i], p);
	}
	double norm = pow(sum, 1.0/p);
	return norm;
}

vvd verif_bottom(vvd a, int m, int n, double dx){
	for(int i = 0; i < n; i++){
		a[m-1][i] = (i*dx)*exp(1);
	}

	return a;
}

vvd verif_top(vvd a, int n, double dx){
	for(int i = 0; i < n; i++){
		a[0][i] = 2.0*exp(i*dx);
	}

	return a;
}

vvd verif_left(vvd a, int m, double dy){
	for(int i = 0; i < m; i++){
		a[i][0] = 0;
	}

	return a;
}

vvd verif_right(vvd a, int m, int n, double dy){
	for(int i = 0; i < m; i++){
		a[i][n-1] = 2*exp(i*dy);
	}

	return a;
}
