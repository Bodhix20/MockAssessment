#include <iostream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

//function prototypes
double f(double x);
double fDerivative(double x);
double forwardSkim(double x, double dx);
double backwardSkim(double x, double dx);
double centralSkim(double x, double dx);

double norm1(vector<double>& v);
double norm2(vector<double>& v);
double norm3(vector<double>& v);

int main(){

	double dx = 0.1 ;
	int steps = 1/dx +1;	

	vector<double> v1;
	vector<double> v2;
	vector<double> v3;

	cout << "x" << "\t"  << "f'(x)" << "\t" << "Fwd" << "\t" << "Error" << "\t" << "Bkwd" << "\t" << "Error" << "\t" << "Cent" << "\t" << "Error" <<"\n" ;
	cout << "************************************************************" << "\n";

	for(int i = 0 ; i<steps ; i++){

		double fwd = forwardSkim(i*dx,dx);
		double bkwd = backwardSkim(i*dx,dx);
		double cent = centralSkim(i*dx,dx);
		double fder = fDerivative(i*dx);

		double fwdError = (fwd - fder) / fder;
		double bkwdError = (bkwd - fder) / fder;
		double centError = (cent - fder) / fder;

		cout << i*dx << "\t" << fder << "\t" << fwd << "\t" << fwdError << "\t" << bkwd << "\t" << bkwdError << "\t" << cent << "\t" << centError << "\n";

		v1.push_back(fwdError);
		v2.push_back(bkwdError);
		v3.push_back(centError);
	} 

  // Norm calculation, comment this part if you dont want it to go in the csv report
	cout << "Error Norm 1 (sum of absolutes) normalized: " << "\n";
	cout << "Forward : " << norm1(v1)/steps << "\n" << "Backward : " << norm1(v2)/steps << "\n" << "Central : " << norm1(v3)/steps << "\n";
	cout << "Error Norm 2 (sqrt of sum of absolutes sqrd) normalized : " << "\n";
	cout << "Forward : " << norm2(v1)/steps << "\n" << "Backward : " << norm2(v2)/steps << "\n" << "Central : " << norm2(v3)/steps << "\n";
	cout << "Error Norm 3 (max value) : " << "\n";
	cout << "Forward : " << norm3(v1) << "\n" << "Backward : " << norm3(v2) << "\n" << "Central : " << norm3(v3) << "\n";

	return 0;
}

double f(double x){
	return 3*pow(x,3)+2*x+1;
	//return 9*pow(x,2)+2;

}

//analytical deriviative of the function
double fDerivative(double x){
	//return 3*pow(x,3)+2*x+1;
	return 9*pow(x,2)+2;
}

double forwardSkim(double x, double dx){	
	return (f(x+dx)-f(x))/dx;
}

double backwardSkim (double x, double dx){
	return (f(x)-f(x-dx))/dx;
}

double centralSkim(double x, double dx){
	return (f(x+dx)-f(x-dx))/(2*dx);
}

//norm1 : sum of absolute values
double norm1(vector<double>& v)
{
	double sum = 0.0;
	for (int i = 0; i < v.size(); i++) sum += fabs(v[i]);
	return sum;
}

//norm 2 : euclidian norm
double norm2(vector<double>& v)
{
	double sum = 0.0;
	for (int i = 0; i < v.size(); i++) sum += v[i] * v[i];
	return sqrt(sum);
}

//norm 3 : uniform norm
double norm3(vector<double>& v)
{
	double max = 0.0;
	for (int i = 0; i < v.size(); i++)
		if (fabs(v[i]) > max) max = v[i];
	return max;
}



