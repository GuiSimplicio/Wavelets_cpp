//class wavelets
#ifndef __wavelets__
#define __wavelets__
#include<complex> //used for imag->complex number i 
#include <vector>
#include <utility>//pair
#include "TGraph2D.h"
#include "TGraph.h"

using namespace std;

class wavelets{
public:
wavelets(){;} // default constructor 
//need to add another constructor with Time_Series as a parameter!
~wavelets(){;} // destructor 

//some wavelets methods using Morlet function
TGraph2D* Get_Wavelet_Transform(vector <pair<double,double>> time_series,double delta_j=0.005,const char* GraphTitle="Wavelet Transform;Period;Frequency;Power");//transforms time series from time space to Fourier space(frequency's), and returns its graph
																		//delta_j here receives a default value, the lower it is the the finer resolution we get! 
																		//this parameter can't be higher than 0.5 for the Morlet wavelet
 					  													// it should be ajustated to provide a smooth picture of wavelet wavelet_power_spectrum
																		//Graph Title is merely the title of the graph and its axis names aswell
vector<complex<double>> DFT(const vector <pair<double,double>> &time_series); // returns the Discrete Fourier Transform of the time series
vector<complex<double>> Wavelet_Transform(const vector <pair<double,double>> &time_series, double s, double delta_t, const vector <complex<double>> &Discrete_Fourier_Transform); //returns the wavelet transform
complex<double> Normalized_Wave(double omega,double s,double step); //returns the complex conjugate of the normalized wave function
												// arguments are omega,scale and time step(dt). (omega=angular frequency)
TGraph* GetConeOfInfluence(); // returns the COI
double GetAutocorrelationCoef(const vector <pair<double,double>> &time_series, double initial_time, double final_time, double time_step); //returns the autocorrelation coefficient of a time series with a certain lag


protected:

	vector<vector<pair<double,complex<double>>>> V_Wavelet_transform; // vector of vector of a pair(set of scale and the wavelet transform)
	TGraph* Cone_of_Influence; // TGraph that contains the Cone of Influence(COI)
};

#endif

