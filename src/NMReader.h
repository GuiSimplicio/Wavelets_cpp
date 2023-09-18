//classe NMReader
#ifndef __NMReader__
#define __NMReader__

//Includes
#include<iostream>// std::cout, std::endl
#include <string>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <vector>
#include<algorithm>
#include<iomanip>
#include<utility>
#include <fstream> // read from file
#include <numeric>
#include <stdlib.h>
#include <stdio.h> //popen
#include <map>
#include<utility>//pair
#include <sstream>

//Includes regarding root
#include"TObject.h"
#include"TFile.h"
#include"TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

using namespace std;

#include "Time_Series.h"// inherits from this class  

class NMReader: public Time_Series {
public:
NMReader(){;} //default constructor
~NMReader(); // destructor

void Read_Data(string station, string start_year, string start_month, string start_day, string end_year, string end_month, string end_day,string resolution=0); //reads data from web and stores it in map_vpair


vector <pair<double,double>> GetDataVector(string); // returns vector with (fractionary date, Neutron Monitor)
TGraph* GetDataGraph(string); // returns plot Neutron Monitor vs fractionary date
void PrintData(string station); // Given station, prints its vector of data
TCanvas* Draw(string stations); // draw plots for all the stations. separate station names with a space (" ")

//some methods to check periodicty of NM:

map <string,vector<pair<double,double>>> GetMovingAverage(int, string); // moving average of a int number of days to average over and a given station(string)
TGraph* GetMovingAverageGraph(string);//returns plot of the moving average of NM vs fractionary date         
TGraph* GetAutocorrelationGraph(double initial_time, double final_time, double time_step, string station); // Returns a pointer to a TGraph of the correlation coeficient R vs  time shift k

protected:
	map <string,vector <pair<double,double>>> map_vpair; //map that contains a key(station) to its data vector(fractionary date, Neutron Monitor)
	map <string,vector <pair<double,double>>> map_moving_average; //map that contains a key(station) to its moving data vector(fractionary date, Neutron Monitor)

};

#endif
