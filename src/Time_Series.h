//class Time_Series
#ifndef __Time_Series__
#define __Time_Series__
#include <vector>
#include <utility>//pair
#include <string>

using namespace std;

class Time_Series{
public:
Time_Series(){;} // default constructor 
~Time_Series(){;} // destructor 

//some Time_Series methods 
virtual vector <pair<double,double>> GetDataVector(string){return Time_Series_vector;} 

protected:

	vector<pair<double,double>> Time_Series_vector; 
};

#endif

