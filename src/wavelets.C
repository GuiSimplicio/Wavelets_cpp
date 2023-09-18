// INCLUDEs
#include "wavelets.h"
#include <cmath> //used for M_PI for example
#include<complex> //used for imag->complex number i 
#include<iostream>// std::cout, std::endl

using namespace std:: complex_literals;

//transforms time series from time space to Fourier space(frequency's)
/*
This function progress through the steps necessary for a wavelet analysis of a given time_series (parameter given).
It returns a Scalogram.
This whole class is based upon: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
  
*/
TGraph2D* wavelets::Get_Wavelet_Transform(vector <pair<double,double>> time_series, double delta_j, const char* GraphTitle){

	int N = time_series.size(); //N is the number of points in the time series

	//***************************************************************************************
	//Step 1: Find the Fourier transform of the (possibly padded) time series

	vector <complex<double>> Discrete_Fourier_Transform =  DFT(time_series); // Equation 3

	
	//Prints to check if DFT is done right
	/*
	for(int i=0;i<N;i++){
		cout << "x[n] is ... " << time_series[i].second << endl;
	}

	for_each( Discrete_Fourier_Transform.begin(), Discrete_Fourier_Transform.end(), [ ](complex<double> m){ cout<<m<<" "<< endl;;}); // check if DFT is done right
	*/
	// Done right !!! Checked with an online calculator

	//***************************************************************************************

	//Step 2: Choose a wavelet function and a set of scales to analyze.

	//Wavelet Function choosen is Morlet
	//To choose an appropriate set of scales one must first discover the delta_t(time interval) of the time series
	//For efficiency purposes it is convenient to calculate the variance(sigma^2 in the article) of the time series in this loop

	double delta_t = time_series[1].first-time_series[0].first;
	cout << "Delta t is: " << delta_t << endl;
	//parameters to calculate variance
	double mean_time_series=0;
	double mean2_time_series=0;
	double variance_time_series=0;

	for(int i=0;i<N-1;i++){

		if(time_series[i+1].first-time_series[i].first > delta_t*1.1){ //here we multiply delta_t for 1.1 to give a tolerance of 10%, one can be more strict in this calculation
			cout << "Warning! Time series is not continuous! Exiting... because of point "<< time_series[i+1].first << endl;
			exit(1);
		}
		mean_time_series += (time_series[i].second);
		mean2_time_series += (time_series[i].second)*(time_series[i].second) ;
	}

	mean2_time_series /= N; 
	mean_time_series /= N;
	variance_time_series=mean2_time_series-mean_time_series*mean_time_series;


	// scale parameters:
	// s[j]=s0*pow(2.,j*delta_j), j=0,1,2,...,J -> Defined through a loop 
	// J= (log2((N*delta_t)/(2*s0))/delta_j) + 1 ;
	
	double w0=6.;//default value so Morlet wave function has a zero mean  
	double scale_period_conversion= 4.*M_PI/(w0+sqrt(2.+w0*w0)); // w0=6
 	//delta_j=default value(0.005) or value called on main; // this parameter can't be higher than 0.5 for the Morlet wavelet
 					  // it should be ajustated to provide a smooth picture of wavelet wavelet_power_spectrum
 	double s0=2*delta_t; //should be ~ Fourier Period 
	int J=(log2((N*delta_t)/(2*s0))/delta_j) + 1 ; // Equation 10 (modified)
	cout << "J is ... " << J << endl;
	double s[J+1]; //scale 

	
	//***************************************************************************************
	//Step 3: For each scale, construct the normalized wavelet function using Equation (6)
	//Step 4: Find the wavelet transform at that scale using Equation (4)
	//Step 5: Determine the cone of influence and the Fourier wavelength at that scale

	vector<vector<complex<double>>> wavelet_transform_vec;// vector of vector od the wavelet transform-> each scale has N transforms

	//Creating the COI with TGraph	
	Cone_of_Influence = new TGraph(2*J+2);
	double e_folding_time[J+1]; // τ_s -> used to calculate Cone of Influence (Table 1) 
	//Loop that creates each scale
	for(int j=0;j<J+1;j++){

                if (j%100 == 0)
                        cout << 100.*j/(J+1) << "%" << endl;

		s[j]=s0*pow(2.,j*delta_j);//Equation 9
		e_folding_time[j]=sqrt(2)*s[j]; // formula on Table 1 -> Calculation for COI

		Cone_of_Influence->SetPoint(j,time_series[0].first+e_folding_time[j],s[j]); 
		Cone_of_Influence->SetPoint(2*J-j+1,time_series[N-1].first-e_folding_time[j],s[j]); 


		//Now Step 3-> Each Scale => Normalized wavelet function
		// Normalized Wavelet Function depends of omega (Table 1), so we must start a double loop that takes care of both step 3 and 4 simultaneously

		wavelet_transform_vec.push_back(Wavelet_Transform(time_series,s[j],delta_t,Discrete_Fourier_Transform));  // Equations 4,5 and 6 inside Wavelet_Transform
		//tambem posso fazer o push_back do vector que tenho na classe!   
	}

	//***************************************************************************************
	// Step 6: After repeating steps 3–5 for all scales,  remove any padding and contour  plot  the  wavelet  power spectrum.

	TGraph2D *graph = new TGraph2D();  

	double power;

	for(int j=0;j<J+1;j++){

		for(int n=0;n<N;n++){
			power=0;
			power=pow(abs(wavelet_transform_vec[j][n]),2); 
			graph->SetPoint(j*N+n,time_series[n].first,s[j],power/variance_time_series); //very important to divide by variance to normalize the power spectrum! 
			graph->SetTitle(GraphTitle);
		}
	}

	//***************************************************************************************
	//Step 7: Assume  a  background  Fourier  power  spectrum(e.g., white or red noise) at each scale,
	// then use the chi-squared distribution to find the 95% confidence(5% significance) contour.

	//Creating a background Fourier Power Spectrum (Equation 16 for red noise)
	//discrete Fourier power spectrum assuming Red-noise:
	int aux_size=N/2 + 1;
	double alpha, alpha_1, alpha_2; 
	double Fourier_red_noise_spectrum[aux_size];

	//Not sure if it is done well !!! *******************************
	
	alpha_1=GetAutocorrelationCoef(time_series,0.,1.,1./365.25); //lag-1 auto correlation  of the time series
	alpha_2=GetAutocorrelationCoef(time_series,0.,2.,1./365.25);//lag-2 auto correlation of the time series
	
	//***************************************************************
	alpha= (alpha_1+sqrt(alpha_2))/2;
	
	for(int k=0;k<aux_size;k++){
		Fourier_red_noise_spectrum[k]=(1-alpha*alpha)/(1+alpha*alpha-2*alpha*cos(2*M_PI*k/N));//equation 16

	}
	
	//now multiply Fourier_red_noise_spectrum[k] by 5.99146/2 and check if it is bigger than power ? the dimensions dont really make sense here!
	
	// ver função das autocorrelações do NMReader e adapatar aqui! 

	//To determine the 95% confidence level (significant at 5%), one multiplies
	//the background spectrum (16) by the 95th percentile value for chi-squared
	//95% percentile value for chi_squared with 2DOF is : 5.99146

	// After finding an appropriate background spectrum and choosing a particular
	//confidence for χ2 such as 95%, one can then calculate (18) at each scale and
	// construct 95% confidence contour lines.

	/*
	Portanto o que eu acho é: basta usar (16) e depois (18) para obter os níveis de significancia!
	Como construir um gráfico daqui? 
	Se for um TGraph une todos os pontos, correto?
	Deixar para depois
	*/
	



	return graph; 
}

//This method returns the Discrete Fourier Transform of the time series 
// Receives as parameter: the time series itself
vector<complex<double>> wavelets::DFT(const vector <pair<double,double>> &time_series){

	int N = time_series.size(); //N is the number of points in the time series
	vector<complex<double>> DFT_Time_Series;
	complex<double> x[N];

	//Equation 3  

	for(int k=0;k<N;k++){

		x[k]=(0.,0.);

		for(int n=0;n<N;n++){

			x[k]+= time_series[n].second * exp(-2*M_PI*1i*(double)k*(double)n/(double)N);
		}

		DFT_Time_Series.push_back(x[k] /= N); 
	}

	return DFT_Time_Series;
}


//This method returns the Inverse Fourier Transform of the wave,ie, the wavelet transform
// Receives as parameter: the time series, the scale, delta_t as its required to calculate the normalized wave function,
// the DFT of the time series (required in formula 4)
vector<complex<double>> wavelets::Wavelet_Transform(const vector <pair<double,double>> &time_series, double s, double delta_t, const vector <complex<double>> &Discrete_Fourier_Transform){
	
	int N = time_series.size(); //N is the number of points in the time series
	double omega; // angular frequency
	complex<double> wave; // will store the normalized wave
	complex<double> W[N]; // auxilary for the wavelet transform
	vector<complex<double>> wavelet_transform;
	
	for(int n=0;n<N;n++){

		W[n]=(0.,0.);

		for(int k=0;k<N;k++) {

			//Equation 5 to determine omega
			if(k<= N/2){

				omega=(2* M_PI*k)/(N*delta_t); 

			}else{
					omega= -(2* M_PI*k)/(N*delta_t);
				}

			wave=Normalized_Wave(omega,s,delta_t);//Equation 6 to determine Normalized wave
			W[n]+= Discrete_Fourier_Transform[k]* wave * exp(1i*omega* (double)n * delta_t);  // Equation 4
			//cout << "W is ... " << W[n] << endl;
			//cout << "Wave is ... " << wave << endl;
			//cout << "Discrete_Fourier_Transform is ... " << Discrete_Fourier_Transform[k] << endl;
			//cout << "Expression is ... " << exp(1i*omega* (double)n * delta_t) << endl;

		}

		wavelet_transform.push_back(W[n]);
	}

	
	return wavelet_transform;

}



//returns the complex conjugate of the normalized wave function
complex<double> wavelets::Normalized_Wave(double omega,double s,double delta_t){

	complex<double> wave;
	complex<double> wave_complex_conjugate;
	double w0=6.0; //default Morlet value

	//wave function (Morlet)

	//Heaviside function  
	if(omega > 0){

		wave= pow(M_PI,-0.25)*exp(-pow((s*omega-w0),2) /2); // Morlet Wave equation (table 1)
	}else {

		wave=(0.0,0.0); 
	}

	//normalized wave function

	wave= sqrt(2*M_PI*s/delta_t)* wave; // Equation 6

	//complex conjugated of normalized wave function

	wave_complex_conjugate=conj(wave); // We require the complex conjugate of the normalized wave function as 
									   // We can see on equation 4


	return wave_complex_conjugate;
}


// returns the COI
TGraph* wavelets::GetConeOfInfluence(){

	Cone_of_Influence->SetMarkerStyle(20);
	Cone_of_Influence->SetMarkerSize(0.8);
	Cone_of_Influence->SetMarkerColor(kRed+2);
	Cone_of_Influence->SetLineColor(kRed+2);
	Cone_of_Influence->SetLineWidth(1.5);


	return Cone_of_Influence;
} 

// Returns the correlation coeficient R 
// initial_time and final_time are related to the interval swept by the time shift(k), ie k sweeps [initial_time,final_time] with a given time step
//k is the lag of the time series! 
double wavelets::GetAutocorrelationCoef(const vector <pair<double,double>> &time_series, double initial_time, double final_time, double time_step) {


  int size=time_series.size();
  double init_year=time_series[0].first;
  double end_year=time_series[size-1].first;


  // initial year of data is  time_series[0].first
  // final year of data is  time_series[size].first

  double k;

  // X [initial_time, final_time-k]
  // X' [initial_time+k, final_time]

  //auxiliar variables to calculate R
  vector<pair<double, double>> X, XX;
  double sumX, sumXX;
  double SUM;
  double sdX, sdXX;

  double R; // correlation coeficient
  int n = 0;

  //k sweeps [initial_time,final_time] with a given time step 
  for (k = initial_time; k <= final_time; k += time_step) {
    //Clearing data
    X.clear();
    XX.clear();
    sumX = 0;
    sumXX = 0;
    SUM = 0;
    sdX = 0;
    sdXX = 0;

    // data vector for a time shift k

    for_each(time_series.begin(), time_series.end(), [&initial_time, &final_time, &k, &X, &XX, &init_year, &end_year](pair<double, double> i){


  		if ((i.first >= init_year) && (i.first <= end_year-k)){
   	 		X.push_back(i);
  		}
  		if ((i.first >= init_year+k) && (i.first <= end_year)){
    		XX.push_back(i);
  		}
     });

    if (X.size() != XX.size()) {
      cout << "Vectors don't have same size in GetAutocorrelationCoef" << endl;
      exit(1);
    }


    // Average

    for (int i = 0; i < X.size(); i++) {
      sumX += X[i].second;
      sumXX += XX[i].second;
    }
    sumX = sumX / X.size();
    sumXX = sumXX / XX.size();

    // Numerator of R

    for (int i = 0; i < X.size(); i++) {
      SUM += ((X[i].second-sumX)*(XX[i].second-sumXX));
    }

    // Standart deviation

    for (int i = 0; i < X.size(); i++) {
      sdX += (X[i].second-sumX)*(X[i].second-sumX);
      sdXX += (XX[i].second-sumXX)*(XX[i].second-sumXX);
    }


    // R ...
    R = SUM / sqrt(sdX * sdXX);

  }

  return R;
}



/*
//returns time series padded with zeroes in the end because Fourier Transform assumes data is cyclic
vector <pair<double,double>> wavelets::Pad_with_zeroes(int size,vector <pair<double,double>> time_series){

// Pad the end of the time series with zeroes before doing the wavelet transform
//the time series is padded with sufficient zeroes to bring the total length N up to the next-higher power of two

	int aux=0;
	while(pow(2,aux))<size{
		aux++;
	} 

	int new_size=pow(2,aux+1);

	vector <pair<double,double>> time_series_with_0(new_size);

	for(int i=size;i<new_size;i++){
		time_series_with_0[i].push_back(make_pair(time_series[i-1].first+delta_t,0));
	}

	return time_series_with_0;

}
*/
