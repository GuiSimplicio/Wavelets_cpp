// INCLUDES regarding classes
#include"NMReader.h"
#include "wavelets.h"

//Include's regarding root
#include "TApplication.h"
 #include "TStyle.h"

  int main(){


  auto app = new TApplication("myapp", nullptr, nullptr);

  vector <pair <double,double>> artificial_data;
  int N_points=1000;
  double delta_t_sampling = 0.1;
  double Signal_Period1=10;
  double Signal_Period2=2;
  double w01 = (2*M_PI)/ Signal_Period1 ;
  double w02 = (2*M_PI)/ Signal_Period2 ;

  double time,wave; 

  TGraph* signal= new TGraph(N_points);

  for(int i=0; i<N_points; i++){

  	time=delta_t_sampling*i;

  	//Defining some waves:
  	//wave=sin(w01*i*delta_t_sampling+0.3);//Simple wave
  	//wave=sin(w01*i*delta_t_sampling+0.7)+sin(w02*i*delta_t_sampling+0.3); // double wave
  	wave=sin(delta_t_sampling*i*(w01+(w01-w02)*delta_t_sampling*i/10)); // chirp
  	//wave=sin(w01*i*delta_t_sampling+0.3); //some other wave??

    artificial_data.push_back(make_pair(time,wave));
    signal->SetPoint(i,time,wave);

  }


  //Check if signal is well built:

  TCanvas* c1= new TCanvas ("signal", "signal");
  signal->SetTitle("Artificial Signal");
  signal->GetXaxis()->SetTitle("t(s)");
  signal->GetYaxis()->SetTitle("x(t)");
  signal->Draw();
  c1->SaveAs("wave_artsig4.png");
  c1->WaitPrimitive();

  //**********************************


  wavelets WV; // creates an objetc
  
  //TGraph2D* graph = WV.Get_Wavelet_Transform(artificial_data);//receives an graph from Wavelet_transform whom needs to receive a time_series as an argument
  TGraph2D* graph = WV.Get_Wavelet_Transform(artificial_data,0.025,"Wavelet Transform;Time(s);Period(s);Power");//one can call the function giving it a non-default value for delta_j (for more information read wavelets.h)
  																			   //and a non-default name for both the axis and the Graph Title.
  																			   //Note that you can define either the delta_j or the graph title separetely, ie, you can call Get_Wavelet_Transform as :
  																			   //Get_Wavelet_Transform(artificial_data,0.025) or Get_Wavelet_Transform(artificial_data,"Wavelet Transform;Period;Frequency;Power") 
  																			   //if you wish to just set one of those two default values  

  auto COI= WV.GetConeOfInfluence(); 

  TCanvas* c2= new TCanvas ("Cwavelet2", "wavelet2");
  TFile *F = new TFile("contour_plot.root","RECREATE");
  gStyle->SetPalette(112);
  graph->Draw("colz");
  COI->Draw("C");
  //c2->SetLogy();
  c2->WaitPrimitive();
  graph->Write();
  COI->Write();
  F->Close();
  app->Run();

  return 0;

  } 

