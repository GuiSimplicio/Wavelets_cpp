// INCLUDES regarding classes
#include"NMReader.h"
#include "wavelets.h"

//Include's regarding root
#include "TApplication.h"
#include "TStyle.h"
	
/*
This main provides examples for the class wavelets, in this case uses NMReader class to get a time-series
*/
  int main(){


  auto app = new TApplication("myapp", nullptr, nullptr);

  string station="OULU";
  string start_day="1";
  string start_month="1";
  string start_year="1950";
  string end_day="1";
  string end_month="9";
  string end_year="2021";
  string resolution="525600"; // this string controls the periodicity we watch(1440=daily,43200=monthly,525600=yearly and so on, in minutes)
  //string resolution="43200"; // this string controls the periodicity we watch(1440=daily,43200=monthly,525600=yearly and so on, in minutes)


  NMReader reader; // simple declaration

  reader.Read_Data(station,start_year,start_month,start_day, end_year,end_month,end_day,resolution);

  // accessing the data vector
   auto data_oulu = reader.GetDataVector(station);
   //auto canvas = reader.Draw("OULU");

    //canvas->WaitPrimitive();


  wavelets WV; // creates an objetc
  
/*

To draw the Fourier Transform in order to try to obtain th period

  auto DFT_vec=WV.DFT(data_oulu);
  TGraph* FT = new TGraph();

  for(int i=0;i<data_oulu.size();i++){
  	FT->SetPoint(i,i,DFT_vec[i].real());
  }

  TCanvas* c1 = new TCanvas ("Fourier Series", "FS");
  FT->GetXaxis()->SetTitle("Frequency(Hz)");
  FT->GetYaxis()->SetTitle("Fourier Transform");


  FT->Draw();

*/
  
  //TGraph2D* graph = WV.Get_Wavelet_Transform(data_oulu);//receives and graph from Wavelet_transform whom needs to receive a time_series as an argument
  TGraph2D* graph = WV.Get_Wavelet_Transform(data_oulu,0.025,"Wavelet Transform;Fractionary Year;Period(Fractionary Year);Power");//one can call the function giving it a non-default value for delta_j (for more information read wavelets.h)
  																			   //and a non-default name for both the axis and the Graph Title.
  																			   //Note that you can define either the delta_j or the graph title separetely, ie, you can call Get_Wavelet_Transform as :
  																			   //Get_Wavelet_Transform(data_oulu,0.025) or Get_Wavelet_Transform(data_oulu,"Wavelet Transform;Period;Frequency;Power") 
  																			   //if you wish to just set one of those two default values    

  auto COI= WV.GetConeOfInfluence(); 


  TCanvas* c2= new TCanvas ("Cwavelet2", "wavelet2");
  TFile *F = new TFile("contour_plot.root","RECREATE");
  gStyle->SetPalette(112);
  graph->Draw("colz");
  COI->Draw("C");
  c2->SetLogy();
  c2->WaitPrimitive();
  graph->Write();
  COI->Write();
  F->Close();
  //c2->SaveAs("OULU_timeframe.png");
  
  app->Run();

  return 0;

  } 
