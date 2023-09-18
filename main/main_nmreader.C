// INCLUDES regarding classes
#include"NMReader.h"

#include "TApplication.h"
 
/*
This main provides examples for the NMReader class
*/

int main()  {


  auto app = new TApplication("myapp", nullptr, nullptr);

  string station="OULU";
  string start_day="1";
  string start_month="1";
  string start_year="1950";
  string end_day="1";
  string end_month="9";
  string end_year="2021";
  //string resolution="43200"; // this string controls the periodicity we watch(1440=daily,43200=monthly,525600=yearly and so on, in minutes)
  string resolution="525600"; // this string controls the periodicity we watch(1440=daily,43200=monthly,525600=yearly and so on, in minutes)


  NMReader reader; // simple declaration

  reader.Read_Data(station,start_year,start_month,start_day, end_year,end_month,end_day,resolution);
  //reader.Read_Data("KERG",start_year,start_month,start_day, end_year,end_month,end_day,resolution);
  //reader.Read_Data("KIEL",start_year,start_month,start_day, end_year,end_month,end_day,resolution);
  reader.Read_Data("JUNG",start_year,start_month,start_day, end_year,end_month,end_day,resolution);
  //reader.Read_Data("HRMS",start_year,start_month,start_day, end_year,end_month,end_day,resolution);
  //reader.Read_Data("MOSC",start_year,start_month,start_day, end_year,end_month,end_day,resolution);
  //reader.Read_Data("NEWK",start_year,start_month,start_day, end_year,end_month,end_day,resolution);
  //reader.Read_Data("THUL",start_year,start_month,start_day, end_year,end_month,end_day,resolution);

  // accessing the data vector
  auto data_oulu = reader.GetDataVector(station); // acesso aos dados

  //auto canvas = reader.Draw("OULU");
  //canvas->WaitPrimitive();
  //auto canvas2 = reader.Draw("OULU KERG KIEL JUNG HRMS MOSC NEWK THUL");
  //canvas2->WaitPrimitive();

  //auto AutoCorrelationGraph = reader.GetAutocorrelationGraph(1.,25.,1./4,station);
  //auto AutoCorrelationGraph1 = reader.GetAutocorrelationGraph(0.,1.,1./365.25,station);
  //auto AutoCorrelationGraph2 = reader.GetAutocorrelationGraph(1.,25.,1./4,"JUNG");
  //auto AutoCorrelationGraph3 = reader.GetAutocorrelationGraph(0.,1.,1./365.25,"JUNG");
  auto Mov_avg= reader.GetMovingAverage(5,station);
  auto Mov_avg2=reader.GetMovingAverageGraph(station);


  app->Run();

  return 0;
}