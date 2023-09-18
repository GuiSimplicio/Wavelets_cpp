// INCLUDES regarding classes
#include"NMReader.h"

#include "TApplication.h"
 
/*
This main provides examples for the NMReader class
*/

int main()  {


  auto app = new TApplication("myapp", nullptr, nullptr);

  string station="OULU";
  string start_day="15";
  string start_month="7";
  string start_year="2021";
  string end_day="17";
  string end_month="7";
  string end_year="2021";

  NMReader reader; // simple declaration

  reader.Read_Data(station,start_year,start_month,start_day, end_year,end_month,end_day);
  reader.Read_Data("AATA",start_year,start_month,start_day, end_year,end_month,end_day);

  // accessing the data vector
   auto data_oulu = reader.GetDataVector(station); // acesso aos dados

  auto canvas = reader.Draw("OULU AATA");
  canvas->WaitPrimitive();


  app->Run();

  return 0;
}