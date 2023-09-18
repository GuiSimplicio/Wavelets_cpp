// INCLUDEs
#include "NMReader.h"
#include "TMultiGraph.h"
#include "TLegend.h"

//destructor
NMReader::~NMReader(){
  map_vpair.clear();
}



void NMReader::Read_Data(string station, string start_year, string start_month, string start_day, string end_year, string end_month, string end_day, string resolution){

  string url;

  url="http://nest.nmdb.eu/draw_graph.php?wget=1&stations[]="
     +station+"&tabchoice=revori&dtype=corr_for_efficiency&tresolution="+resolution+"&yunits=0&date_choice=bydate&start_day="
     +start_day+"&start_month="+start_month+"&start_year="+start_year+"&start_hour=0&start_min=0&end_day="
     +end_day+"&end_month="+end_month+"&end_year="+end_year+"&end_hour=0&end_min=0&output=ascii&asciifrac=1";

  //url explanation...
  // http://nest.nmdb.eu/draw_graph.php?wget=1 - wget needs to be =1
  // &stations[]=KERG - station
  // &tabchoice=revori - Revised Original
  // &dtype=corr_for_efficiency - corrections
  // &tresolution=0 - best time resolution(?)
  // &yunits=0 - counts/s
  // &date_choice=bydate -
  // &start_day=17 -
  // &start_month=7 -
  // &start_year=2021 -
  // &start_hour=0 -
  // &start_min=0 -
  // &end_day=22 -
  // &end_month=7 -
  // &end_year=2021 -
  // &end_hour=23 -
  // &end_min=59 -
  // &output=ascii - .txt
  // &asciifrac=1 - fractionary date

  string cmd;
  cmd="wget -O - -o /dev/null \""+url+"\"";

  FILE *fp = popen(cmd.c_str(),"r");

  if(fp==NULL){

    cout<<"Error... File of station "<<station <<" doesn't exist... exiting..." << endl;
    pclose(fp);
    exit(1);

  }else{
    while(!feof(fp)){

      double Neutron_Monitor,frac_date;
      int size=128;
      char aux[size];
      char aux1;

      if(fgets(aux,size,fp)!=NULL){
        if(sscanf(aux,"%lf%c%lf",&frac_date,&aux1,&Neutron_Monitor)==3){ // only reads data that matters
          //stores vector on map with the key of the station
          // map_vpair[station];
          map_vpair[station].push_back(make_pair(frac_date,Neutron_Monitor));
        }
      }
    }
  }

  cout << "Reading data from station " << station << ". Loaded " << map_vpair[station].size() << " entries." << endl;

  pclose(fp);
}


//Prints vector of data of certain station:
void NMReader::PrintData(string station){
    for(int i=0;i<map_vpair[station].size();i++){

    cout << "Station: " << station << endl;

    cout<< setprecision(7) <<map_vpair[station][i].first << " " << map_vpair[station][i].second <<endl;

  }
}

//Get Data vector
vector <pair<double,double>> NMReader::GetDataVector(string st){
  if(! map_vpair.count(st)){
    cout << "Error !!! Station: "<<st << " not found ... closing program" << endl;
    exit(1);
  }

  return map_vpair[st];
}

TCanvas* NMReader::Draw(string stations) {
  auto mgraph = new TMultiGraph();
  auto legend = new TLegend();


  stringstream ss(stations);
  string station;
  const char* station_char;

  while(ss >> station) {
    mgraph->Add(GetDataGraph(station),"PL");
    station_char=station.c_str();
    legend->AddEntry(GetDataGraph(station),station_char,"lp");
  }

  auto canvas = new TCanvas("nmreader_canvas", "NMReader canvas");
  mgraph->Draw("AP PLC PMC");
  legend->Draw();
  mgraph->SetTitle("Neutron Monitor");
  mgraph->GetYaxis()->SetTitle("Counts/s");
  mgraph->GetXaxis()->SetTitle("Date(Year)");
  canvas->Update();
  return canvas;
}

TGraph* NMReader::GetDataGraph(string st){

  if(! map_vpair.count(st)){
    printf("[%s] Station %s not found.\n", __PRETTY_FUNCTION__, st.c_str());
    return nullptr;
  }

  TGraph* Gdata = new TGraph;

  for (int i = 0; i < map_vpair[st].size(); i++)
  {
      Gdata->SetPoint(i, map_vpair[st][i].first, map_vpair[st][i].second);
  }

  int legend_index;

  if(legend_index==0){
  	Gdata->SetMarkerStyle(20);
  	Gdata->SetMarkerSize(0.4);
  	Gdata->SetMarkerColor(kBlue+1);
	}else if(legend_index==1){
		Gdata->SetMarkerStyle(20);
  		Gdata->SetMarkerSize(0.4);
  		Gdata->SetMarkerColor(kRed+1);
	}else if(legend_index==2){
		Gdata->SetMarkerStyle(20);
  		Gdata->SetMarkerSize(0.4);
  		Gdata->SetMarkerColor(kGreen+1);
	}else{
		Gdata->SetMarkerStyle(20);
  Gdata->SetMarkerSize(0.4);
  Gdata->SetMarkerColor(kYellow+1);
	}
  

  legend_index++;

  return Gdata;
}

map <string,vector<pair<double,double>>> NMReader::GetMovingAverage(int avg_window, string station){
  int size=map_vpair[station].size();
  int c = (avg_window+2)/2;
  double sum=0;

  map_moving_average[station];

  for(int i=0;i<c;i++){
    map_moving_average[station].push_back(make_pair(map_vpair[station][i].first,map_vpair[station][i].second));
  }

  if(avg_window % 2 == 0){ // check if avg_window is an even integer

    for (int i = c; i < size-c+2; i++) {
    sum = 0;
    for (int j = i-c; j < i+c-1; j++) {
      sum += map_vpair[station][j].second;
    }
    map_moving_average[station].push_back(make_pair(map_vpair[station][i].first, sum/avg_window));
  }

  for (int i = size-c+2; i < size; i++){
    map_moving_average[station].push_back(make_pair(map_vpair[station][i].first, map_vpair[station][i].second));
  }

  }else {

    for (int i = c; i < size-c+3; i++) {
      sum = 0;
      for (int j = i-c; j < i+c-2; j++) {
        sum += map_vpair[station][j].second;
      }
      map_moving_average[station].push_back(make_pair(map_vpair[station][i].first, sum/avg_window));
    }
    for (int i = size-c+3; i < size; i++){
      map_moving_average[station].push_back(make_pair(map_vpair[station][i].first, map_vpair[station][i].second));
    }
  }

  return map_moving_average;
}

TGraph* NMReader::GetMovingAverageGraph(string station){


  TGraph* Gdata = new TGraph;
  TCanvas *c3 = new TCanvas("c3", "Multigraph" ,800, 600);

  for (int i = 0; i < map_moving_average[station].size(); i++)
  {
      Gdata->SetPoint(i, map_moving_average[station][i].first, map_moving_average[station][i].second);
  }

  Gdata->SetMarkerStyle(20);
  Gdata->SetMarkerSize(0.4);

  Gdata->SetTitle("Neutron Monitor moving average");
  Gdata->GetYaxis()->SetTitle("NM");
  Gdata->GetXaxis()->SetTitle("Date(Year)");
  Gdata->SetMarkerColor(kYellow);
  Gdata->Draw();
  //c3->SaveAs("11_day.pdf");

  //This last two lines are comented to speed up program
  //They are needed if I want to see the graph being drawn

  return Gdata;

}

// Returns a pointer to a TGraph of the correlation coeficient R vs  time shift k
// initial_time and final_time are related to the interval swept by the time shift(k), ie k sweeps [initial_time,final_time] with a given time step 
TGraph* NMReader::GetAutocorrelationGraph(double initial_time, double final_time, double time_step, string station) {



  TGraph* Gdata = new TGraph();
  int size=map_vpair[station].size();
  double init_year=map_vpair[station][0].first;
  double end_year=map_vpair[station][size-1].first;


  // initial year of data is  map_vpair[station][0].first
  // final year of data is  map_vpair[station][size].first

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

    for_each(map_vpair[station].begin(), map_vpair[station].end(), [&initial_time, &final_time, &k, &X, &XX, &init_year, &end_year](pair<double, double> i){


  		if ((i.first >= init_year) && (i.first <= end_year-k)){
   	 		X.push_back(i);
  		}
  		if ((i.first >= init_year+k) && (i.first <= end_year)){
    		XX.push_back(i);
  		}
     });

    if (X.size() != XX.size()) {
      cout << "Vectors don't have same size in GetAutocorrelationGraph" << endl;
      exit(1);
    }

    /*
    for(int n=0;n<size;n++){
    	cout << X[n].first << " " << X[n].second << endl; 
    	cout << XX[n].first << " " << XX[n].second << endl; 

    }

	*/

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
    // Insert point in Graph
    Gdata->SetPoint(n, k, R);
    //cout << "Sum is " << SUM << endl;
    //cout << "sqrt(sdX * sdXX) is " << sqrt(sdX * sdXX) << endl;
    n++;

  }

  //Graph related ...


    Gdata->SetMarkerStyle(22);
    Gdata->SetMarkerSize(0.4);
    Gdata->SetMarkerColor(kRed+2);
  	Gdata->SetMarkerColor(0);
    Gdata->SetNameTitle(station.c_str(), station.c_str());
    Gdata->SetLineWidth(5);
    Gdata->SetLineColor(kBlue+4);



    //Gata->GetXaxis->SetTitle("");
    //Gata->GetYaxis->SetTitle("");
  TCanvas *c4 = new TCanvas("c4", "Auto Correlation Graph");


  Gdata->Draw();
  //c4->SaveAs("Gcorr_detailed.pdf");






  return Gdata;
}
