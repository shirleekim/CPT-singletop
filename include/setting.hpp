#include <iostream>
using namespace std;

//In case of any values changed here, do 'make clean' then 'make' on terminal
//Mad Graph generated simulation type(0,1,2,3)
  int mg = 0;

  string version ; // neccessary for b_mu2020, sensibility, Xcarre, syst_error
  string simulation_t ; // neccessary for b_mu2020
  string simulation_tbar; // neccessary for b_mu2020
  string b_mu ; // neccessary for b_mu2020, sensibility
  double bmu; // neccessary for Xcarre, syst_error
  int nbevent;

  double range ; // neccessary for Xcarre, syst_error
  double stat_error; // neccessary for sensibility
  int modul_range;
  int times; // for signal amplification (for modulation graph) //necessary for sensibility

void setting (){
do{
  cout << "Enter the simulation type generated by MadGraph (possible values : 0,1,2,3)" << endl;
  cin >> mg;
}while ((mg !=0) and (mg !=1) and (mg !=2) and (mg!=3));


if (mg==0){

  version = "_mg0";
  simulation_t = "1Mt13TeV";//MG0 = 1Mt13TeV ,
  simulation_tbar = "1Mtbar13TeV"; // MG0 = 1Mtbar13TeV
  nbevent = 1e6; // mg0 : 1e6, // mg 1,2,3: 1e5
  b_mu = "1";
  bmu = 1;
  range = 6;
  stat_error = 4.61;
  modul_range = 1e5;
  times = 10;

}

if (mg==1){

  version = "_mg1";
  simulation_t = "MG1";//MG0 = 1Mt13TeV ,
  simulation_tbar = "MG1"; // MG0 = 1Mtbar13TeV
  nbevent = 1e5; // mg0 : 1e6, // mg 1,2,3: 1e5
  b_mu = "1";
  bmu = 1;
  range = 10;
  stat_error = 4.67;
  modul_range = 1e5;
  times = 10;

}

if (mg==2){

  version = "_mg2";
  simulation_t = "MG2";//MG0 = 1Mt13TeV ,
  simulation_tbar = "MG2"; // MG0 = 1Mtbar13TeV
  nbevent = 1e5;
  b_mu = "1";
  bmu = 1;
  range = 10;
  stat_error = 3.55;
  modul_range = 12e4;
  times = 1;


}

if (mg==3){

  version = "_mg3";
  simulation_t = "MG3";//MG0 = 1Mt13TeV ,
  simulation_tbar = "MG3"; // MG0 = 1Mtbar13TeV
  nbevent = 1e5;
  b_mu = "1";
  bmu=1;
  range = 7;
  stat_error = 3.55;
  modul_range = 12e4;
  times = 1;

}
}
