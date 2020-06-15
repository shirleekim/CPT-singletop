#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include "../include/Analyzetw.hpp"

#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TH2F.h>
#include <TLegend.h>
#include <sstream>
#include "../include/setting.hpp"

using namespace std;

int main (){

/* included in setting.cpp
  string version = "";
  string simulation_t = "1Mt13TeV";//MG0 = 1Mt13TeV ,
  string simulation_tbar = "1Mtbar13TeV"; // MG0 = 1Mtbar13TeV

  int nbevent = 1e6; // mg0 : 1e6, // mg 1,2,3: 1e5
*/

  setting();

  TString diagram;
  cout<<"We consider an energy in the center-of-mass of 13 TeV for both t-channel and tw-channel"<<endl;
  do {
    cout<<"Choose which channel you wish to study within the 4 diagrams available :"<<endl
    <<"t - t-channel through b and q diffusion"<<endl
    <<"tw - tW-channel"<<endl
    <<"tbar - tbar-channel through antib and q diffusion"<<endl
    <<"tbarw - tbarW-channel"<<endl;
    cin>>diagram;}
  while (diagram !="t" && diagram !="tw" && diagram !="tbar" && diagram !="tbarw" );






  // Some needed variables
  int choice;

  int nbtqprecut = 0;
  int nbtqbarprecut = 0;
  int nbtq = 0;
  int nbtqbar =0;
  int nbtw = 0;
  double SM = 0;
  double ratio =0;
  TLorentzVector SME (0. ,0. ,0. ,0.);
  TLorentzVector Wtq (0., 0., 0., 0.);
  TLorentzVector Wtqbar (0., 0., 0., 0.);
  TLorentzVector Wtw (0., 0., 0., 0.);
  TLorentzVector Wt (0., 0., 0., 0.);
  TLorentzVector WtAv(0., 0., 0., 0.);



  //Choose the channel and diagram we want to study.

  //Extract the corresponding rootfile
  if (diagram == "t" || diagram == "tbar")
  {
    root_simu* go;
    if (diagram == "t")
    go = new root_simu(simulation_t+".root");
    if (diagram == "tbar")
    go = new root_simu(simulation_tbar+".root");
    go->Loop();
    int Events = go->fChain->GetEntriesFast();

    //Declare the variable to sum the values of the matrix elements for each iteration
    for (int i=0; i<Events; i++){

      if (go->nature!=0) {
        MatrixSingleTop *M = new MatrixSingleTop(go->pMother1[i], go->b[i], go->pt[i], go->p3[i], go->pMother1_PID[i], go->p3_PID[i], go->b_PID[i], go->nature[i], go->pTmu[i], go->etamu[i], go->pTelec[i], go->etaElec[i], go->pTb[i], go->pTnon_b[i], go->eta_b[i], go->eta_non_b[i]);
        choice = go->nature[i];
          //Treat t-channel in the case of tq and tbarq
          if (choice == 1) {
            nbtqprecut++;
            if (abs(go->eta_non_b[i])>=4.7 || go->pTnon_b[i]<=40) continue; //eliminates quark non b
            if (abs(go->eta_b[i])>=2.4 ||  go->pTb[i]<=40) continue; // eliminates quark
            if ((abs(go->etaElec[i])<1.479 && go->pTelec[i]>35) || (abs(go->etamu[i])<2.4 && go->pTmu[i]>26)){ //select electron or muon

              if (diagram == "t")
              {
                SME = M->getMbq_mu();
              }
              if (diagram == "tbar")
              {
                SME = -1.0*M->getMbq_mu();
              }

              SM  = M->getMbq();
              ratio = 1/SM;
              Wtq+=SME*ratio;//Sum each calculated change ratio at each iteration

              //cout << " SME " << SME[2] << endl;
              nbtq++;         //Count the number of time we calculate a charge ratio to calculate the mean value later
            }

        }

          //Treat t-channel in the case of tqbar and tbarqbar
          if (choice == 2 ){
            nbtqbarprecut++;
            if (abs(go->eta_non_b[i])>=4.7 || go->pTnon_b[i]<=40) continue; //same as above
            if (abs(go->eta_b[i])>=2.4 ||  go->pTb[i]<=40) continue; // same as above
            if ((abs(go->etaElec[i])<1.479 && go->pTelec[i]>35) || (abs(go->etamu[i])<2.4 && go->pTmu[i]>26)){ // same as above
              if (diagram == "t")
              {
                SME = M->getMbqbar_mu();
              }
              if (diagram == "tbar")
              {
                SME = -1.0*M->getMbqbar_mu();
              }
              SM  = M->getMbqbar();
              ratio = 1/SM;
              Wtqbar+=SME*ratio;
              nbtqbar++;
              }

          }
      }

      //delete M;

    }
  }

/*
  if (diagram == "tw" || diagram == "tbarw")
  {
    Analyzetw* go;
    if (diagram == "tw")
    go = new Analyzetw("1Mtw13TeV.root");
    if (diagram == "tbarw")
    go = new Analyzetw("tbarwchan_13TeV_fixed.root");
    go->Loop();
    int Events = go->fChain->GetEntriesFast();

    //Declare the variable to sum the values of the matrix elements for each iteration
    for (int i=0; i<Events; i++){


      if (go->nature!=0) {
        MatrixSingleTop *M = new MatrixSingleTop(go->pMother1[i], go->b[i], go->pt[i], go->p3[i], go->pMother1_PID[i], go->p3_PID[i], go->b_PID[i], go->nature[i], go->pTmu[i], go->etamu[i], go->pTelec[i], go->etaElec[i]);
        choice = go->nature[i];

          //Ttreat tW-channel in the case of t and tbar
          if(choice == 3){
            if (abs(go->etaElec[i])<2.4 && go->pTelec[i]>20){
              if (diagram == "tbarw")
              {
                SME = -1.0*M->getMbg_mu();
              }
              if (diagram == "tw")
              {
                SME = M->getMbg_mu();
              }
              SM  = M->getMbg();
              ratio = 1/SM;
              Wtw+=SME*ratio;
              nbtw++;
            }

            else if (abs(go->etamu[i])<2.4 && go->pTmu[i]>20){
              if (diagram == "tbarw")
              {
                SME = -1.0*M->getMbg_mu();
              }
              if (diagram == "tw")
              {
                SME = M->getMbg_mu();
              }
              SM  = M->getMbg();
              ratio = 1/SM;
              Wtw+=SME*ratio;
              nbtw++;
            }

            else if (abs(go->eta_b[i])<2.4 && go->pTb[i]>40){
              if (diagram == "t")
              {
                SME = M->getMbq_mu();
              }
              if (diagram == "tbar")
              {
                SME = -1.0*M->getMbq_mu();
              }
             SM  = M->getMbq();
             ratio = 1/SM;
             Wtq+=SME*ratio;
             nbtq++;
           }
           else if (abs(go->eta_non_b[i])<4.7 && go->pTnon_b[i]>40){
             if (diagram == "t")
             {
               SME = M->getMbq_mu();
             }
             if (diagram == "tbar")
             {
               SME = -1.0*M->getMbq_mu();
             }
            SM  = M->getMbq();
            ratio = 1/SM;
            Wtq+=SME*ratio;
            nbtq++;
          }
          }

      delete M;

    }
  }
  */
//cout<<1<<endl;
//Calculate the ratio considering we have different diagram for a single generation
//When we simule tq and so tqbar, we have different cross-sect for tq and tqbar. The following expression weights each calculated matrix element.
double ratiotq = 1.0/nbtq;
double ratiotqbar=1.0/nbtqbar;
double ratiotw = 1.0/nbtw;
double pourctq = 1.0*nbtqprecut/nbevent;
double pourctqbar = 1.0*nbtqbarprecut/nbevent;
//If Wtw != 0 then Wtq and Wtqbar should be 0


if (diagram == "t"){
  cout << "nbtq" << nbtq << endl;
  cout << "nbtqbar" << nbtqbar << endl;
  cout << "nbevent" << nbevent << endl;


  WtAv = pourctq*Wtq*ratiotq + pourctqbar*Wtqbar*ratiotqbar;

  //WtAv = Wtq*ratiotq;
  //WtAv = Wtqbar*ratiotqbar;
  ofstream tchan("results/matrix_elements/t-channel"+version+".txt");
  for (int i=0; i<4; i++) {
    cout << "Wtq [" << i << "] " <<Wtq[i] << endl;
    cout << "Wtqbar [" << i << "] " << Wtqbar[i]<<endl;

    tchan<<WtAv[i]<<endl;
  }

  tchan.close();
}

if (diagram == "tw" ){
  WtAv =Wtw*ratiotw;
  ofstream twchan("results/matrix_elements/tw-channel"+version+".txt");
  for (int i=0; i<4; i++) {twchan<<WtAv[i]<<endl;}
  twchan.close();
}
if (diagram == "tbar" ){
  WtAv =pourctq*Wtq*ratiotq +pourctqbar*Wtqbar*ratiotqbar;

  ofstream twchan("results/matrix_elements/tbar-channel"+version+".txt");
  for (int i=0; i<4; i++) {twchan<<WtAv[i]<<endl;}
  twchan.close();
}
if (diagram == "tbarw" ){
  WtAv =Wtw*ratiotw;

  ofstream tbarwchan("results/matrix_elements/tbarw-channel"+version+".txt");
  for (int i=0; i<4; i++) {tbarwchan<<WtAv[i]<<endl;}
  tbarwchan.close();
}

//________________Macro creating the f(t) graph using _______________//

//string b_mu = "1"; included in setting.cpp

double b_muvalue = atof(b_mu.c_str());
string nameFile = "results/modulation/Modulation_Temporelle_13TeV_b"+b_mu+""+version+".root";
TLorentzVector b_mu0 (b_muvalue, 0., 0., 0.);//x
TLorentzVector b_mu1 (0., 0., 0., 0.); //y  //test when b_muvalue ==0
TLorentzVector b_mu2 (0., 0., 0., 0.); //z //test when b_muvalue ==0
TLorentzVector b_mu3 (0., 0., 0., b_muvalue); //t

TFile* File = new TFile (nameFile.c_str(), "UPDATE");



auto legend = new TLegend(0.70,0.75,0.95,0.95);
TGraph* graphx = new TGraph (48);
TGraph* graphy = new TGraph (48);
TGraph* graphz = new TGraph (48);
TGraph* grapht = new TGraph (48);
TGraph* graphSM = new TGraph(48);

double rad = M_PI/180.;
double theta = 101.2790 * rad;
double lambda = 46.309 * rad;
double w = 2*M_PI/24;
double ft = 0;
double ft_x, ft_y, ft_z, ft_t;
double sl_cth = sin(lambda)*cos(theta); // eta 1 dans le rapport de Sane
double sth = sin(theta); // eta 2
double cl_cth = cos(lambda)*cos(theta); // eta 3
double cth = cos(theta);
double sl_sth = sin(lambda)*sin(theta);
double sl = sin(lambda);
double cl_sth = cos(lambda)*sin(theta);
double cl = cos(lambda);


//Compute the temporal modulation of the top cross section

double min,max;
if (diagram == "t" || diagram == "tbar")
{min = -0.01; max = 0.01;}

if (diagram == "tw" || diagram == "tbarw")
{min = -0.01; max =0.01;}

int i = 0;
for (double t=0; t<=24 ; t+=0.5)
{// WeightAverage
  //delta P0/P = WtAv[3], delta P1/P = WtAv[0], delta P2/P = WtAv[1],  delta P3/P = WtAv[2]
  // b_mu0[0] = b^1, b_mu1[1] = b^2, b_mu2[2] = b^3, b_mu3[3] = b^0
  ft_x = b_mu0[0]*((cos(w*t)*sl_sth-sin(w*t)*cth)*WtAv[0]-(cos(w*t)*cl)*WtAv[1]+(cos(w*t)*sl_cth-sin(w*t)*sth)*WtAv[2]);
  graphx->SetPoint(i, t, ft);
  //cout << "matrix elements : " << WtAv[0] << " " << WtAv[1] << " " << WtAv[2] << endl;
  //cout << " terms : " << (cos(w*t)*sl_sth-sin(w*t)*cth)*WtAv[0] << "    " << (cos(w*t)*cl)*WtAv[1] << "    " << (cos(w*t)*sl_cth-sin(w*t)*sth)*WtAv[2] << endl;
  //cout<< "f(t) x : " << ft << endl;
  if (ft>max) {do {max += 0.01;} while (ft>max);};
  if (ft<min) {do {min -=0.01;} while (ft<min);};

  ft_y = b_mu1[1]*((sin(w*t)*sl_sth+cos(w*t)*cth)*WtAv[0]-(sin(w*t)*cl)*WtAv[1]+(sin(w*t)*sl_cth+cos(w*t)*sth)*WtAv[2]);
  graphy->SetPoint(i, t, ft);
  if (ft>max) {do {max += 0.01;} while (ft>max);};
  if (ft<min) {do {min -=0.01;} while (ft<min);};

  ft_z = -b_mu2[2]*((cl_sth)*WtAv[0]+(sl)*WtAv[2]+(cl_cth)*WtAv[2]);
  graphz->SetPoint(i, t, ft);
  if (ft>max) {do {max += 0.01;} while (ft>max);};
  if (ft<min) {do {min -=0.01;} while (ft<min);};

  ft_t = b_mu3[3]*WtAv[3];
 grapht->SetPoint(i, t, ft);
  if (ft>max) {do {max += 0.01;} while (ft>max);};
  if (ft<min) {do {min -=0.01;} while (ft<min);};

  ft = ft_t+ft_x+ft_y+ft_z;
  graphx->SetPoint(i, t, ft);


  i++;

}


//Draw the graph of temporal modulation and save the data in rootfile

string namex = "b_{x} = "+b_mu+" GeV";
string namey = "b_{y} = "+b_mu+" GeV";
string namez = "b_{z} = "+b_mu+" GeV";
string namet = "b_{t} = "+b_mu+" GeV";

graphx->SetLineWidth(3);
graphy->SetLineWidth(3);
graphz->SetLineWidth(3);
grapht->SetLineWidth(3);

graphx->SetLineColor(kOrange);
graphy->SetLineColor(kRed);
graphz->SetLineColor(kBlue);
grapht->SetLineColor(8);

legend->AddEntry(graphx,namex.c_str(), "l");
legend->AddEntry(graphy,namey.c_str(), "l");
legend->AddEntry(graphz,namez.c_str(), "l");
legend->AddEntry(grapht,namet.c_str(), "l");


File->cd();
graphx->Write(diagram+"X");
graphy->Write(diagram+"Y");
graphz->Write(diagram+"Z");
grapht->Write(diagram+"T");

File->Close();

/*for (double t=0; t<=24 ; t+=0.5)
{
  ft = 0;
  graphSM->SetPoint(i, t, ft);
  i++;
}*/

//graphSM->SetLineColor(kBlack);
//graphSM->SetLineWidth(3);
//graphSM->Write(diagram+"SM");
//legend->AddEntry(graphSM,"SM", "l");

TCanvas* modul = new TCanvas("Modulation", "", 600 , 600);
modul->SetLeftMargin(0.14);
modul->SetRightMargin(0.03);
TH2F* axe = new TH2F (" ", " ",49 ,0 ,24 ,49 ,-0.01 ,0.01 );

axe->SetStats(kFALSE);
axe->GetYaxis()->SetTitleOffset(2.0);
axe->GetXaxis()->SetTitle("t (h)");
axe->GetYaxis()->SetTitle("f(t) = #sigma_{SME}/#sigma_{SM} - 1 ");
//axe->GetYaxis()->SetTitle("Charge asymmetry");
axe->GetXaxis()->CenterTitle(kTRUE);
axe->GetYaxis()->CenterTitle(kTRUE);
axe->Draw("SAME");
graphx->Draw("SAME");
graphy->Draw("SAME");
graphz->Draw("SAME");
grapht->Draw("SAME");
graphSM->Draw("SAME");
legend->Draw();


modul->SaveAs("results/graph/Modulation_"+diagram+version+".pdf");


return 0;
}
