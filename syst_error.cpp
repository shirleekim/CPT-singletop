#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
//#include <ofstream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLine.h>
#include <sstream>

using namespace std;

double calcul(double systttbar, double systmulti, double systsingletop, double systSingleTopbar, double systelectro, double lumi){

  string version = "_mg2";
  double b_mu=100;


  string nameFile1 = "results/histo/hist_tchannel"+version+".root";
  string nameFile2 = "results/histo/hist_tbarchannel"+version+".root";


    TFile* Filet = TFile::Open(nameFile1.c_str());

    TFile* Filetbar = TFile::Open(nameFile2.c_str());

    double range = 100; //*********************** 7 for mg0, mg1, and  4000 for mg2, mg3

    int Nbin = 1e7;
    double pas =  2.0*range/(double) Nbin;

  //signal = modulation *b
    TH1F* datt = (TH1F*)Filet->Get("Asimov");
    TH1F* Sigt = (TH1F*)Filet->Get("Modulation_t");
    TH1F* ttbart = (TH1F*)Filet->Get("ttbar");
    TH1F* electroweakt = (TH1F*)Filet->Get("Electroweak");
    TH1F* multijett = (TH1F*)Filet->Get("Multijet");
    TH1F* SingleTopSMt = (TH1F*)Filet->Get("t_channel");
    TH1F* chiCarre = new TH1F ("Chi Carre", "" , Nbin, -range, range);
    TF1* fithist = new TF1 ("pol2","pol2", -range,range);


    TH1F* dattbar = (TH1F*)Filetbar->Get("Asimov");
    TH1F* Sigtbar = (TH1F*)Filetbar->Get("Modulation_tbar");
    TH1F* ttbartbar = (TH1F*)Filetbar->Get("ttbar"); //comes from hist_tchannel
    TH1F* electroweaktbar = (TH1F*)Filetbar->Get("Electroweak");
    TH1F* SingleTopSMtbar = (TH1F*)Filetbar->Get("tbar_channel");
    TH1F* multijettbar = (TH1F*)Filetbar->Get("Multijet");

    double asimovt = datt->GetBinContent(1);
    double ttbartbkg = ttbart->GetBinContent(1)*systttbar;
    double multijettbkg = multijett->GetBinContent(1)*systmulti;
    double SingleTopSMtbkg = SingleTopSMt->GetBinContent(1)*systsingletop;
    double electroweaktbkg = electroweakt->GetBinContent(1)*systelectro; //multipile par 1

    double asimovtbar = dattbar->GetBinContent(1);
    double ttbartbarbkg = ttbartbar->GetBinContent(1)*systttbar;
    double multijettbarbkg = multijettbar->GetBinContent(1)*systmulti;
    double electroweaktbarbkg = electroweaktbar->GetBinContent(1)*systelectro;
    double SingleTopSMtbarbkg = SingleTopSMtbar->GetBinContent(1)*systSingleTopbar;
    double sigmatott = sqrt(lumi*(ttbartbkg+multijettbkg+electroweaktbkg+SingleTopSMtbkg));
    double sigmatottbar = sqrt(lumi*(ttbartbarbkg+multijettbarbkg+electroweaktbarbkg+SingleTopSMtbarbkg));


    double dividet = 1.0/(sigmatott*sigmatott*24);
    double dividetbar = 1.0/(sigmatottbar*sigmatottbar*24);
    int bin = 1;
 //
    for (double i = 0 ; i < Nbin; i++)
    {
      double b = -range + (double) i * pas;
      double tchannel = 0;
      double tbarchannel = 0;

      for (int p=1; p<26; p++)
      {//X^2
        tchannel += (asimovt-lumi*((b*(Sigt->GetBinContent(p))+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg))*(asimovt-lumi*((b*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg));
        tbarchannel += (asimovtbar-lumi*((b*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg))*(asimovtbar-lumi*((b*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg));
      }

      chiCarre->SetBinContent(bin,tchannel*dividet+tbarchannel*dividetbar);
      bin++;
    }
    chiCarre->Fit(fithist,"R"); //minimisation
    //cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;
  //**** pour bu (stat+syst) ^^^^^^^^si on met les syst != 1.0
return b_mu*fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1);
}

int main ()
{
  //int nProcessus = 6;
  double statsyst1[6];
  double statsyst2[6];
  string process[6] = {"ttbar", "multi jet", "t channel", "tbar channel", "electroweak", "L"};
                    //ttbar, multi jet, t channel, tbar channel, electro, L
  double per[6] = {0.009, 0.017, 0.033, 0.055, 0.05, 0.025};

/*
  cout << "ttabr" << calcul(1.0+per[0], 1.0, 1.0, 1.0, 1.0, 1.0) << " " << calcul(1.0-per[0], 1.0, 1.0, 1.0, 1.0, 1.0) << endl;; //ttbar
  cout << "multi jet" << calcul(1.0, 1.0+per[1], 1.0, 1.0, 1.0, 1.0) << " " << calcul(1.0, 1.0-per[1], 1.0, 1.0, 1.0, 1.0) << endl; //multi jet
  cout << "t channel" << calcul(1.0, 1.0, 1.0+per[2], 1.0+per[3], 1.0, 1.0) << " " << calcul(1.0, 1.0, 1.0-per[2], 1.0-per[3], 1.0, 1.0)<< endl;//single t
  cout << "t bar channel" << calcul(1.0, 1.0, 1.0, 1.0+per[3], 1.0, 1.0) << " " << calcul(1.0, 1.0, 1.0, 1.0-per[3], 1.0, 1.0)<< endl;//single tbar
  cout << "electroweak" <<calcul(1.0, 1.0, 1.0, 1.0, 1.0+per[4], 1.0) << " " << calcul(1.0, 1.0, 1.0, 1.0, 1.0-per[4], 1.0)<< endl;//electro
  cout << "L" << calcul(1.0, 1.0, 1.0, 1.0, 1.0, 1.0+per[5]) << " " << calcul(1.0, 1.0, 1.0, 1.0, 1.0, 1.0-per[5])<< endl;//L
*/

    double stat = calcul(1.0,1.0,1.0,1.0,1.0, 1.0);

    statsyst1[0] = calcul(1.0+per[0], 1.0, 1.0, 1.0, 1.0, 1.0);//ttbar
    statsyst1[1] = calcul(1.0, 1.0+per[1], 1.0, 1.0, 1.0, 1.0); //multi jet
    statsyst1[2] = calcul(1.0, 1.0, 1.0+per[2], 1.0+per[3], 1.0, 1.0);//single t
    statsyst1[3] = calcul(1.0, 1.0, 1.0+per[2], 1.0-per[3], 1.0, 1.0);//single tbar
    statsyst1[4] = calcul(1.0, 1.0, 1.0, 1.0, 1.0+per[4], 1.0);//electro
    statsyst1[5] = calcul(1.0, 1.0, 1.0, 1.0, 1.0, 1.0+per[5]);//L


    statsyst2[0] = calcul(1.0-per[0], 1.0, 1.0, 1.0, 1.0, 1.0); //ttbar
    statsyst2[1] = calcul(1.0, 1.0-per[1], 1.0, 1.0, 1.0, 1.0); //multi jet
    statsyst2[2] = calcul(1.0, 1.0, 1.0-per[2], 1.0-per[3], 1.0, 1.0);//single t
    statsyst2[3] = calcul(1.0, 1.0, 1.0-per[2], 1.0+per[3], 1.0, 1.0);//single tbar
    statsyst2[4] = calcul(1.0, 1.0, 1.0, 1.0, 1.0-per[4], 1.0);//electro
    statsyst2[5] = calcul(1.0, 1.0, 1.0, 1.0, 1.0, 1.0-per[5]);//L

double syst[6];
double max_statsyst[6];
double squaredSum=0;

    for(int i=0; i<6; i++)
    {
        max_statsyst[i] = max(abs(statsyst1[i]),abs(statsyst2[i]));

        if(i==3)
        {
          process[i] = process[i] + " " + process[i-1];
          max_statsyst[i] = max(max_statsyst[i],max_statsyst[i-1]);
        }

        syst[i] = sqrt(abs(max_statsyst[i]*max_statsyst[i] - stat*stat));
        squaredSum += syst[i]*syst[i];

      if(i==2)
      {
        } else
        {
          cout << process[i] << " : " << syst[i] << endl;
        }

}
cout<<"statical error of b_mu when delta Chi^2 at 1 GeV = "<<stat<<endl;

cout << "total systematic error of b mu : " << sqrt(squaredSum);


return 0;

}
