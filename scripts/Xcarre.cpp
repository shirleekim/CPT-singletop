#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLine.h>
#include <sstream>
#include "../include/setting.hpp"

using namespace std;
int main ()
{
  setting();

  string nameFile1 = "results/histo/hist_tchannel"+version+".root";
  string nameFile2 = "results/histo/hist_tbarchannel"+version+".root";

  TFile* Filet = TFile::Open(nameFile1.c_str(),"UPDATE");
  TCanvas* feuille = new TCanvas ("", "", 450, 450);
  TFile* Filetbar = TFile::Open(nameFile2.c_str(),"UPDATE");

  feuille->SetLeftMargin(0.14);
  feuille->SetRightMargin(0.03);
  auto legend = new TLegend(0.65,0.75,0.90,0.85);

  //adapt for each simulation
  //double range = 7; // see setting.cpp 7
  int Nbin = 1e7;
  double pas = 2.0*range/(double) Nbin;

  //change the value to 1 + % of incertitude to get the value of Xcarre with syst_error
  double systttbar = 1.0;
  double systmulti = 1.0;
  double systsingletop = 1.0;
  double systSingleTopbar = 1.0;
  double systelectro = 1.0;
  double lumi = 1.0;
//signal = modulation *b
  TH1F* datt = (TH1F*)Filet->Get("Asimov");
  TH1F* Sigt = (TH1F*)Filet->Get("Modulation_t"); // composant x de f(t) intervient ("Paricle" dans sensibility.cpp)
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

  double asimovt = datt->GetBinContent(1); //somme de toutes les contributions en background + t-channel
  double ttbartbkg = ttbart->GetBinContent(1)*systttbar;
  double multijettbkg = multijett->GetBinContent(1)*systmulti;
  double SingleTopSMtbkg = SingleTopSMt->GetBinContent(1)*systsingletop;
  double electroweaktbkg = electroweakt->GetBinContent(1)*systelectro; //multipile par 1

  double asimovtbar = dattbar->GetBinContent(1);
  double ttbartbarbkg = ttbartbar->GetBinContent(1)*systttbar;
  double multijettbarbkg = multijettbar->GetBinContent(1)*systmulti;
  double electroweaktbarbkg = electroweaktbar->GetBinContent(1)*systelectro;
  double SingleTopSMtbarbkg = SingleTopSMtbar->GetBinContent(1)*systSingleTopbar;
  //cout << "ttbartbarbkg : " << ttbartbarbkg << endl; //inserted
  double sigmatott = sqrt(lumi*(ttbartbkg+multijettbkg+electroweaktbkg+SingleTopSMtbkg));
  double sigmatottbar = sqrt(lumi*(ttbartbarbkg+multijettbarbkg+electroweaktbarbkg+SingleTopSMtbarbkg));

  chiCarre->GetXaxis()->SetTitle("b_{#mu} (GeV)");
  chiCarre->GetYaxis()->SetTitle("#delta#chi^{2}");
  chiCarre->GetXaxis()->CenterTitle(kTRUE);
  chiCarre->GetYaxis()->CenterTitle(kTRUE);
  chiCarre->SetStats(kFALSE);

  double dividet = 1.0/(sigmatott*sigmatott*24);
  double dividetbar = 1.0/(sigmatottbar*sigmatottbar*24);
  int bin = 1;


  for (double i = 0 ; i < Nbin; i++)
  {
    double b = -range + (double) i * pas;
    double tchannel = 0;
    double tbarchannel = 0;

    for (int p=1; p<25; p++)
    {//calcu; de X^2
      tchannel += (asimovt-lumi*((b*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg))*(asimovt-lumi*((b*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg));
      tbarchannel += (asimovtbar-lumi*((b*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg))*(asimovtbar-lumi*((b*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg));
    }
    chiCarre->SetBinContent(bin,(tchannel*dividet+tbarchannel*dividetbar));
    bin++;
  }
  string nameXroot = "XCarre"+version+".root";
  TFile* Xroot = new TFile(nameXroot.c_str(),"RECREATE");
  Xroot->cd();
  chiCarre->Fit(fithist,"R");
  chiCarre->Write("chiCarre");

  //**** It gives delta bu (stat+syst) when we put syst != 1.0
  //See setting.cpp for bmu
  cout<<"statistical error of b_mu when delta Chi^2 at 1 GeV = "<<bmu*fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;//modifier

  chiCarre->SetLineColor(kRed);
  chiCarre->Draw();
  legend->AddEntry(chiCarre, "#chi_{tot}^{2} = #chi_{t}^{2} +#chi_{#bar{t}}^{2}", "l");


  /*chiCarre = new TH1F ("Chi Carre", "" , range*2/0.0001, -range, range);
  lumi=1.025;
  dividet = 1.0/(sigmatott*sigmatott*24);
  dividetbar = 1.0/(sigmatottbar*sigmatottbar*24);

  double buff=0;
  bin = 1;

  for (double i = -range ; i<range-0.0001; i+=0.0001)
  {

    double tchannel = 0;
    double tbarchannel = 0;
    for (int p=1; p<26; p++)
    {
      tchannel += (asimovt-lumi*((i*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg))*(asimovt-lumi*((i*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg));
      tbarchannel += (asimovtbar-lumi*((i*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg))*(asimovtbar-lumi*((i*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg));
    }
    chiCarre->SetBinContent(bin,tchannel*dividet+tbarchannel*dividetbar);
    //cout<<chiCarre->GetBinContent(bin)<<endl;
    bin++;
  }
  //cout<<chiCarre->GetBinContent(chiCarre->GetMinimumBin())<<endl;
  //cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;

  double min = chiCarre->GetBinContent(chiCarre->GetMinimumBin());
  cout<<"min="<<min<<endl;
  for (int bin = 1 ; bin<=140001; bin++)
  {
    buff = chiCarre->GetBinContent(bin);
    chiCarre->SetBinContent(bin,buff-min);

  }

  chiCarre->SetLineColor(kBlue);
  chiCarre->Draw("SAME");
  chiCarre->Write("XCarre+");
  legend->AddEntry(chiCarre, "Luminosity = 154 fb^{-1}", "l");
  chiCarre = new TH1F ("Chi Carre", "" , range*2/0.0001, -range, range);
  buff=0;

  lumi = 0.975;
  dividet = 1.0/(sigmatott*sigmatott*24);
  dividetbar = 1.0/(sigmatottbar*sigmatottbar*24);
  bin=1;
  for (double i = -range ; i<range-0.0001; i+=0.0001)
  {
    double tchannel = 0;
    double tbarchannel = 0;

    for (int p=1; p<26; p++)
    {
      tchannel += (asimovt-lumi*((i*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg))*(asimovt-lumi*((i*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg));
      tbarchannel += (asimovtbar-lumi*((i*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg))*(asimovtbar-lumi*((i*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg));
    }
    chiCarre->SetBinContent(bin,tchannel*dividet+tbarchannel*dividetbar);
    bin++;
  }
  cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;
  min = chiCarre->GetBinContent(chiCarre->GetMinimumBin());

  for (int bin = 1 ; bin<=140001; bin++)
  {
    buff = chiCarre->GetBinContent(bin);
    chiCarre->SetBinContent(bin,buff-min);
  }

  chiCarre->SetLineColor(kRed);
  chiCarre->Draw("SAME");
  chiCarre->Write("XCarre-");
  legend->AddEntry(chiCarre, "Luminosity = 146 fb^{-1}", "l");*/
  legend->Draw();
  TLine* deltaChi = new TLine (-range, 1, range, 1);
  TLine* deltaChiPlus = new TLine (-abs(fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)), 0, -abs(fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)) , 1);
  TLine* deltaChiMoins = new TLine (abs(fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)), 0, abs(fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)), 1);
  //deltaChi->SetLineColor(kBlue);
  deltaChi->Draw();
  //legend->AddEntry(deltaChi, "#delta #chi^{2}", "l");
  deltaChiPlus->Draw();
  deltaChiMoins->Draw();
  feuille->SaveAs(("results/graph/XCarre"+version+".pdf").c_str());
  Xroot->Close();
  return 0;
}
