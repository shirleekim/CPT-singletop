#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>
#include <sstream>


using namespace std;


int main (){
//////////////change for each simulation /////////////////////////
  int nbevent = 100000; // Number of events simulation takes

  string version = "_mg2";
  string b_mu = "100";

  double b_mu_chiSquare = 5.3; // change when temporal modulation histogram is needed. (it doesn't affect Xcarre.cpp or syst_error.cpp)
//////////////////////////////////////////////////////////////////

  string nameFile = "results/modulation/Modulation_Temporelle_13TeV_b"+b_mu+version+".root";

  TFile* File = TFile::Open(nameFile.c_str());



  auto legend = new TLegend(0.70,0.75,0.98,0.95);

  TCanvas* Signalt = new TCanvas ("", "", 450, 450);
  TGraph* Parti = (TGraph*)File->Get("tX");
  TGraph* AntiParti = (TGraph*)File->Get("tbarX");
  double* Particle = Parti->GetY();
  double* AntiParticle = AntiParti->GetY();



 TH1F* h1sigt = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammat = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbart = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultit = new TH1F ("","", 24, 0, 24);
  TH1F* h5ftsignt = new TH1F ("","",24 ,0,24);
  TH1F* h6errort = new TH1F ("","",24,0,24);

  h1sigt->SetStats(kFALSE);
  h2bggammat->SetStats(kFALSE);
  h3bgttbart->SetStats(kFALSE);
  h4bgMultit->SetStats(kFALSE);
  h5ftsignt->SetStats(kFALSE);
  h6errort->SetStats(kFALSE);
  h1sigt->GetXaxis()->SetTitle("t (h)");
  h1sigt->GetYaxis()->SetTitle("Events");
  h1sigt->GetXaxis()->CenterTitle(kTRUE);
  h1sigt->GetYaxis()->CenterTitle(kTRUE);
  //h1sigt->GetYaxis()->SetRangeUser(0,80000); //pourquoi commentaire????


  //Value of backgrounds like ttbar,EW, Multijet
  double backgroundt_1=8446; // Multijet
  double backgroundt_2=18298; //EW
  double backgroundt_3=39392; //ttbar
  double backgroundt_4=9048; // signal t

  double sigmatott = sqrt(backgroundt_1+backgroundt_2+backgroundt_3+backgroundt_4);

  int k=0;
//
  for (double i = 1 ; i<=24 ; i++)
  {

    h1sigt->SetBinContent(i,backgroundt_1+backgroundt_2+backgroundt_3+backgroundt_4);
    h3bgttbart->SetBinContent(i,backgroundt_1+backgroundt_2+backgroundt_3);
    h4bgMultit->SetBinContent(i,backgroundt_1);
    h2bggammat->SetBinContent(i,backgroundt_1+backgroundt_2);
    h5ftsignt->SetBinContent(i, (backgroundt_1+backgroundt_2+backgroundt_3+backgroundt_4+backgroundt_4*Particle[k]*10*b_mu_chiSquare)); //signal x 10
    //h1sigt->SetBinError(i,sigmatott);
    h6errort->SetBinError(i, sigmatott);
    h6errort->SetBinContent(i, backgroundt_1+backgroundt_2+backgroundt_3+backgroundt_4);
    k+=2;
  }

  h1sigt->SetFillColor(kRed+1);
  h2bggammat->SetFillColor(kBlue+1);
  h3bgttbart->SetFillColor(kGreen+2);
  h4bgMultit->SetFillColor(kOrange-3);
  h5ftsignt->SetLineColor(kViolet-3);
  legend->AddEntry(h1sigt,"t channel", "f");
  legend->AddEntry(h2bggammat, "W/Z/#gamma + jets", "f");
  legend->AddEntry(h3bgttbart,"t#bar{t}/tW", "f");
  legend->AddEntry(h4bgMultit,"Multijet", "f");
  legend->AddEntry(h5ftsignt,"10x expected signal b_{x} = 5.37 GeV","l");
  h1sigt->GetYaxis()->SetRangeUser(0,100000); //commentaire? ********************************



  h1sigt->Draw("SAME");
  h3bgttbart->Draw("SAME");
  h2bggammat->Draw("SAME");
  h4bgMultit->Draw("SAME");
  h6errort->Draw("SAME");
  h5ftsignt->Draw("SAME"); //5.3
  legend->Draw("SAME");


  //Tcanvas
  Signalt->SaveAs(("results/graph/Histo_t"+version+".pdf").c_str());


  string name = "results/histo/hist_tchannel"+version+".root";
//______________Save t_channel to rootFile____________//
  TFile* tchannel = new TFile (name.c_str(),"RECREATE");

  TH1F* h1sigthist = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammathist = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbarthist = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultithist = new TH1F ("","", 24, 0, 24);
  TH1F* h5ftsignthist = new TH1F ("","",24 ,0,24);
  TH1F* asimovt = new TH1F ("","",24,0,24);

  double error1 = sqrt(backgroundt_4);
  double error2 = sqrt (backgroundt_2);
  double error3 = sqrt (backgroundt_3);
  double error4 = sqrt (backgroundt_1);
  double error0 = sqrt(backgroundt_2+backgroundt_1+backgroundt_3+backgroundt_4);

  tchannel->cd();

  k=0;
  for (double i = 1 ; i<=24 ; i++)
  { //TH1F
    asimovt->SetBinContent(i,backgroundt_4+backgroundt_2+backgroundt_3+backgroundt_1);
    h1sigthist->SetBinContent(i,backgroundt_4);
    h2bggammathist->SetBinContent(i,backgroundt_2);
    h3bgttbarthist->SetBinContent(i,backgroundt_3);
    h4bgMultithist->SetBinContent(i,backgroundt_1);
    h5ftsignthist->SetBinContent(i,Particle[k]);
    h1sigthist->SetBinError(i,error1);
    h3bgttbarthist->SetBinError(i,error3);
    h4bgMultithist->SetBinError(i,error4);
    h2bggammathist->SetBinError(i,error2);
    asimovt->SetBinError(i,error0); // 0 : sum of every bacgrounds
    k+=2;
  }
  asimovt->Write("Asimov");
  h1sigthist->Write("t_channel");
  h3bgttbarthist->Write("ttbar");
  h4bgMultithist->Write("Multijet");
  h2bggammathist->Write("Electroweak");
  h5ftsignthist->Write("Modulation_t");

  tchannel->Close();


  string name2 = "results/histo/hist_tbarchannel"+version+".root";
//_______________tbar channel sensibility plot_________________//
  TFile* tbarchannel = new TFile (name2.c_str(),"RECREATE");

  TCanvas* Signaltbar = new TCanvas("","", 450, 450);
  TH1F* h1sigtbar = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammatbar = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbartbar = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultitbar = new TH1F ("","", 24 , 0, 24);
  TH1F* h5ftsigntbar = new TH1F ("","", 24,0 ,24);
  TH1F* h6errortbar = new TH1F ("","",24,0,24);
  auto legendtbar = new TLegend(0.70,0.73,0.98,0.98);


  h1sigtbar->SetStats(kFALSE);
  h2bggammatbar->SetStats(kFALSE);
  h3bgttbartbar->SetStats(kFALSE);
  h4bgMultitbar->SetStats(kFALSE);
  h5ftsigntbar->SetStats(kFALSE);
  h6errortbar->SetStats(kFALSE);

  h1sigtbar->GetXaxis()->SetTitle("t (h)");
  h1sigtbar->GetYaxis()->SetTitle("Events");
  h1sigtbar->GetXaxis()->CenterTitle(kTRUE);
  h1sigtbar->GetYaxis()->CenterTitle(kTRUE);

  double backgroundtw_1 = 8404;
  double backgroundtw_2 = 16233;
  double backgroundtw_3 = 39618;
  double backgroundtw_4 = 5739; //signal de t

  double sigmatottbar = sqrt(backgroundtw_1+backgroundtw_2+backgroundtw_3+backgroundtw_4);
  k=0;


  for (double i = 1 ; i<=24 ; i++)
  {
  //  h5tfsingtbar->Scale(10);//*********************************
    h1sigtbar->SetBinContent(i,backgroundtw_1+backgroundtw_2+backgroundtw_3+backgroundtw_4); //asimov
    h3bgttbartbar->SetBinContent(i,backgroundtw_1+backgroundtw_2+backgroundtw_3);
    h4bgMultitbar->SetBinContent(i,backgroundtw_1);
    h2bggammatbar->SetBinContent(i,backgroundtw_1+backgroundtw_2);
    //antiparticle[k] : f(t) centree en 0 au
    h5ftsigntbar->SetBinContent(i,(backgroundtw_1+backgroundtw_2+backgroundtw_3+backgroundtw_4)+backgroundtw_4*AntiParticle[k]*5.37*10);//commentaire
    h6errortbar->SetBinContent(i,backgroundtw_1+backgroundtw_2+backgroundtw_3+backgroundtw_4);
    h6errortbar->SetBinError(i,sigmatottbar);
    k+=2;

  }
//  h5ftsigntbar->Scale(10);//*********************************
  h1sigtbar->SetFillColor(kRed+1);
  h2bggammatbar->SetFillColor(kBlue+1); //electroweak
  h3bgttbartbar->SetFillColor(kGreen+2);
  h4bgMultitbar->SetFillColor(kOrange-3); //multijet
  h5ftsigntbar->SetLineColor(kViolet-3); //commentaire

  legendtbar->AddEntry(h1sigt,"#bar{t} channel", "f");
  legendtbar->AddEntry(h3bgttbart,"t#bar{t}/tW", "f");
  legendtbar->AddEntry(h2bggammat, "W/Z/#gamma* + jets", "f");
  legendtbar->AddEntry(h4bgMultit,"Multijet", "f");
  legendtbar->AddEntry(h5ftsignt,"10x expected signal b_{x} = 5.37 GeV","l"); //commentaire?

  h1sigtbar->GetYaxis()->SetRangeUser(0,nbevent);//*********************************
  h1sigtbar->SetFillColor(kRed+1);
  h2bggammatbar->SetFillColor(kBlue+1);
  h3bgttbartbar->SetFillColor(kGreen+2);
  h4bgMultitbar->SetFillColor(kOrange-3);
  h5ftsignt->SetLineColor(kViolet-3); //commentaire?

  h1sigtbar->Draw("SAME");

  h3bgttbartbar->Draw("SAME");

  h2bggammatbar->Draw("SAME");
  h4bgMultitbar->Draw("SAME");
  h6errortbar->Draw("SAME");
  h5ftsigntbar->Draw("SAME");
  legendtbar->Draw();
  Signaltbar->SaveAs(("results/graph/Histo_tbar"+version+".pdf").c_str());



//______________Save tbar_channel to rootFile_____________//
TH1F* h1sigthistbar = new TH1F ( "", "", 24 , 0, 24);
TH1F* h2bggammathistbar = new TH1F ( "", "", 24 , 0, 24);
TH1F* h3bgttbarthistbar = new TH1F ("","",24 , 0, 24);
TH1F* h4bgMultithistbar = new TH1F ("","", 24, 0, 24);
TH1F* h5ftsignthistbar = new TH1F ("","",24 ,0,24);
TH1F* asimovtbar = new TH1F ("","",24 ,0,24);

  tbarchannel->cd();

  k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    asimovtbar->SetBinContent(i,backgroundtw_4+backgroundtw_3+backgroundtw_1+backgroundtw_2);//asimov
    h1sigthistbar->SetBinContent(i,backgroundtw_4);
    h3bgttbarthistbar->SetBinContent(i,backgroundtw_3);
    h4bgMultithistbar->SetBinContent(i,backgroundtw_1);
    h2bggammathistbar->SetBinContent(i,backgroundtw_2);
    h5ftsignthistbar->SetBinContent(i, AntiParticle[k]);
    k+=2;
  }
  asimovtbar->Write("Asimov");
  h1sigthistbar->Write("tbar_channel");
  h3bgttbarthistbar->Write("ttbar");
  h4bgMultithistbar->Write("Multijet");
  h2bggammathistbar->Write("Electroweak");
  h5ftsignthistbar->Write("Modulation_tbar");
  tbarchannel->Close();

  return 0;

}
