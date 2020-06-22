#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>
using namespace std;

int main (){
  double asym;
  string b_mu;
  cout<<"Give the value of b_mu you used for the data you want to study"<<endl;
  cin>>b_mu;
  TFile* File = TFile::Open("results/modulation/Modulation_Temporelle_13TeV_b1.root");
  TFile* Save = new TFile ("Asymmetry.root","RECREATE");
  TCanvas * modul = new TCanvas("Asymmetry", "Simple graph", 450 , 450);
  TH2F* axe = new TH2F (" ", " ", 49 ,0,24,100,0.265,0.285 );
  auto legend = new TLegend(0.7,0.75,0.95,0.95);
  modul->SetLeftMargin(0.14);
  modul->SetRightMargin(0.03);


  axe->SetStats(kFALSE);
  axe->Draw();
  axe->GetXaxis()->SetTitle("t (h)");
  axe->GetYaxis()->SetTitle("Charge Asymmetry");
  axe->GetYaxis()->CenterTitle(kTRUE);
  axe->GetXaxis()->CenterTitle(kTRUE);
  axe->GetYaxis()->SetTitleOffset(2.0);
  Save->cd();

  TGraph* Parti = (TGraph*)File->Get("tX");
  TGraph* AntiParti = (TGraph*)File->Get("tbarX");
  TGraph* Asym = new TGraph (48);
  double* Particle = Parti->GetY();
  double* AntiParticle = AntiParti->GetY();

  double Xsect = 12.26;
  double Xsectbar = 7.004;

  /*cout<<"Give the Xsect of the top production events"<<endl;
  cin>>Xsect;
  cout<<"Give the Xsect of the anti-top production events"<<endl;
  cin>>Xsectbar;*/
  double l = 0;
for (int i = 0; i<=48; i++ )
  {
    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    cout<<"Asymx="<<asym<<endl;
    l+=0.5;
  }
  string name = "b_{x} = "+b_mu+" GeV";
  Asym->SetLineColor(kOrange);
  Asym->SetLineWidth(3);
  Asym->Write("tX");
  legend->AddEntry(Asym,name.c_str(), "l");
  Asym->Draw("Same");

  Parti = (TGraph*)File->Get("tY");
  AntiParti = (TGraph*)File->Get("tbarY");
  Particle = Parti->GetY();
  AntiParticle = AntiParti->GetY();
  Asym = new TGraph (48);

  l=0;
  for (int i = 0; i <=48; i++)
  {
    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(kRed);
  Asym->SetLineWidth(3);
  Asym->Write("tY");
  name = "b_{y} = "+b_mu+" GeV";
  legend->AddEntry(Asym,name.c_str(), "l");
  Asym->Draw("Same");


  Parti = (TGraph*)File->Get("tZ");
  AntiParti = (TGraph*)File->Get("tbarZ");
  Particle = Parti->GetY();
  AntiParticle = AntiParti->GetY();
  Asym = new TGraph (48);

l=0;
  for (int i = 0; i <=48; i++ )
  {

    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(kBlue);
  Asym->SetLineWidth(3);
  Asym->Write("tZ");
  name = "b_{z} = "+b_mu+" GeV";
  legend->AddEntry(Asym,name.c_str(), "l");
  Asym->Draw("Same");


  Parti = (TGraph*)File->Get("tT");
  AntiParti = (TGraph*)File->Get("tbarT");
  Particle = Parti->GetY();
  AntiParticle = AntiParti->GetY();
  Asym = new TGraph (48);

l=0;
  for (int i = 0; i <= 48; i++)
  {
    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(8);
  Asym->SetLineWidth(3);
  Asym->Write("tT");
  name = "b_{t} = "+b_mu+" GeV";
  legend->AddEntry(Asym,name.c_str(), "l");
  legend->Draw("Same");
  Asym->Draw("Same");

  Parti = (TGraph*)File->Get("tT");
  Asym = new TGraph (48);
  l=0;
  for (int i = 0; i <= 48; i++)
  {
    asym = (Xsect - Xsectbar)/(Xsect + Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(kBlack);
  Asym->SetLineWidth(3);
  Asym->Write("SM");
  legend->AddEntry(Asym,"SM", "l");
  legend->Draw("Same");
  Asym->Draw("Same");

modul->SaveAs("Asymmetry_t.pdf");

  return 0;
}
