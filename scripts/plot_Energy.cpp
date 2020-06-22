#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TFile.h>
#include <TAttMarker.h>
#include <TAttPad.h>
using namespace std;

int main ()
{
  TCanvas * energy = new TCanvas("", "Simple graph", 450 , 450);
  energy->SetLogx();

  TH2F* axe = new TH2F (" ", "", 7 ,2,150,7, 0,0.1);
  auto legend = new TLegend(0.6,0.15,0.85,0.25);
  axe->SetMarkerStyle(22);
  axe->SetStats(kFALSE);

  axe->GetXaxis()->SetTitle("Center-of-Mass energy (TeV)");
  axe->GetXaxis()->CenterTitle(kTRUE);
  axe->GetYaxis()->SetTitle("Asymmetry (E)");
  axe->GetYaxis()->CenterTitle(kTRUE);
  axe->GetYaxis()->SetTitleOffset(1.35);
  axe->Draw();

  TGraph* Depend = new TGraph (7);
  //tw dependency
  Depend->SetMarkerSize(1);
  Depend->SetPoint(0,2,0.00892);
  Depend->SetPoint(1,7,0.0211);
  Depend->SetPoint(2,13,0.0284);
  Depend->SetPoint(3,14,0.0296);
  Depend->SetPoint(4,27,0.0371);
  Depend->SetPoint(5,100,0.0496);
  Depend->SetPoint(6,150,0.0539942);
  //TGraph* Depend = new TGraph (7);
  //t Energy_Dependency
  /*Depend->SetMarkerSize(1);
  Depend->SetPoint(0,2,0.276621);
  Depend->SetPoint(1,7,0.278291);
  Depend->SetPoint(2,13,0.278911);
  Depend->SetPoint(3,14,0.279102);
  Depend->SetPoint(4,27,0.279197);
  Depend->SetPoint(5,100,0.278734);
  Depend->SetPoint(6,150,0.27855);*/



  Depend->SetLineColor(kRed);
  legend->AddEntry(Depend, "b_{x} = 1 GeV", "l");
  Depend->Draw("Same");
  Depend->Draw("*");
  legend->Draw("Same");
  energy->SaveAs("Energy_Dependency.pdf");
  /*double x[25];
  double y[25];
  double h = 0;
  TH2F* axe = new TH2F (" ", "(5 + 20*f(t)/(20+5*f(t))", 24 ,0,24,24, -2,2);
  axe->GetXaxis()->SetTitle("t");
  axe->GetXaxis()->CenterTitle(kTRUE);
  axe->GetYaxis()->SetTitle("f(t)");
  axe->GetYaxis()->CenterTitle(kTRUE);
  axe->GetYaxis()->SetTitleOffset(1.35);
  axe->Draw();
  axe->SetStats(kFALSE);
  for (int i = 0 ; i<=24 ; i++)
  {
    x[i]=h;
    y[i]= (5 + 20*(sin(2*M_PI/24*i) + cos(2*M_PI/24*i) ))/(20 + 5*(cos(2*M_PI/24*i)+sin(2*M_PI/24*i)));
    h++;
  }
    TGraph* Sinus = new TGraph (24, x , y);
    Sinus->Draw("Same");
    energy->SaveAs("Cosinus.pdf");*/
  return 0;
}
