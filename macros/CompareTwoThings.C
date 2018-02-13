// 'CompareTwoThings.C'
// Derek Anderson
// 06.28.2017
//
// Use this to compare a couple of histograms.


#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// i/o parameters
static const TString sOutput("dataXembedding.eTwrForEp.evtScaled.d17m8y2017.root");
static const TString sInputA("../JetData/pp200r9.fullSample.EtwrCheck.hotTowersRemoved.d16m8y2017.root");
static const TString sInputB("../JetData/pp200r12.pt5_-1.EtwrCheck.hotTowersRemoved.d16m8y2017.root");
static const TString sHistA("EoverP/hEpTwrEne");
static const TString sHistB("EoverP/hEpTwrEne");
// histogram parameters
static const TString sNameA("hTwrEne_data");
static const TString sNameB("hTwrEne_embed");
static const TString sTitle("E_{twr} (for E/p calculation)");
static const TString sTitleX("E_{twr}");
static const TString sTitleY("(1/N_{evt}) dN_{twr}/dE_{twr}");
// legend parameters
static const TString sLegendA("data");
static const TString sLegendB("embedding");
// label parameters
static const TString sLabel1("Requiring one match");
static const TString sLabel2("Scaled by N_{evt}");
static const TString sLabel3("E_{twr} > 0.2 GeV/c, E_{twr}^{corr} > 0.2 GeV/c, |#eta^{twr}| < 1.0");
static const TString sLabel4("Anti-k_{T}, R = 0.3, full jets");
// canvas parameters
static const TString sCanvas("cTwrEne");



void CompareTwoThings() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  if (!fOutput || !fInputA || !fInputB) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  if (!hInputA || !hInputB) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
  }
  cout << "    Histograms grabbed." << endl;


  const Int_t    cA  = 890;
  const Int_t    cB  = 810;
  const Int_t    mA  = 4;
  const Int_t    mB  = 1;
  const Int_t    txt = 42;
  const Int_t    cnt = 1;
  const Double_t lab = 0.02;
  hInputA -> SetLineColor(cA);
  hInputA -> SetMarkerColor(cA);
  hInputA -> SetMarkerStyle(mA);
  hInputA -> SetTitleFont(txt);
  hInputA -> SetTitle(sTitle.Data());
  hInputA -> SetName(sNameA.Data());
  hInputA -> GetXaxis() -> SetLabelSize(lab);
  hInputA -> GetXaxis() -> CenterTitle(cnt);
  hInputA -> GetXaxis() -> SetTitleFont(txt);
  hInputA -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputA -> GetYaxis() -> SetLabelSize(lab);
  hInputA -> GetYaxis() -> CenterTitle(cnt);
  hInputA -> GetYaxis() -> SetTitleFont(txt);
  hInputA -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputB -> SetLineColor(cB);
  hInputB -> SetMarkerColor(cB);
  hInputB -> SetMarkerStyle(mB);
  hInputB -> SetTitleFont(txt);
  hInputB -> SetTitle(sTitle.Data());
  hInputB -> SetName(sNameB.Data());
  hInputB -> GetXaxis() -> SetLabelSize(lab);
  hInputB -> GetXaxis() -> CenterTitle(cnt);
  hInputB -> GetXaxis() -> SetTitleFont(txt);
  hInputB -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputB -> GetYaxis() -> SetLabelSize(lab);
  hInputB -> GetYaxis() -> CenterTitle(cnt);
  hInputB -> GetYaxis() -> SetTitleFont(txt);
  hInputB -> GetYaxis() -> SetTitle(sTitleY.Data());
  cout << "    Styles set." << endl;


  // scale histograms
  const Bool_t   scaleA = true;
  const Bool_t   scaleB = true;
  const Double_t aScale = 1429435.;
  const Double_t bScale = 1336827.;
  if (scaleA) {
    hInputA -> Scale(1. / aScale);
    cout << "    Input A scaled." << endl;
  }
  if (scaleB) {
    hInputB -> Scale(1. / bScale);
    cout << "    Input B scaled." << endl;
  }


  const Int_t    cL  = 0;
  const Int_t    fL  = 0;
  const Int_t    sL  = 0;
  const Double_t x1L = 0.1;
  const Double_t x2L = 0.3;
  const Double_t y1L = 0.1;
  const Double_t y2L = 0.3;
  TLegend *lLegend = new TLegend(x1L, y1L, x2L, y2L);
  lLegend -> SetFillColor(cL);
  lLegend -> SetFillStyle(sL);
  lLegend -> SetLineColor(cL);
  lLegend -> SetLineStyle(sL);
  lLegend -> SetTextFont(txt);
  lLegend -> AddEntry(hInputA, sLegendA.Data());
  lLegend -> AddEntry(hInputB, sLegendB.Data());
  cout << "    Legend created." << endl;

  const Double_t x1P = 0.3;
  const Double_t x2P = 0.5;
  const Double_t y1P = 0.1;
  const Double_t y2P = 0.3;
  TPaveText *pLabel = new TPaveText(x1P, y1P, x2P, y2P, "NDC NB");
  pLabel -> SetFillColor(cL);
  pLabel -> SetFillStyle(sL);
  pLabel -> SetLineColor(cL);
  pLabel -> SetLineStyle(sL);
  pLabel -> SetTextFont(txt);
  pLabel -> AddText(sLabel1.Data());
  pLabel -> AddText(sLabel2.Data());
  //pLabel -> AddText(sLabel3.Data());
  //pLabel -> AddText(sLabel4.Data());
  cout << "    Label created." << endl;


  // make plots
  fOutput -> cd();

  const Int_t wC  = 800;
  const Int_t hC  = 800;
  const Int_t grd = 0;
  const Int_t log = 0;
  TCanvas *cPlot = new TCanvas(sCanvas.Data(), "", wC, hC);
  cPlot   -> SetGrid(grd, grd);
  cPlot   -> SetLogy(log);
  hInputA -> Draw();
  hInputB -> Draw("same");
  lLegend -> Draw();
  //pLabel  -> Draw();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Plot drawn." << endl;


  fOutput -> cd();
  hInputA -> Write();
  hInputB -> Write();
  fOutput -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputA -> cd();
  fInputA -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
