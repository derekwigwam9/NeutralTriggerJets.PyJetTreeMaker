// 'CompareBeforeAndAfter.C'
// Derek Anderson
// 05.16.2017
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
static const TString sOutput("pp200r9.beforeXafterPi0.dFjet.r03a02rm1chrg.d16m5y2017.root");
static const TString sBefore("raw/pp200r9.before.r03a02rm1chrg.plots.d15m5y2017.root");
static const TString sAfter("raw/pp200r9.after.r03a02rm1chrg.plots.d15m5y2017.root");
static const TString sHistB("Pi0/hAllDeltaPhiP");
static const TString sHistA("Pi0/hAllDeltaPhiP");
// histogram parameters
static const TString sNameB("hJetDfP_Before");
static const TString sNameA("hJetDfP_After");
static const TString sTitle("Jet #Delta#varphi, #pi^{0} trigger");
static const TString sTitleX("#Delta#varphi = #varphi_{jet} - #varphi_{trg}");
static const TString sTitleY("(1/N_{trg}) dN_{jet}/(d#Delta#varphi d#eta^{jet})");
// legend parameters
static const TString sLegendB("Before fix");
static const TString sLegendA("After fix");
// label parameters
static const TString sLabel1("pp-collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9, 30) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} #in (0.2, 30) GeV/c, |#eta^{trk}| < 1.0");
static const TString sLabel4("Anti-k_{T}, R = 0.3, charged jet");
// canvas parameters
static const TString sCanvas("hJetDfPi0");



void CompareBeforeAndAfter() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fBefore = new TFile(sBefore.Data(), "read");
  TFile *fAfter  = new TFile(sAfter.Data(), "read");
  if (!fOutput || !fBefore || !fAfter) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fBefore);
    assert(fAfter);
  }
  cout << "    Files opened." << endl;

  TH1D *hBefore = (TH1D*) fBefore -> Get(sHistB.Data());
  TH1D *hAfter  = (TH1D*) fAfter  -> Get(sHistA.Data());
  if (!hBefore || !hAfter) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hBefore);
    assert(hAfter);
  }
  cout << "    Histograms grabbed." << endl;


  const Int_t    cB  = 890;
  const Int_t    cA  = 810;
  const Int_t    mB  = 20;
  const Int_t    mA  = 24;
  const Int_t    txt = 42;
  const Int_t    cnt = 1;
  const Double_t lab = 0.02;
  hBefore -> SetLineColor(cB);
  hBefore -> SetMarkerColor(cB);
  hBefore -> SetMarkerStyle(mB);
  hBefore -> SetTitleFont(txt);
  hBefore -> SetTitle(sTitle.Data());
  hBefore -> SetName(sNameB.Data());
  hBefore -> GetXaxis() -> SetLabelSize(lab);
  hBefore -> GetXaxis() -> CenterTitle(cnt);
  hBefore -> GetXaxis() -> SetTitleFont(txt);
  hBefore -> GetXaxis() -> SetTitle(sTitleX.Data());
  hBefore -> GetYaxis() -> SetLabelSize(lab);
  hBefore -> GetYaxis() -> CenterTitle(cnt);
  hBefore -> GetYaxis() -> SetTitleFont(txt);
  hBefore -> GetYaxis() -> SetTitle(sTitleY.Data());
  hAfter  -> SetLineColor(cA);
  hAfter  -> SetMarkerColor(cA);
  hAfter  -> SetMarkerStyle(mA);
  hAfter  -> SetTitleFont(txt);
  hAfter  -> SetTitle(sTitle.Data());
  hAfter  -> SetName(sNameA.Data());
  hAfter  -> GetXaxis() -> SetLabelSize(lab);
  hAfter  -> GetXaxis() -> CenterTitle(cnt);
  hAfter  -> GetXaxis() -> SetTitleFont(txt);
  hAfter  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hAfter  -> GetYaxis() -> SetLabelSize(lab);
  hAfter  -> GetYaxis() -> CenterTitle(cnt);
  hAfter  -> GetYaxis() -> SetTitleFont(txt);
  hAfter  -> GetYaxis() -> SetTitle(sTitleY.Data());
  cout << "    Styles set." << endl;

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
  lLegend -> AddEntry(hBefore, sLegendB.Data());
  lLegend -> AddEntry(hAfter, sLegendA.Data());
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
  pLabel -> AddText(sLabel3.Data());
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
  hBefore -> Draw();
  hAfter  -> Draw("same");
  lLegend -> Draw();
  pLabel  -> Draw();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Plot drawn." << endl;


  fOutput -> cd();
  hBefore -> Write();
  hAfter  -> Write();
  fOutput -> Close();
  fBefore -> cd();
  fBefore -> Close();
  fAfter  -> cd();
  fAfter -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
