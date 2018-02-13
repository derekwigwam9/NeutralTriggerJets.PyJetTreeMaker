// 'CompareThreeThrings.C'
// Derek Anderson
// 05.16.2017
//
// Use this to compare three histograms.


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
static const TString sOutput("pythiaXdata.pTtrkPi0.eff88.r03a02rm1chrg.d5m7y2017.root");
static const TString sInputA("output/EfficiencyStudy/pythiaPP.r03a02rm1chrg.d19m5y2017.root");
static const TString sInputB("output/EfficiencyStudy/pythiaPD.eff88.r03rm1chrg.d2m7y2017.root");
static const TString sInputC("../JetData/pp200r9.r03rm1chrg.d28m6y2017.root");
static const TString sHistA("QA/hPtTrk");
static const TString sHistB("QA/hPtTrk");
static const TString sHistC("QA/Pi0/hTrkPtP");
// histogram parameters
static const TString sNameA("hTrkPtP_par");
static const TString sNameB("hTrkPtP_det");
static const TString sNameC("hTrkPtP_data");
static const TString sTitle("Track p_{T}, #pi^{0} trigger");
static const TString sTitleX("p_{T}^{trk}");
static const TString sTitleY("(1/N_{trg}) dN_{trk}/dp_{T}^{trk}");
// legend parameters
static const TString sLegendA("Pythia, particle");
static const TString sLegendB("Pythia, w/ #epsilon_{0} = 88%");
static const TString sLegendC("Run 9");
// label parameters
static const TString sLabel1("pp-collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9, 30) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} #in (0.2, 30) GeV/c, |#eta^{trk}| < 1.0");
static const TString sLabel4("Anti-k_{T}, R = 0.3, charged jet");
// canvas parameters
static const TString sCanvas("cTrkPtPi0");



void CompareThreeThings() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  TFile *fInputC = new TFile(sInputC.Data(), "read");
  if (!fOutput || !fInputA || !fInputB || !fInputC) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
    assert(fInputC);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  TH1D *hInputC = (TH1D*) fInputC -> Get(sHistC.Data());
  if (!hInputA || !hInputB || !hInputC) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
    assert(hInputC);
  }
  cout << "    Histograms grabbed." << endl;


  const Int_t    cA  = 890;
  const Int_t    cB  = 850;
  const Int_t    cC  = 810;
  const Int_t    mA  = 20;
  const Int_t    mB  = 24;
  const Int_t    mC  = 7;
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
  hInputC -> SetLineColor(cC);
  hInputC -> SetMarkerColor(cC);
  hInputC -> SetMarkerStyle(mC);
  hInputC -> SetTitleFont(txt);
  hInputC -> SetTitle(sTitle.Data());
  hInputC -> SetName(sNameC.Data());
  hInputC -> GetXaxis() -> SetLabelSize(lab);
  hInputC -> GetXaxis() -> CenterTitle(cnt);
  hInputC -> GetXaxis() -> SetTitleFont(txt);
  hInputC -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputC -> GetYaxis() -> SetLabelSize(lab);
  hInputC -> GetYaxis() -> CenterTitle(cnt);
  hInputC -> GetYaxis() -> SetTitleFont(txt);
  hInputC -> GetYaxis() -> SetTitle(sTitleY.Data());
  cout << "    Styles set." << endl;

  // rebin input
  const UInt_t nRebinA = 5;
  const UInt_t nRebinB = 5;
  const UInt_t nRebinC = 5;
  const UInt_t rebinA  = true;
  const UInt_t rebinB  = true;
  const Bool_t rebinC  = false;
  if (rebinA) {
    hInputA -> Rebin(nRebinA);
    hInputA -> Scale(1. / (Double_t) nRebinA);
    cout << "    Rebinned input A." << endl;
  }
  if (rebinB) {
    hInputB -> Rebin(nRebinB);
    hInputB -> Scale(1. / (Double_t) nRebinB);
    cout << "    Rebinned input B." << endl;
  }
  if (rebinC) {
    hInputC -> Rebin(nRebinC);
    hInputC -> Scale(1. / (Double_t) nRebinC);
    cout << "    Rebinned input C." << endl;
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
  lLegend -> AddEntry(hInputC, sLegendC.Data());
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
  pLabel -> AddText(sLabel4.Data());
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
  hInputC -> Draw("same");
  lLegend -> Draw();
  pLabel  -> Draw();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Plot drawn." << endl;


  fOutput -> cd();
  hInputA -> Write();
  hInputB -> Write();
  hInputC -> Write();
  fOutput -> Close();
  fInputA -> cd();
  fInputA -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputC -> cd();
  fInputC -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
