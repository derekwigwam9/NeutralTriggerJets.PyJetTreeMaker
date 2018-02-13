// 'CompareFourThrings.C'
// Derek Anderson
// 07.25.2017
//
// Use this to compare four histograms.


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
static const TString sOutput("dataXembedding.pTtrk.lessThan5.d16m8y2017.root");
static const TString sInputA("../JetData/TowerResponseStudy/Run9/pp200r9.fullSample.hotTowersRemoved.d14m8y2017.root");
static const TString sInputB("../JetData/pp200r12.pt5_-1.hotTowersRemoved.d15m8y2017.root");
static const TString sInputC("../JetData/TowerResponseStudy/Run9/pp200r9.fullSample.hotTowersRemoved.d14m8y2017.root");
static const TString sInputD("../JetData/pp200r12.pt5_-1.hotTowersRemoved.d15m8y2017.root");
static const TString sHistA("EoverP/hEpTrkPtH0");
static const TString sHistB("EoverP/hEpTrkPtH0");
static const TString sHistC("EoverP/hEpTrkPtE0");
static const TString sHistD("EoverP/hEpTrkPtE0");
// histogram parameters
static const TString sNameA("hPtHadr_data");
static const TString sNameB("hPtHadr_embed");
static const TString sNameC("hPtElec_data");
static const TString sNameD("hPtElec_embed");
static const TString sTitle("Track p_{T}, e^{#pm} vs. h^{#pm}");
static const TString sTitleX("p_{T}");
static const TString sTitleY("arb.");
// legend parameters
static const TString sLegendA("h^{#pm}, data");
static const TString sLegendB("h^{#pm}, embedding");
static const TString sLegendC("e^{#pm}, data");
static const TString sLegendD("e^{#pm}, embedding");
// label parameters
static const TString sLabel1("Requiring one match");
static const TString sLabel2("E_{twr} < 5 GeV");
static const TString sLabel3("p_{T}^{e,h} #in (0,100) GeV/c");
static const TString sLabel4("Anti-k_{T}, R = 0.3, charged jet");
// canvas parameters
static const TString sCanvas("cPtEleVsHad");



void CompareFourThings() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  TFile *fInputC = new TFile(sInputC.Data(), "read");
  TFile *fInputD = new TFile(sInputD.Data(), "read");
  if (!fOutput || !fInputA || !fInputB || !fInputC || !fInputD) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
    assert(fInputC);
    assert(fInputD);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  TH1D *hInputC = (TH1D*) fInputC -> Get(sHistC.Data());
  TH1D *hInputD = (TH1D*) fInputD -> Get(sHistD.Data());
  if (!hInputA || !hInputB || !hInputC || !hInputD) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
    assert(hInputC);
    assert(hInputD);
  }
  cout << "    Histograms grabbed." << endl;


  const Int_t    cA  = 810;
  const Int_t    cB  = 800;
  const Int_t    cC  = 890;
  const Int_t    cD  = 880;
  const Int_t    mA  = 1;
  const Int_t    mB  = 4;
  const Int_t    mC  = 1;
  const Int_t    mD  = 4;
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
  hInputD -> SetLineColor(cD);
  hInputD -> SetMarkerColor(cD);
  hInputD -> SetMarkerStyle(mD);
  hInputD -> SetTitleFont(txt);
  hInputD -> SetTitle(sTitle.Data());
  hInputD -> SetName(sNameD.Data());
  hInputD -> GetXaxis() -> SetLabelSize(lab);
  hInputD -> GetXaxis() -> CenterTitle(cnt);
  hInputD -> GetXaxis() -> SetTitleFont(txt);
  hInputD -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputD -> GetYaxis() -> SetLabelSize(lab);
  hInputD -> GetYaxis() -> CenterTitle(cnt);
  hInputD -> GetYaxis() -> SetTitleFont(txt);
  hInputD -> GetYaxis() -> SetTitle(sTitleY.Data());
  cout << "    Styles set." << endl;


  // rebin input
  const UInt_t nRebinA = 5;
  const UInt_t nRebinB = 5;
  const UInt_t nRebinC = 5;
  const UInt_t nRebinD = 5;
  const UInt_t rebinA  = false;
  const UInt_t rebinB  = false;
  const Bool_t rebinC  = false;
  const Bool_t rebinD  = false;
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
  if (rebinD) {
    hInputD -> Rebin(nRebinD);
    hInputD -> Scale(1. / (Double_t) nRebinD);
    cout << "    Rebinned input D." << endl;
  }

  // scale input
  const Bool_t scaleByEntries = true;
  if (scaleByEntries) {
    const UInt_t nA = hInputA -> GetEntries();
    const UInt_t nB = hInputB -> GetEntries();
    const UInt_t nC = hInputC -> GetEntries();
    const UInt_t nD = hInputD -> GetEntries();
    hInputA -> Scale(1. / (Double_t) nA);
    hInputB -> Scale(1. / (Double_t) nB);
    hInputC -> Scale(1. / (Double_t) nC);
    hInputD -> Scale(1. / (Double_t) nD);
    cout << "    Scaled input." << endl;
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
  lLegend -> AddEntry(hInputD, sLegendD.Data());
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
  hInputC -> Draw("same");
  hInputD -> Draw("same");
  lLegend -> Draw();
  pLabel  -> Draw();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Plot drawn." << endl;


  fOutput -> cd();
  hInputA -> Write();
  hInputB -> Write();
  hInputC -> Write();
  hInputD -> Write();
  fOutput -> Close();
  fInputA -> cd();
  fInputA -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputC -> cd();
  fInputC -> Close();
  fInputD -> cd();
  fInputD -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
