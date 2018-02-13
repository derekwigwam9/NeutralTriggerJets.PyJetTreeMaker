// 'MakeRatiosOfFourThings.C'
// Derek Anderson
// 06.05.2017
//
// Does what's on the tin.


#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

using namespace std;


static const Int_t  nThings = 4;
static const Int_t  nRatios = 6;
static const Int_t  nCanvas = 3;
static const Int_t  nRebin  = 10;
static const Bool_t doRebin = true;
static const Bool_t doScale = true; 
// i/o parameters
static const TString sOutput("pp200r12.ratios.smallerBins.d6m6y2017.root");
static const TString sInput[nThings]  = {"input/pp200r12.tupleHistograms.d2m6y2017.root", "input/pp200r12.tupleHistograms.d2m6y2017.root", "input/cut.gam.firstEvtsExcluded.d30m5y2017.root", "input/cut.gam.allEvts.d5m6y2017.root"};
static const TString sHists[nThings]  = {"hPtPy", "hPtGe", "hPtTrk_NoCuts", "FitCut2/hPtTrkG_n2d0"};
// histogram parameters
static const TString  sNamesH[nThings] = {"hTrkPtP", "hTrkPtG", "hTrkPtU1", "hTrkPtU2"};
static const TString  sLegH[nThings]   = {"Pythia", "Geant", "MuDst (before cuts)", "MuDst (after cuts)"};
static const Double_t scale[nThings]   = {1., 1., (4424872. / 4531611.), (4424872. / 4531611.)};
// ratio parameters
static const TString sNamesR[nRatios] = {"hRatioGP", "hRatioU1P", "hRatioU2P", "hRatioU1G", "hRatioU2G", "hRatioU2U1"};
static const TString sLegR[nRatios]   = {"Geant / Pythia", "MuDst (before) / Pythia", "MuDst(after) / Pythia", "MuDst (before) / Geant", "MuDst (after) / Geant", "MuDst (after) / MuDst (before)"};
static const TString sCanR[nCanvas]   = {"cRatiosPy", "cRatiosGe", "cRatiosU1"};



void MakeRatiosOfFourThings() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning ratio-making script..." << endl;

  // global constants
  const TString sTitleXH("p_{T}^{trk}");
  const TString sTitleXR("p_{T}^{trk}");
  const TString sTitleYH("dN_{trk}/dp_{T}^{trk}");
  const TString sTitleYR("ratio");


  // open files
  TFile *fInput[nThings];
  for (Int_t iFiles = 0; iFiles < nThings; iFiles++) {
    fInput[iFiles] = new TFile(sInput[iFiles].Data(), "read");
    if (!fInput[iFiles]) {
      cerr << "PANIC: couldn't open file no. " << iFiles << endl;
      assert(fInput[iFiles]);
    }
  }
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  cout << "    Files opened." << endl;

  // grab histograms
  TH1D *hInput[nThings];
  for (Int_t iHist = 0; iHist < nThings; iHist++) {
    hInput[iHist] = (TH1D*) fInput[iHist] -> Get(sHists[iHist].Data());
    if (!hInput[iHist]) {
      cerr << "PANIC: couldn't grab histogram no. " << iHist << endl;
      assert(hInput[iHist]);
    }
  }
  cout << "    Histograms grabbed." << endl;


  // rebin histograms
  if (doRebin) {
    for (Int_t iHist = 0; iHist < nThings; iHist++) {
      hInput[iHist] -> Rebin(nRebin);
    }
    cout << "    Histograms rebinned." << endl;
  }


  // scale histograms
  if (doScale) {
    for (Int_t iHist = 0; iHist < nThings; iHist++) {
      hInput[iHist] -> Scale(1. / scale[iHist]);
    }
    cout << "    Histograms scaled." << endl;
  }


  // all histograms should have same dimensions
  Bool_t nBinsIsSame       = true;
  Bool_t xBin1IsSame       = true;
  Bool_t xBin2IsSame       = true;
  Bool_t dimensionsAreSame = true;
  for (Int_t iHist = 0; iHist < nThings; iHist++) {
    const Int_t    nBinsI = hInput[iHist] -> GetNbinsX();
    const Double_t xBin1I = hInput[iHist] -> GetBinLowEdge(1);
    const Double_t xBin2I = hInput[iHist] -> GetBinLowEdge(nBinsI + 1);
    for(Int_t jHist = 0; jHist < nThings; jHist++) {
      const Int_t    nBinsJ = hInput[jHist] -> GetNbinsX();
      const Double_t xBin1J = hInput[jHist] -> GetBinLowEdge(1);
      const Double_t xBin2J = hInput[jHist] -> GetBinLowEdge(nBinsJ + 1);
      if (nBinsI != nBinsJ) nBinsIsSame = false;
      if (xBin1I != xBin1J) xBin1IsSame = false;
      if (xBin2I != xBin2J) xBin2IsSame = false;
      if (!nBinsIsSame || !xBin1IsSame || !xBin2IsSame) {
        dimensionsAreSame = false;
        break;
      }
    }
  }
  if (!dimensionsAreSame) {
    cerr << "PANIC: histogram dimensions not the same!" << endl;
    assert(dimensionsSame);
  }
  const Int_t    nBins = hInput[0] -> GetNbinsX();
  const Double_t xBin1 = hInput[0] -> GetBinLowEdge(1);
  const Double_t xBin2 = hInput[0] -> GetBinLowEdge(nBins + 1);

  // calculate ratios
  TH1D *hRatios[nRatios];
  for (Int_t iRatio = 0; iRatio < nRatios; iRatio++) {
    hRatios[iRatio] = new TH1D(sNamesR[iRatio].Data(), "", nBins, xBin1, xBin2);
    switch (iRatio) {
      case 0:
        hRatios[iRatio] -> Divide(hInput[1], hInput[0], 1., 1.);
        break;
      case 1:
        hRatios[iRatio] -> Divide(hInput[2], hInput[0], 1., 1.);
        break;
      case 2:
        hRatios[iRatio] -> Divide(hInput[3], hInput[0], 1., 1.);
        break;
      case 3:
        hRatios[iRatio] -> Divide(hInput[2], hInput[1], 1., 1.);
        break;
      case 4:
        hRatios[iRatio] -> Divide(hInput[3], hInput[1], 1., 1.);
        break;
      case 5:
        hRatios[iRatio] -> Divide(hInput[3], hInput[2], 1., 1.);
        break;
    }
  }
  cout << "    Ratios calculated." << endl;


  // set styles
  const Int_t    cHist[nThings] = {810, 850, 890, 910};
  const Int_t    mHist[nThings] = {7, 4, 7, 4};
  const Int_t    fHist[nThings] = {3017, 3018, 3017, 3018};
  const Int_t    tHist[nThings] = {42, 42, 42, 42};
  const Int_t    eHist[nThings] = {1, 1, 1, 1};
  const Double_t aHist[nThings] = {0.02, 0.02, 0.02, 0.02};
  for (Int_t iHist = 0; iHist < nThings; iHist++) {
    hInput[iHist] -> SetName(sNamesH[iHist].Data());
    hInput[iHist] -> SetTitle("");
    hInput[iHist] -> SetTitleFont(tHist[iHist]);
    hInput[iHist] -> SetLineColor(cHist[iHist]);
    hInput[iHist] -> SetMarkerColor(cHist[iHist]);
    hInput[iHist] -> SetMarkerStyle(mHist[iHist]);
    hInput[iHist] -> GetXaxis() -> SetLabelSize(aHist[iHist]);
    hInput[iHist] -> GetXaxis() -> SetTitleFont(tHist[iHist]);
    hInput[iHist] -> GetXaxis() -> SetTitle(sTitleXH.Data());
    hInput[iHist] -> GetXaxis() -> CenterTitle(eHist[iHist]);
    hInput[iHist] -> GetYaxis() -> SetLabelSize(aHist[iHist]);
    hInput[iHist] -> GetYaxis() -> SetTitleFont(tHist[iHist]);
    hInput[iHist] -> GetYaxis() -> SetTitle(sTitleYH.Data());
    hInput[iHist] -> GetYaxis() -> CenterTitle(eHist[iHist]);
  }

  const Int_t    cRatios[nRatios] = {810, 850, 890, 810, 890, 910};
  const Int_t    mRatios[nRatios] = {7, 4, 7, 4, 7, 4};
  const Int_t    fRatios[nRatios] = {3017, 3018, 3017, 3018, 3017, 3018};
  const Int_t    tRatios[nRatios] = {42, 42, 42, 42, 42, 42};
  const Int_t    eRatios[nRatios] = {1, 1, 1, 1, 1, 1};
  const Double_t aRatios[nRatios] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  for (Int_t iRatio = 0; iRatio < nRatios; iRatio++) {
    hRatios[iRatio] -> SetTitle("");
    hRatios[iRatio] -> SetTitleFont(tRatios[iRatio]);
    hRatios[iRatio] -> SetLineColor(cRatios[iRatio]);
    hRatios[iRatio] -> SetMarkerColor(cRatios[iRatio]);
    hRatios[iRatio] -> SetMarkerStyle(mRatios[iRatio]);
    hRatios[iRatio] -> GetXaxis() -> SetLabelSize(aRatios[iRatio]);
    hRatios[iRatio] -> GetXaxis() -> SetTitleFont(tRatios[iRatio]);
    hRatios[iRatio] -> GetXaxis() -> SetTitle(sTitleXR.Data());
    hRatios[iRatio] -> GetXaxis() -> CenterTitle(eRatios[iRatio]);
    hRatios[iRatio] -> GetYaxis() -> SetLabelSize(aRatios[iRatio]);
    hRatios[iRatio] -> GetYaxis() -> SetTitleFont(tRatios[iRatio]);
    hRatios[iRatio] -> GetYaxis() -> SetTitle(sTitleYR.Data());
    hRatios[iRatio] -> GetYaxis() -> CenterTitle(eRatios[iRatio]);
  }
  cout << "    Styles set." << endl;


  // make legends
  const Int_t    cLeg = 0;
  const Int_t    fLeg = 0;
  const Int_t    lLeg = 0;
  const Int_t    tLeg = 42;
  const Double_t xL1  = 0.1;
  const Double_t yL1  = 0.1;
  const Double_t xL2  = 0.3;
  const Double_t yL2  = 0.3;
  TLegend *lInput  = new TLegend(xL1, yL1, xL2, yL2);
  lInput -> SetFillColor(cLeg);
  lInput -> SetFillStyle(fLeg);
  lInput -> SetLineColor(cLeg);
  lInput -> SetLineStyle(lLeg);
  lInput -> SetTextFont(tLeg);
  for (iHist = 0; iHist < nThings; iHist++) {
    lInput -> AddEntry(hInput[iHist], sLegH[iHist].Data());
  }

  TLegend *lRatio[nCanvas];
  for (Int_t iCan = 0; iCan < nCanvas; iCan++) {
    lRatio[iCan] = new TLegend(xL1, yL1, xL2, yL2);
    lRatio[iCan] -> SetFillColor(cLeg);
    lRatio[iCan] -> SetFillStyle(fLeg);
    lRatio[iCan] -> SetLineColor(cLeg);
    lRatio[iCan] -> SetLineStyle(lLeg);
    lRatio[iCan] -> SetTextFont(tLeg);
    switch (iCan) {
      case 0:
        lRatio[iCan] -> AddEntry(hRatios[0], sLegR[0].Data());
        lRatio[iCan] -> AddEntry(hRatios[1], sLegR[1].Data());
        lRatio[iCan] -> AddEntry(hRatios[2], sLegR[2].Data());
        break;
      case 1:
        lRatio[iCan] -> AddEntry(hRatios[3], sLegR[3].Data());
        lRatio[iCan] -> AddEntry(hRatios[4], sLegR[4].Data());
        break;
      case 2:
        lRatio[iCan] -> AddEntry(hRatios[5], sLegR[5].Data());
        break;
    }
  }
  cout << "    Legends created." << endl;


  // make plots
  const Int_t wCan = 800;
  const Int_t hCan = 800;
  const Int_t grid = 0;
  const Int_t log  = 1;
  TCanvas *cPlotI = new TCanvas("cDistributions", "", wCan, hCan);
  cPlotI    -> SetGrid(grid, grid);
  cPlotI    -> SetLogy(log);
  hInput[0] -> Draw();
  for (Int_t iHist = 1; iHist < nThings; iHist++) {
    hInput[iHist] -> Draw("same");
  }
  lInput -> Draw();
  cPlotI -> Write();
  cPlotI -> Close();

  TCanvas *cPlotR[nCanvas];
  for (Int_t iCan = 0; iCan < nCanvas; iCan++) {
    cPlotR[iCan] = new TCanvas(sCanR[iCan].Data(), "", wCan, hCan);
    cPlotR[iCan] -> SetGrid(grid, grid);
    switch (iCan) {
      case 0:
        hRatios[0]   -> Draw();
        hRatios[1]   -> Draw("same");
        hRatios[2]   -> Draw("same");
        lRatio[iCan] -> Draw();
        break;
      case 1:
        hRatios[3]   -> Draw();
        hRatios[4]   -> Draw("same");
        lRatio[iCan] -> Draw();
        break;
      case 2:
        hRatios[5]   -> Draw();
        lRatio[iCan] -> Draw();
        break;
    }
    cPlotR[iCan] -> Write();
    cPlotR[iCan] -> Close();
  }
  cout << "    Plots made." << endl;


  // close files
  fOutput -> cd();
  for (Int_t iHist = 0; iHist < nThings; iHist++) {
    hInput[iHist] -> Write();
  }
  for (Int_t iRatio = 0; iRatio < nRatios; iRatio++) {
    hRatios[iRatio] -> Write();
  }
  fOutput -> Close();
  for (Int_t iFiles = 0; iFiles < nThings; iFiles++) {
    fInput[iFiles] -> cd();
    fInput[iFiles] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
