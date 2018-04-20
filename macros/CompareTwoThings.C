// 'CompareTwoThings.C'
// Derek Anderson
// 06.28.2017
//
// Use this to compare a couple of histograms.


#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// i/o parameters
static const TString sOutput("resolution.pTtrkCheck.c5p1.d20m4y2019.root");
static const TString sInputA("pp200r9.resCheck.r03rm1chrg.d20m4y2018.root");
static const TString sInputB("pp200py.resTestC5P1.eTtrg920pi0det.r03rm1chrg.d19m4y2018.root");
static const TString sHistA("QA/Pi0/hTrkPtP");
static const TString sHistB("QA/hPtTrk");
// histogram parameters
static const TString sNameA("hData");
static const TString sNameB("hPythia");
static const TString sNameR("hRatio");
static const TString sTitle("Track p_{T}");
static const TString sTitleX("p_{T}^{trk}");
static const TString sTitleY("(1/N^{trg}) dN^{trk}/dp_{T}^{trk}");
static const TString sTitleR("data / pythia");
static const Float_t xRange[2] = {-0.5, 20.5};
// legend parameters
static const TString sLegendA("data");
static const TString sLegendB("pythia8");
// label parameters
static const TString sLabel1("pp-collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("#pi^{0} trigger");
static const TString sLabel3("E^{trg}_{T} #in (9, 20) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel4("#varsigma = 0.05 #oplus 0.01");
// canvas parameters
static const TString sCanvas("cTrkPt");
static const TString sPad[2] = {"pRatio", "pTracks"};



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


  // scale histograms
  const Bool_t  scaleA = false;
  const Bool_t  scaleB = false;
  const Float_t aScale = 1429435.;
  const Float_t bScale = 1336827.;
  if (scaleA) {
    hInputA -> Scale(1. / aScale);
    cout << "    Input A scaled." << endl;
  }
  if (scaleB) {
    hInputB -> Scale(1. / bScale);
    cout << "    Input B scaled." << endl;
  }


  // rebin histograms
  const Bool_t  rebinA  = false;
  const Bool_t  rebinB  = true;
  const UInt_t  nRebinA = 2;
  const UInt_t  nRebinB = 5;

  Float_t oldBin(0.);
  Float_t newBin(0.);
  if (rebinA) {
    oldBin = hInputA -> GetBinWidth(17);
    hInputA -> Rebin(nRebinA);
    newBin = hInputA -> GetBinWidth(17);
    hInputA -> Scale(oldBin / newBin);
    cout << "    Input A rebinned." << endl;
  }
  if (rebinB) {
    oldBin = hInputB -> GetBinWidth(17);
    hInputB -> Rebin(nRebinB);
    newBin = hInputB -> GetBinWidth(17);
    hInputB -> Scale(oldBin / newBin);
    cout << "    InputB rebinned." << endl;
  }


  // calculate ratio
  const UInt_t  nBins   = hInputA -> GetNbinsX();
  const Float_t bin1    = hInputA -> GetBinLowEdge(1);
  const Float_t bin2    = hInputA -> GetBinLowEdge(nBins + 1);
  const Float_t weightA = 1.;
  const Float_t weightB = 1.;

  TH1D *hRatio = new TH1D("hRatio", "", nBins, bin1, bin2);
  hRatio -> Sumw2();
  hRatio -> Divide(hInputA, hInputB, weightA, weightB);
  cout << "    Calculated ratio." << endl;


  const Int_t   cA   = 896;
  const Int_t   cB   = 856;
  const Int_t   cR   = 876;
  const Int_t   mA   = 7;
  const Int_t   mB   = 4;
  const UInt_t  mR   = 8;
  const Int_t   txt  = 42;
  const Int_t   cnt  = 1;
  const Float_t off  = 1.;
  const Float_t labI = 0.02;
  const Float_t labR = 0.045;
  hInputA -> SetLineColor(cA);
  hInputA -> SetMarkerColor(cA);
  hInputA -> SetMarkerStyle(mA);
  hInputA -> SetTitleFont(txt);
  hInputA -> SetTitle(sTitle.Data());
  hInputA -> SetName(sNameA.Data());
  hInputA -> GetXaxis() -> SetLabelSize(labI);
  hInputA -> GetXaxis() -> CenterTitle(cnt);
  hInputA -> GetXaxis() -> SetTitleFont(txt);
  hInputA -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputA -> GetXaxis() -> SetTitleOffset(off);
  hInputA -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hInputA -> GetYaxis() -> SetLabelSize(labI);
  hInputA -> GetYaxis() -> CenterTitle(cnt);
  hInputA -> GetYaxis() -> SetTitleFont(txt);
  hInputA -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputB -> SetLineColor(cB);
  hInputB -> SetMarkerColor(cB);
  hInputB -> SetMarkerStyle(mB);
  hInputB -> SetTitleFont(txt);
  hInputB -> SetTitle(sTitle.Data());
  hInputB -> SetName(sNameB.Data());
  hInputB -> GetXaxis() -> SetLabelSize(labI);
  hInputB -> GetXaxis() -> CenterTitle(cnt);
  hInputB -> GetXaxis() -> SetTitleFont(txt);
  hInputB -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputB -> GetXaxis() -> SetTitleOffset(off);
  hInputB -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hInputB -> GetYaxis() -> SetLabelSize(labI);
  hInputB -> GetYaxis() -> CenterTitle(cnt);
  hInputB -> GetYaxis() -> SetTitleFont(txt);
  hInputB -> GetYaxis() -> SetTitle(sTitleY.Data());
  hRatio  -> SetLineColor(cR);
  hRatio  -> SetMarkerColor(mR);
  hRatio  -> SetName(sNameR.Data());
  hRatio  -> GetXaxis() -> SetLabelSize(labR);
  hRatio  -> GetXaxis() -> CenterTitle(cnt);
  hRatio  -> GetXaxis() -> SetTitleFont(txt);
  hRatio  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatio  -> GetXaxis() -> SetTitleOffset(off);
  hRatio  -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hRatio  -> GetYaxis() -> SetLabelSize(labR);
  hRatio  -> GetYaxis() -> CenterTitle(cnt);
  hRatio  -> GetYaxis() -> SetTitleFont(txt);
  hRatio  -> GetYaxis() -> SetTitle(sTitleR.Data());
  cout << "    Styles set." << endl;


  const Int_t   cL  = 0;
  const Int_t   fL  = 0;
  const Int_t   sL  = 0;
  const Float_t x1L = 0.1;
  const Float_t x2L = 0.3;
  const Float_t y1L = 0.1;
  const Float_t y2L = 0.3;
  TLegend *lLegend = new TLegend(x1L, y1L, x2L, y2L);
  lLegend -> SetFillColor(cL);
  lLegend -> SetFillStyle(sL);
  lLegend -> SetLineColor(cL);
  lLegend -> SetLineStyle(sL);
  lLegend -> SetTextFont(txt);
  lLegend -> AddEntry(hInputA, sLegendA.Data());
  lLegend -> AddEntry(hInputB, sLegendB.Data());
  cout << "    Legend created." << endl;

  const UInt_t  align = 12;
  const Float_t x1P   = 0.3;
  const Float_t x2P   = 0.5;
  const Float_t y1P   = 0.1;
  const Float_t y2P   = 0.3;
  TPaveText *pLabel = new TPaveText(x1P, y1P, x2P, y2P, "NDC NB");
  pLabel -> SetFillColor(cL);
  pLabel -> SetFillStyle(sL);
  pLabel -> SetLineColor(cL);
  pLabel -> SetLineStyle(sL);
  pLabel -> SetTextFont(txt);
  pLabel -> SetTextAlign(align);
  pLabel -> AddText(sLabel1.Data());
  pLabel -> AddText(sLabel2.Data());
  pLabel -> AddText(sLabel3.Data());
  pLabel -> AddText(sLabel4.Data());
  cout << "    Label created." << endl;


  // calculate chi2
  TPaveText *pChi2;

  const Bool_t  calculateChi2 = true;
  const Float_t calcRange[2]  = {0.2, 20.};
  if (calculateChi2) {

    // determine where to start and stop
    const Int_t   minA = hInputA -> FindFirstBinAbove(0.);
    const Int_t   maxA = hInputA -> FindLastBinAbove(0.);
    const Int_t   minB = hInputB -> FindFirstBinAbove(0.);
    const Int_t   maxB = hInputB -> FindLastBinAbove(0.);
    const Int_t   iMin = TMath::Max(minA, minB);
    const Int_t   iMax = TMath::Min(maxA, maxB);
    const Float_t xMin = hInputA -> GetBinCenter(iMin);
    const Float_t xMax = hInputA -> GetBinCenter(iMax);

    // loop over bins
    UInt_t   nChi = 0;
    Double_t sum2 = 0.;
    Double_t chi2 = 0.;
    for (UInt_t iBinA = 1; iBinA < (nBins + 1); iBinA++) {

      const Float_t xA = hInputA -> GetBinCenter(iBinA);
      const Float_t yA = hInputA -> GetBinContent(iBinA);
      const Float_t eA = hInputA -> GetBinError(iBinA);
      if ((xA < xMin) || (xA < calcRange[0])) continue;
      if ((xA > xMax) || (xA > calcRange[1])) continue;

      const UInt_t  iBinB = hInputB -> FindBin(xA);
      const Float_t xB    = hInputB -> GetBinCenter(iBinB);
      const Float_t yB    = hInputB -> GetBinContent(iBinB);
      const Float_t eB    = hInputB -> GetBinError(iBinB);
      if ((yA <= 0.) || (yB <= 0.)) continue;
      if ((eA <= 0.) || (eB <= 0.)) continue;

      const Double_t dErr = TMath::Sqrt((eA * eA) + (eB * eB));
      const Double_t num  = TMath::Power(yA - yB, 2.);
      const Double_t den  = TMath::Power(dErr, 2.);
      const Double_t val  = num / den;

      // calculate chi2
      sum2 += val;
      nChi++;

    }  // end bin loop

    if (nChi > 0) {
      sum2 /= (Double_t) nChi;
      chi2 = sum2;
    }

    // make text
    TString sChi2("#chi^{2} = ");
    sChi2 += chi2;

    const Float_t xyChi[4] = {0.7, 0.7, 0.9, 0.9};
    pChi2 = new TPaveText(xyChi[0], xyChi[1], xyChi[2], xyChi[3], "NDC NB");
    pChi2 -> SetFillColor(cL);
    pChi2 -> SetFillStyle(sL);
    pChi2 -> SetLineColor(cL);
    pChi2 -> SetLineStyle(sL);
    pChi2 -> SetTextFont(txt);
    pChi2 -> SetTextAlign(align);
    pChi2 -> AddText(sChi2.Data());
    cout << "    Calculate chi2:\n"
         << "      chi2 = " << chi2
         << endl;

  }  // end chi2 calculation


  // make plots
  fOutput -> cd();

  const Int_t   wC      = 800;
  const Int_t   hC      = 800;
  const Int_t   grd     = 0;
  const Int_t   log     = 0;
  const Float_t mar     = 0.;
  const Float_t pad1[4] = {0., 0., 1., 0.35};
  const Float_t pad2[4] = {0., 0.35, 1., 1.};
  TCanvas *cPlot = new TCanvas(sCanvas.Data(), "", wC, hC);
  TPad    *pPad1 = new TPad(sPad[0].Data(), "", pad1[0], pad1[1], pad1[2], pad1[3]);
  TPad    *pPad2 = new TPad(sPad[1].Data(), "", pad2[0], pad2[1], pad2[2], pad2[3]);
  pPad1   -> SetGrid(grd, grd);
  pPad1   -> SetTopMargin(mar);
  pPad2   -> SetGrid(grd, grd);
  pPad2   -> SetLogy(log);
  pPad2   -> SetBottomMargin(mar);
  cPlot   -> cd();
  pPad1   -> Draw();
  pPad2   -> Draw();
  pPad1   -> cd();
  hRatio  -> Draw();
  if (calculateChi2)
    pChi2 -> Draw();
  pPad2   -> cd();
  hInputA -> Draw();
  hInputB -> Draw("same");
  lLegend -> Draw();
  pLabel  -> Draw();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Plot drawn." << endl;


  fOutput -> cd();
  hInputA -> Write();
  hInputB -> Write();
  hRatio  -> Write();
  fOutput -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputA -> cd();
  fInputA -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
