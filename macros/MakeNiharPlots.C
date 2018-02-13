// 'MakeNiharPlots.C'
//
// A quick macro to make plots similar to Nihar's.

#include <TSystem>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


void MakeNiharPlots() {

  const Double_t xMin = -3.;
  const Double_t xMax = 40.;
  const Double_t yMin = 0.00001;
  const Double_t yMax = 5.;
  const Double_t cWF  = 700.;
  const Double_t cHF  = 500.;
  const Double_t cWP  = 500.;
  const Double_t cHP  = 500.;
  const TString  fTitle("#gamma_{dir} triggered jets; #Delta#varphi; (1/N_{Trig})dN^{jet}/d#Delta#varphi");
  const TString  pTitle("#gamma_{dir} triggered recoil jets; Recoil-Jet p_{T}^{Rec}-#rhoA [GeV/c]; (1/N_{Trig})dN^{jet}/dp_{T}");
  const TString  beam("Pythia:p+p 200 GeV");
  const TString  trig("Trigger-#gamma_{dir}");
  const TString  jets("anti-k_{T}, R=0.3");
  const TString  eTrg("9<p_{T}^{trig}<20 GeV/c");

  TFile *oFile   = new TFile("Pythia22g.r03a02rm1.gPlots.Aug30.root", "recreate");
  TFile *iFile   = new TFile("Pythia22p.r03a02rm1.g.Aug30.root");
  TH1D  *hDfF    = (TH1D*) iFile -> Get("QA/hDfJet_full");
  TH1D  *hDfC    = (TH1D*) iFile -> Get("QA/hDfJet_chrg");
  TH1D  *hPtF    = (TH1D*) iFile -> Get("QA/hPtRE_full");
  TH1D  *hPtC    = (TH1D*) iFile -> Get("QA/hPtRE_chrg");
  TH1D  *hDfFull = hDfF -> Clone("hDfJet_full");
  TH1D  *hDfChrg = hDfC -> Clone("hDfJet_chrg");
  TH1D  *hPtFull = hPtF -> Clone("hPtCorr_full");
  TH1D  *hPtChrg = hPtC -> Clone("hPtCorr_chrg");

  // compute dF residuals
  const Int_t    nDf = hDfFull -> GetNbinsX();
  const Double_t dF1 = hDfFull -> GetBinLowEdge(1);
  const Double_t dF2 = hDfFull -> GetBinLowEdge(nDf+1);
  TH1D  *hDfRes      = new TH1D("hDfRes", "#Delta#varphi residuals [#Delta#varphi(charged)/#Delta#varphi(full)", nDf, dF1, dF2);
  hDfRes -> Divide(hDfChrg, hDfFull, 1., 1.);
  hDfRes -> SetMarkerColor(kTeal+10);
  hDfRes -> SetLineColor(kTeal+10);

  // create label
  TPaveText *pt = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  pt -> AddText(beam);
  pt -> AddText(trig);
  pt -> AddText(jets);
  pt -> AddText(eTrg);
  pt -> SetFillColor(kWhite);


  // full dF style
  hDfFull -> SetMarkerColor(kOrange+10);
  hDfFull -> SetLineColor(kOrange+10);
  hDfFull -> SetTitle(fTitle);
  // charged dF style
  hDfChrg -> SetLineColor(kViolet+10);
  hDfChrg -> SetLineColor(kViolet+10);
  hDfChrg -> SetTitle(fTitle);
  oFile   -> cd();

  // create dF legend
  TLegend *lgf = new TLegend(0.1, 0.1, 0.3, 0.3);
  lgf -> AddEntry(hDfFull, "Full Jet");
  lgf -> AddEntry(hDfChrg, "Charged Jet");
  lgf -> SetLineColor(kWhite);
  lgf -> SetFillColor(kWhite);

  TCanvas *cf  = new TCanvas("cFullVsChrgDf", "Full vs. Charged #Delta#varphi", cWF, cHF);
  cf      -> SetCanvasSize(cWF, cHF);
  TPad    *pF1 = new TPad("pF1", "Residuals", 0, 0, 1, 0.3);
  TPad    *pF2 = new TPad("pF2", "Data", 0, 0.3, 1, 1);
  pF1     -> SetFillStyle(4000);
  pF2     -> SetFillStyle(4000);
  pF1     -> Draw();
  pF2     -> Draw();
  pF1     -> cd();
  pF1     -> SetGrid(0, 0);
  hDfRes  -> Draw();
  pF2     -> cd();
  pF2     -> SetGrid(0, 0);
  hDfFull -> Draw();
  hDfChrg -> Draw("same");
  lgf     -> Draw();
  pt      -> Draw();
  cf      -> Write();
  cf      -> Close();


  // full pT style
  hPtFull -> SetMarkerStyle(8);
  hPtFull -> SetMarkerColor(kGreen+2);
  hPtFull -> SetLineColor(kGreen+2);
  hPtFull -> GetXaxis() -> SetRangeUser(xMin, xMax);
  hPtFull -> GetYaxis() -> SetRangeUser(yMin, yMax);
  hPtFull -> SetTitle(pTitle);
  hPtFull -> SetName("hPtCorr_full");
  // charged pT style
  hPtChrg -> SetMarkerStyle(33);
  hPtChrg -> SetMarkerColor(kMagenta+2);
  hPtChrg -> SetLineColor(kMagenta+2);
  hPtChrg -> GetXaxis() -> SetRangeUser(xMin, xMax);
  hPtChrg -> GetYaxis() -> SetRangeUser(yMin, yMax);
  hPtChrg -> SetTitle(pTitle);
  hPtChrg -> SetName("hPtCorr_chrg");

  // create pT legend
  TLegend *lgp = new TLegend(0.1, 0.1, 0.3, 0.3);
  lgp -> AddEntry(hPtFull, "Full Jet");
  lgp -> AddEntry(hPtChrg, "Charged Jet");
  lgp -> SetLineColor(kWhite);
  lgp -> SetFillColor(kWhite);
  
  TCanvas *cp  = new TCanvas("cFullVsChrgPt", "Full vs. Charged p_{T}^{corr}", cWP, cHP);
  cp      -> SetCanvasSize(cWP, cHP);
  cp      -> SetGrid(0, 0);
  cp      -> SetLogy();
  hPtFull -> Draw();
  hPtChrg -> Draw("same");
  lgp     -> Draw();
  pt      -> Draw();
  cp      -> Write();
  cp      -> Close();


  hDfFull -> Write();
  hDfChrg -> Write();
  hDfRes  -> Write();
  hPtFull -> Write();
  hPtChrg -> Write();
  oFile   -> Close();
  iFile   -> cd();
  iFile   -> Close();

}

// End ------------------------------------------------------------------------
