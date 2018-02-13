// 'PlotRatios.C'
// 
// Use this to plot some ratios.

#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


void PlotRatios() {

  gStyle -> SetOptStat(0);

  // parameters
  const Bool_t   doEtaNorm = false;
  const Double_t rJet      = 0.3;
  const Double_t gPurityPP = 0.4;
  const Double_t pPurityPP = 0.6;
  const Double_t gPurityAA = 0.7;
  const Double_t pPurityAA = 0.3;


  // open files
  TFile *fGam = new TFile("pythiaCheckG.r03a02rm1full.d27m1y2017.root", "read");
  TFile *fPi0 = new TFile("pythiaCheckP.r03a02rm1full.d27m1y2017.root", "read");
  TFile *fOut = new TFile("pythiaCheckR.r03a02rm1full.d28m1y2017.root", "recreate");

  // grab and make distributions
  const Double_t hNorm = 2. * (1. - rJet);
  TH1D *hDir   = (TH1D*) fGam -> Get("QA/hPtRE");
  TH1D *hPi0   = (TH1D*) fPi0 -> Get("QA/hPtRE");
  TH1D *hRchPP = (TH1D*) hDir -> Clone();
  TH1D *hRchAA = (TH1D*) hDir -> Clone();
  TH1D *hRd    = (TH1D*) hDir -> Clone();
  TH1D *hRrPP  = (TH1D*) hDir -> Clone();
  TH1D *hRrAA  = (TH1D*) hDir -> Clone();
  hDir   -> SetNameTitle("hGammaDirect", "Recoil jet p_{T}^{corr}, #gamma_{dir}");
  if (doEtaNorm) hDir -> Scale(1. / hNorm);
  hPi0   -> SetNameTitle("hPi0", "Recoil jet p_{T}^{corr}, #pi^{0}");
  if (doEtaNorm) hPi0 -> Scale(1. / hNorm);
  hRchPP -> SetNameTitle("hGammaRichPP", "Recoil jet p_{T}^{corr}, #gamma_{rich} (pp)");
  hRchPP -> Add(hPi0, hDir, pPurityPP, gPurityPP);
  hRchAA -> SetNameTitle("hGammaRichAA", "Recoil jet p_{T}^{corr}, #gamma_{rich} (AA)");
  hRchAA -> Add(hPi0, hDir, pPurityAA, gPurityAA);
  hRd    -> SetNameTitle("hRatioDirect", "");
  hRd    -> Divide(hPi0, hDir, 1, 1);
  hRrPP  -> SetNameTitle("hRatioRichPP", "");
  hRrPP  -> Divide(hPi0, hRchPP, 1, 1);
  hRrAA  -> SetNameTitle("hRatioRichAA", "");
  hRrAA  -> Divide(hPi0, hRchAA, 1, 1);


  // set styles
  hDir   -> GetXaxis() -> SetRangeUser(-2., 30.);
  hDir   -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{jet}/dp_{T}^{corr}d#eta^{jet}");
  hDir   -> GetYaxis() -> CenterTitle();
  hDir   -> SetMarkerStyle(4);
  hDir   -> SetMarkerColor(kRed);
  hDir   -> SetLineColor(kRed);
  hPi0   -> GetXaxis() -> SetRangeUser(-2., 30.);
  hPi0   -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{jet}/dp_{T}^{corr}d#eta^{jet}");
  hPi0   -> GetYaxis() -> CenterTitle();
  hPi0   -> SetMarkerStyle(29);
  hPi0   -> SetMarkerColor(kBlue);
  hPi0   -> SetLineColor(kBlue);
  hRchPP -> GetXaxis() -> SetRangeUser(-2., 30.);
  hRchPP -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{jet}/dp_{T}^{corr}d#eta^{jet}");
  hRchPP -> GetYaxis() -> CenterTitle();
  hRchPP -> SetMarkerStyle(27);
  hRchPP -> SetMarkerColor(kGreen);
  hRchPP -> SetLineColor(kGreen);
  hRchAA -> GetXaxis() -> SetRangeUser(-2., 30.);
  hRchAA -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{jet}/dp_{T}^{corr}d#eta^{jet}");
  hRchAA -> GetYaxis() -> CenterTitle();
  hRchAA -> SetMarkerStyle(27);
  hRchAA -> SetMarkerColor(kGreen);
  hRchAA -> SetLineColor(kGreen);
  hRd    -> GetXaxis() -> SetRangeUser(-2., 30.);
  hRd    -> GetYaxis() -> SetRangeUser(0., 10.);
  hRd    -> GetXaxis() -> SetTitle("p_{T}^{corr} = p_{T}^{jet} - #rhoA");
  hRd    -> GetXaxis() -> CenterTitle();
  hRd    -> GetYaxis() -> SetTitle("#pi^{0}/#gamma");
  hRd    -> GetYaxis() -> CenterTitle();
  hRd    -> SetMarkerStyle(8);
  hRd    -> SetMarkerColor(kViolet);
  hRd    -> SetLineColor(kViolet);
  hRrPP  -> GetXaxis() -> SetRangeUser(-2., 30.);
  hRrPP  -> GetYaxis() -> SetRangeUser(0., 10.);
  hRrPP  -> GetXaxis() -> SetTitle("p_{T}^{corr} = p_{T}^{jet} - #rhoA");
  hRrPP  -> GetXaxis() -> CenterTitle();
  hRrPP  -> GetYaxis() -> SetTitle("#pi^{0}/#gamma");
  hRrPP  -> GetYaxis() -> CenterTitle();
  hRrPP  -> SetMarkerStyle(33);
  hRrPP  -> SetMarkerColor(kTeal);
  hRrPP  -> SetLineColor(kTeal);
  hRrAA  -> GetXaxis() -> SetRangeUser(-2., 30.);
  hRrAA  -> GetYaxis() -> SetRangeUser(0., 10.);
  hRrAA  -> GetXaxis() -> SetTitle("p_{T}^{corr} = p_{T}^{jet} - #rhoA");
  hRrAA  -> GetXaxis() -> CenterTitle();
  hRrAA  -> GetYaxis() -> SetTitle("#pi^{0}/#gamma");
  hRrAA  -> GetYaxis() -> CenterTitle();
  hRrAA  -> SetMarkerStyle(33);
  hRrAA  -> SetMarkerColor(kTeal);
  hRrAA  -> SetLineColor(kTeal);


  // make labels
  TLegend *lg1 = new TLegend(0.5, 0.5, 0.7, 0.7);
  lg1 -> SetLineColor(kWhite);
  lg1 -> SetFillColor(kWhite);
  lg1 -> AddEntry(hDir, "#gamma^{dir} trigger");
  lg1 -> AddEntry(hRchPP, "#gamma^{rich} trigger (pp)");
  lg1 -> AddEntry(hPi0, "#pi^{0} trigger");

  TLegend *lg2 = new TLegend(0.5, 0.5, 0.7, 0.7);
  lg2 -> SetLineColor(kWhite);
  lg2 -> SetFillColor(kWhite);
  lg2 -> AddEntry(hRd, "#pi^{0}/#gamma^{dir}");
  lg2 -> AddEntry(hRrPP, "#pi^{0}/#gamma^{rich}(pp)");

  TLegend *lg3 = new TLegend(0.5, 0.5, 0.7, 0.7);
  lg3 -> SetLineColor(kWhite);
  lg3 -> SetFillColor(kWhite);
  lg3 -> AddEntry(hDir, "#gamma^{dir} trigger");
  lg3 -> AddEntry(hRchAA, "#gamma^{rich} trigger (AA)");
  lg3 -> AddEntry(hPi0, "#pi^{0} trigger");

  TLegend *lg4 = new TLegend(0.5, 0.5, 0.7, 0.7);
  lg4 -> SetLineColor(kWhite);
  lg4 -> SetFillColor(kWhite);
  lg4 -> AddEntry(hRd, "#pi^{0}/#gamma^{dir}");
  lg4 -> AddEntry(hRrAA, "#pi^{0}/#gamma^{rich}(AA)");

  TPaveText *pt = new TPaveText(0.5, 0.3, 0.7, 0.5, "NDC NB");
  pt -> SetLineColor(kWhite);
  pt -> SetFillColor(kWhite);
  pt -> AddText("Pythia, #sqrt{s} = 200 GeV");
  pt -> AddText("Recoil jets");
  pt -> AddText("anti-k_{T}, R = 0.5");
  pt -> AddText("A_{jet} > 0.65, N_{rm} = 1");


  // draw plots
  TCanvas *cRatPP = new TCanvas("cRatioPP", "Ratio of pi0-triggered to gamma-triggered (pp) jet distributions", 500, 500);
  TPad    *pRatPP = new TPad("pRatioPP", "Ratio (pp)", 0, 0, 1, 0.3);
  TPad    *pDisPP = new TPad("pDistributionPP", "Jet distributions (pp)", 0, 0.3, 1, 1);
  pRatPP -> SetFillStyle(4000);
  pRatPP -> SetTopMargin(0);
  pRatPP -> SetBottomMargin(0.2);
  pDisPP -> SetFillStyle(4000);
  pDisPP -> SetBottomMargin(0);
  pRatPP -> Draw();
  pDisPP -> Draw();
  pRatPP -> cd();
  pRatPP -> SetGrid(0, 0);
  hRd    -> Draw();
  hRrPP  -> Draw("same");
  lg2    -> Draw();
  pDisPP -> cd();
  pDisPP -> SetLogy(1);
  pDisPP -> SetGrid(0, 0);
  hDir   -> Draw();
  hRchPP -> Draw("same");
  hPi0   -> Draw("same");
  lg1    -> Draw();
  pt     -> Draw();
  cRatPP -> cd();
  cRatPP -> Write();
  cRatPP -> Close();

  TCanvas *cRatAA = new TCanvas("cRatioAA", "Ratio of pi0-triggered to gamma-triggered (AA) jet distributions", 500, 500);
  TPad    *pRatAA = new TPad("pRatioAA", "Ratio", 0, 0, 1, 0.3);
  TPad    *pDisAA = new TPad("pDistributionAA", "Jet distributions", 0, 0.3, 1, 1);
  pRatAA -> SetFillStyle(4000);
  pRatAA -> SetTopMargin(0);
  pRatAA -> SetBottomMargin(0.2);
  pDisAA -> SetFillStyle(4000);
  pDisAA -> SetBottomMargin(0);
  pRatAA -> Draw();
  pDisAA -> Draw();
  pRatAA -> cd();
  pRatAA -> SetGrid(0, 0);
  hRd    -> Draw();
  hRrAA  -> Draw("same");
  lg4    -> Draw();
  pDisAA -> cd();
  pDisAA -> SetLogy(1);
  pDisAA -> SetGrid(0, 0);
  hDir   -> Draw();
  hRchAA -> Draw("same");
  hPi0   -> Draw("same");
  lg3    -> Draw();
  pt     -> Draw();
  cRatAA -> cd();
  cRatAA -> Write();
  cRatAA -> Close();


  // save and close files
  fOut   -> cd();
  hDir   -> Write();
  hRchPP -> Write();
  hRchAA -> Write();
  hPi0   -> Write();
  hRd    -> Write();
  hRrPP  -> Write();
  hRrAA  -> Write();
  fOut   -> Close();

  fGam -> cd();
  fGam -> Close();
  fPi0 -> cd();
  fPi0 -> Close();

}

// End ------------------------------------------------------------------------
