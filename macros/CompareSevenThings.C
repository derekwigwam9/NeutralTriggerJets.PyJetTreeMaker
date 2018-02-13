// 'CompareSevenThrings.C'
// Derek Anderson
// 06.07.2017
//
// Use this to compare seven histograms.


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
static const TString sOutput("pythiaXembedding.shapeCheck.pTtrk.d7m6y2017.root");
static const TString sInputA("input/pythia.particle.d22m5y2017.root");
static const TString sInputB("input/pp200r12.tupleHistograms.d2m6y2017.root");
static const TString sInputC("input/pp200r12.tupleHistograms.d2m6y2017.root");
static const TString sInputD("input/cut.gam.allEvts.d5m6y2017.root");
static const TString sInputE("input/cut.gam.allEvts.d5m6y2017.root");
static const TString sInputF("input/ptTrks2_pi0trig.root");
static const TString sInputG("input/pp200r9.after.r03rm1chrg.d15m5y2017.root");
static const TString sHistA("QA/hPtTrk");
static const TString sHistB("hPtPy");
static const TString sHistC("hPtGe");
static const TString sHistD("hPtTrk_NoCuts");
static const TString sHistE("FitCut2/hPtTrkG_n2d0");
static const TString sHistF("htrks2_pi0trig");
static const TString sHistG("hTrkPtP");
// histogram parameters
static const TString sNameA("hTrkPt_r12u_cut");
static const TString sNameB("hTrkPt_r12u_uncut");
static const TString sNameC("hTrkPt_r12g");
static const TString sNameD("hTrkPt_r12p");
static const TString sNameE("hTrkPt_Pyth");
static const TString sNameF("hTrkPt_DataHad");
static const TString sNameG("hTrkPt_DataJet");
static const TString sTitle("Primary track p_{T}");
static const TString sTitleX("p_{T}^{trk}");
static const TString sTitleY("(1/N_{trg}) dN_{trk}/dp_{T}^{trk}");
// legend parameters
static const TString sLegendA("Pythia, #pi^{0} trigger");
static const TString sLegendB("Run 12, Pythia");
static const TString sLegendC("Run 12, Geant");
static const TString sLegendD("Run 12, MuDst (before cuts)");
static const TString sLegendE("Run 12, MuDst (after cuts)");
static const TString sLegendF("Run 9, #gamma-h^{#pm}");
static const TString sLegendG("Run 9, #gamma-jet");
// label parameters
static const TString sLabel1("pp-collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9, 30) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} #in (0.2, 30) GeV/c, |#eta^{trk}| < 1.0");
static const TString sLabel4("Anti-k_{T}, R = 0.3, charged jet");
// canvas parameters
static const TString sCanvasH("hTrkPt");
static const TString sCanvasR("hRatio");

// rebin histograms if necesary
static const Int_t  rebinA  = 5;
static const Int_t  rebinB  = 50;
static const Int_t  rebinC  = 50;
static const Int_t  rebinD  = 50;
static const Int_t  rebinE  = 50;
static const Int_t  rebinF  = 1;
static const Int_t  rebinG  = 5;
static const Bool_t doRebin = true;

// scale histograms if necesary
static const Double_t scaleA1 = 1.;
static const Double_t scaleB1 = 4424872.;
static const Double_t scaleC1 = 4424872.;
static const Double_t scaleD1 = 4531611.;
static const Double_t scaleE1 = 4531611.;
static const Double_t scaleF1 = 1.;
static const Double_t scaleG1 = 1.;
static const Double_t scaleA2 = 5.;
static const Double_t scaleB2 = 0.5;
static const Double_t scaleC2 = 0.5;
static const Double_t scaleD2 = 0.5;
static const Double_t scaleE2 = 0.5;
static const Double_t scaleF2 = 1.;
static const Double_t scaleG2 = 5.;
static const Bool_t   doScale = true;
static const Bool_t   doInt   = true;



void CompareSevenThings() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  TFile *fInputC = new TFile(sInputC.Data(), "read");
  TFile *fInputD = new TFile(sInputD.Data(), "read");
  TFile *fInputE = new TFile(sInputE.Data(), "read");
  TFile *fInputF = new TFile(sInputF.Data(), "read");
  TFile *fInputG = new TFile(sInputG.Data(), "read");
  if (!fOutput || !fInputA || !fInputB || !fInputC || !fInputD || !fInputE || !fInputF || !fInputG) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
    assert(fInputC);
    assert(fInputD);
    assert(fInputE);
    assert(fInputF);
    assert(fInputG);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  TH1D *hInputC = (TH1D*) fInputC -> Get(sHistC.Data());
  TH1D *hInputD = (TH1D*) fInputD -> Get(sHistD.Data());
  TH1D *hInputE = (TH1D*) fInputE -> Get(sHistE.Data());
  TH1D *hInputF = (TH1D*) fInputF -> Get(sHistF.Data());
  TH1D *hInputG = (TH1D*) fInputG -> Get(sHistG.Data());
  if (!hInputA || !hInputB || !hInputC || !hInputD || !hInputE || !hInputF || !hInputG) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
    assert(hInputC);
    assert(hInputD);
    assert(hInputE);
    assert(hInputF);
    assert(hInputG);
  }
  cout << "    Histograms grabbed." << endl;

  // scale histograms
  if (doScale) {
    if (rebinA != 1) hInputA -> Rebin(rebinA);
    if (rebinB != 1) hInputB -> Rebin(rebinB);
    if (rebinC != 1) hInputC -> Rebin(rebinC);
    if (rebinD != 1) hInputD -> Rebin(rebinD);
    if (rebinE != 1) hInputE -> Rebin(rebinE);
    if (rebinF != 1) hInputF -> Rebin(rebinF);
    if (rebinG != 1) hInputG -> Rebin(rebinG);
  }


  // scale histograms
  if (doScale) {
    hInputA -> Scale(1. / scaleA1);
    hInputA -> Scale(1. / scaleA2);
    hInputB -> Scale(1. / scaleB1);
    hInputB -> Scale(1. / scaleB2);
    hInputC -> Scale(1. / scaleC1);
    hInputC -> Scale(1. / scaleC2);
    hInputD -> Scale(1. / scaleD1);
    hInputD -> Scale(1. / scaleD2);
    hInputE -> Scale(1. / scaleE1);
    hInputE -> Scale(1. / scaleE2);
    hInputF -> Scale(1. / scaleF1);
    hInputF -> Scale(1. / scaleF2);
    hInputG -> Scale(1. / scaleG1);
    hInputG -> Scale(1. / scaleG2);
    if (doInt) {
      const Double_t iA = hInputA -> Integral();
      const Double_t iB = hInputB -> Integral();
      const Double_t iC = hInputC -> Integral();
      const Double_t iD = hInputD -> Integral();
      const Double_t iE = hInputE -> Integral();
      const Double_t iF = hInputF -> Integral();
      const Double_t iG = hInputG -> Integral();
      hInputA -> Scale(1. / iA);
      hInputB -> Scale(1. / iB);
      hInputC -> Scale(1. / iC);
      hInputD -> Scale(1. / iD);
      hInputE -> Scale(1. / iE);
      hInputF -> Scale(1. / iF);
      hInputG -> Scale(1. / iG);
    }
  }


  // set histogram styles
  const Int_t    cA  = 810;
  const Int_t    cB  = 830;
  const Int_t    cC  = 850;
  const Int_t    cD  = 870;
  const Int_t    cE  = 890;
  const Int_t    cF  = 910;
  const Int_t    cG  = 1;
  const Int_t    mA  = 1;
  const Int_t    mB  = 1;
  const Int_t    mC  = 1;
  const Int_t    mD  = 1;
  const Int_t    mE  = 1;
  const Int_t    mF  = 1;
  const Int_t    mG  = 1;
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
  hInputE -> SetLineColor(cE);
  hInputE -> SetMarkerColor(cE);
  hInputE -> SetMarkerStyle(mE);
  hInputE -> SetTitleFont(txt);
  hInputE -> SetTitle(sTitle.Data());
  hInputE -> SetName(sNameE.Data());
  hInputE -> GetXaxis() -> SetLabelSize(lab);
  hInputE -> GetXaxis() -> CenterTitle(cnt);
  hInputE -> GetXaxis() -> SetTitleFont(txt);
  hInputE -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputE -> GetYaxis() -> SetLabelSize(lab);
  hInputE -> GetYaxis() -> CenterTitle(cnt);
  hInputE -> GetYaxis() -> SetTitleFont(txt);
  hInputE -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputF -> SetLineColor(cF);
  hInputF -> SetLineColor(cF);
  hInputF -> SetMarkerColor(cF);
  hInputF -> SetMarkerStyle(mF);
  hInputF -> SetTitleFont(txt);
  hInputF -> SetTitle(sTitle.Data());
  hInputF -> SetName(sNameF.Data());
  hInputF -> GetXaxis() -> SetLabelSize(lab);
  hInputF -> GetXaxis() -> CenterTitle(cnt);
  hInputF -> GetXaxis() -> SetTitleFont(txt);
  hInputF -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputF -> GetYaxis() -> SetLabelSize(lab);
  hInputF -> GetYaxis() -> CenterTitle(cnt);
  hInputF -> GetYaxis() -> SetTitleFont(txt);
  hInputF -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputG -> SetLineColor(cG);
  hInputG -> SetLineColor(cG);
  hInputG -> SetMarkerColor(cG);
  hInputG -> SetMarkerStyle(mG);
  hInputG -> SetTitleFont(txt);
  hInputG -> SetTitle(sTitle.Data());
  hInputG -> SetName(sNameG.Data());
  hInputG -> GetXaxis() -> SetLabelSize(lab);
  hInputG -> GetXaxis() -> CenterTitle(cnt);
  hInputG -> GetXaxis() -> SetTitleFont(txt);
  hInputG -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputG -> GetYaxis() -> SetLabelSize(lab);
  hInputG -> GetYaxis() -> CenterTitle(cnt);
  hInputG -> GetYaxis() -> SetTitleFont(txt);
  hInputG -> GetYaxis() -> SetTitle(sTitleY.Data());
  cout << "    Styles set." << endl;


  const Int_t    cL  = 0;
  const Int_t    fL  = 0;
  const Int_t    sL  = 0;
  const Double_t x1L = 0.1;
  const Double_t x2L = 0.3;
  const Double_t y1L = 0.1;
  const Double_t y2L = 0.3;
  TLegend *lLegendH = new TLegend(x1L, y1L, x2L, y2L);
  lLegendH -> SetFillColor(cL);
  lLegendH -> SetFillStyle(sL);
  lLegendH -> SetLineColor(cL);
  lLegendH -> SetLineStyle(sL);
  lLegendH -> SetTextFont(txt);
  lLegendH -> AddEntry(hInputA, sLegendA.Data());
  lLegendH -> AddEntry(hInputB, sLegendB.Data());
  lLegendH -> AddEntry(hInputC, sLegendC.Data());
  lLegendH -> AddEntry(hInputD, sLegendD.Data());
  lLegendH -> AddEntry(hInputE, sLegendE.Data());
  lLegendH -> AddEntry(hInputF, sLegendF.Data());
  lLegendH -> AddEntry(hInputG, sLegendG.Data());
  cout << "    Legends created." << endl;

  const Double_t x1H = 0.3;
  const Double_t x2H = 0.5;
  const Double_t y1H = 0.1;
  const Double_t y2H = 0.3;
  TPaveText *pLabelH = new TPaveText(x1H, y1H, x2H, y2H, "NDC NB");
  pLabelH -> SetFillColor(cL);
  pLabelH -> SetFillStyle(sL);
  pLabelH -> SetLineColor(cL);
  pLabelH -> SetLineStyle(sL);
  pLabelH -> SetTextFont(txt);
  pLabelH -> AddText(sLabel1.Data());
  pLabelH -> AddText(sLabel2.Data());
  pLabelH -> AddText(sLabel3.Data());
  cout << "    Labels created." << endl;


  // make plots
  fOutput -> cd();

  const Int_t wC  = 800;
  const Int_t hC  = 800;
  const Int_t grd = 0;
  const Int_t log = 0;
  TCanvas *cPlotH = new TCanvas(sCanvasH.Data(), "", wC, hC);
  cPlotH   -> SetGrid(grd, grd);
  cPlotH   -> SetLogy(log);
  hInputA  -> Draw();
  hInputB  -> Draw("same");
  hInputC  -> Draw("same");
  hInputD  -> Draw("same");
  hInputE  -> Draw("same");
  hInputF  -> Draw("same");
  hInputG  -> Draw("same");
  lLegendH -> Draw();
  pLabelH  -> Draw();
  cPlotH   -> Write();
  cPlotH   -> Close();


  fOutput  -> cd();
  hInputA  -> Write();
  hInputB  -> Write();
  hInputC  -> Write();
  hInputD  -> Write();
  hInputE  -> Write();
  hInputF  -> Write();
  hInputG  -> Write();
  fOutput  -> Close();
  fInputA  -> cd();
  fInputA  -> Close();
  fInputB  -> cd();
  fInputB  -> Close();
  fInputC  -> cd();
  fInputC  -> Close();
  fInputD  -> cd();
  fInputD  -> Close();
  fInputE  -> cd();
  fInputE  -> Close();
  fInputF  -> cd();
  fInputF  -> Close();
  fInputG  -> cd();
  fInputG  -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
