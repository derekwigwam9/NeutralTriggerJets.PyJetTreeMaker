// 'CompareNineThings.C'
// Derek Anderson
// 05.25.2017
//
// Use this to compare nine histograms.


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
static const TString sOutput("dataXembedding.eTtrg.pi0vsHighTwrs.d3m9y2017.root");
static const TString sInputA("../JetData/pp200r12pt5.highTwr.d3m9y2017.root");
static const TString sInputB("../JetData/pp200r12pt7.highTwr.d3m9y2017.root");
static const TString sInputC("../JetData/pp200r12pt9.highTwr.d3m9y2017.root");
static const TString sInputD("../JetData/pp200r12pt11.highTwr.d3m9y2017.root");
static const TString sInputE("../JetData/pp200r12pt15.highTwr.d3m9y2017.root");
static const TString sInputF("../JetData/pp200r12pt20.highTwr.d3m9y2017.root");
static const TString sInputG("../JetData/pp200r12pt25.highTwr.d3m9y2017.root");
static const TString sInputH("../JetData/pp200r12pt35.highTwr.d3m9y2017.root");
static const TString sInputI("../JetData/pp200r9.pi0.d3m9y2017.root");
static const TString sHistA("hTrgEt");
static const TString sHistB("hTrgEt");
static const TString sHistC("hTrgEt");
static const TString sHistD("hTrgEt");
static const TString sHistE("hTrgEt");
static const TString sHistF("hTrgEt");
static const TString sHistG("hTrgEt");
static const TString sHistH("hTrgEt");
static const TString sHistI("hTrgEt");
// histogram parameters
static const TString sNameA("hTrgEt_pt57");
static const TString sNameB("hTrgEt_pt79");
static const TString sNameC("hTrgEt_pt911");
static const TString sNameD("hTrgEt_pt1115");
static const TString sNameE("hTrgEt_pt1520");
static const TString sNameF("hTrgEt_pt2025");
static const TString sNameG("hTrgEt_pt2535");
static const TString sNameH("hTrgEt_pt35-1");
static const TString sNameI("hTrgEt_data");
static const TString sTitle("Triger energy");
static const TString sTitleX("E^{trg}_{T}");
static const TString sTitleY("dN_{trg}/dE^{trg}_{T}");
// legend parameters
static const TString sLegendA("p_{T}^{part}#in(5, 7) GeV/c");
static const TString sLegendB("p_{T}^{part}#in(7, 9) GeV/c");
static const TString sLegendC("p_{T}^{part}#in(9, 11) GeV/c");
static const TString sLegendD("p_{T}^{part}#in(11, 15) GeV/c");
static const TString sLegendE("p_{T}^{part}#in(15, 20) GeV/c");
static const TString sLegendF("p_{T}^{part}#in(20, 25) GeV/c");
static const TString sLegendG("p_{T}^{part}#in(25, 35) GeV/c");
static const TString sLegendH("p_{T}^{part}#in(35, -1) GeV/c");
static const TString sLegendI("Run 9 data");
// label parameters
static const TString sLabel1("pp-collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} > 6 GeV, |#eta^{trg}| < 1.0");
static const TString sLabel3("Hot towers removed");
static const TString sLabel4("Anti-k_{T}, R = 0.3, charged jet");
// canvas parameters
static const TString sCanvasH("cTrgEt");

// rebin histograms if necesary
static const Int_t  rebinA  = 10;
static const Int_t  rebinB  = 10;
static const Int_t  rebinC  = 10;
static const Int_t  rebinD  = 10;
static const Int_t  rebinE  = 10;
static const Int_t  rebinF  = 10;
static const Int_t  rebinG  = 10;
static const Int_t  rebinH  = 10;
static const Int_t  rebinI  = 10;
static const Bool_t doRebin = false;
// scale histograms if necesary
static const Double_t scaleA1 = 17.;
static const Double_t scaleB1 = 178.;
static const Double_t scaleC1 = 799.;
static const Double_t scaleD1 = 1574.;
static const Double_t scaleE1 = 3796.;
static const Double_t scaleF1 = 8229.;
static const Double_t scaleG1 = 16952.;
static const Double_t scaleH1 = 10421.;
static const Double_t scaleI1 = 871483.;
static const Double_t scaleA2 = 1.;
static const Double_t scaleB2 = 1.;
static const Double_t scaleC2 = 1.;
static const Double_t scaleD2 = 1.;
static const Double_t scaleE2 = 1.;
static const Double_t scaleF2 = 1.;
static const Double_t scaleG2 = 1.;
static const Double_t scaleH2 = 1.;
static const Double_t scaleI2 = 1.;
static const Bool_t   doScale = false;
static const Bool_t   doInt   = false;



void CompareNineThings() {

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
  TFile *fInputH = new TFile(sInputH.Data(), "read");
  TFile *fInputI = new TFile(sInputI.Data(), "Read");
  if (!fOutput || !fInputA || !fInputB || !fInputC || !fInputD || !fInputE || !fInputF || !fInputG || !fInputH || !fInputI) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
    assert(fInputC);
    assert(fInputD);
    assert(fInputE);
    assert(fInputF);
    assert(fInputG);
    assert(fInputH);
    assert(fInputI);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  TH1D *hInputC = (TH1D*) fInputC -> Get(sHistC.Data());
  TH1D *hInputD = (TH1D*) fInputD -> Get(sHistD.Data());
  TH1D *hInputE = (TH1D*) fInputE -> Get(sHistE.Data());
  TH1D *hInputF = (TH1D*) fInputF -> Get(sHistF.Data());
  TH1D *hInputG = (TH1D*) fInputG -> Get(sHistG.Data());
  TH1D *hInputH = (TH1D*) fInputH -> Get(sHistH.Data());
  TH1D *hInputI = (TH1D*) fInputI -> Get(sHistI.Data());
  if (!hInputA || !hInputB || !hInputC || !hInputD || !hInputE || !hInputF || !hInputG || !hInputH || !hInputI) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
    assert(hInputC);
    assert(hInputD);
    assert(hInputE);
    assert(hInputF);
    assert(hInputG);
    assert(hInputH);
    assert(hInputI);
  }
  cout << "    Histograms grabbed." << endl;

  // rebin histograms
  if (doRebin) {
    if (rebinA != 1) hInputA -> Rebin(rebinA);
    if (rebinB != 1) hInputB -> Rebin(rebinB);
    if (rebinC != 1) hInputC -> Rebin(rebinC);
    if (rebinD != 1) hInputD -> Rebin(rebinD);
    if (rebinE != 1) hInputE -> Rebin(rebinE);
    if (rebinF != 1) hInputF -> Rebin(rebinF);
    if (rebinG != 1) hInputG -> Rebin(rebinG);
    if (rebinH != 1) hInputH -> Rebin(rebinH);
    if (rebinI != 1) hInputI -> Rebin(rebinI);
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
    hInputH -> Scale(1. / scaleH1);
    hInputH -> Scale(1. / scaleH2);
    hInputI -> Scale(1. / scaleI1);
    hInputI -> Scale(1. / scaleI2);
    if (doInt) {
      const Double_t iA = hInputA -> Integral();
      const Double_t iB = hInputB -> Integral();
      const Double_t iC = hInputC -> Integral();
      const Double_t iD = hInputD -> Integral();
      const Double_t iE = hInputE -> Integral();
      const Double_t iF = hInputF -> Integral();
      const Double_t iG = hInputG -> Integral();
      const Double_t iH = hInputH -> Integral();
      const Double_t iI = hInputI -> Integral();
      hInputA -> Scale(1. / iA);
      hInputB -> Scale(1. / iB);
      hInputC -> Scale(1. / iC);
      hInputD -> Scale(1. / iD);
      hInputE -> Scale(1. / iE);
      hInputF -> Scale(1. / iF);
      hInputG -> Scale(1. / iG);
      hInputH -> Scale(1. / iH);
      hInputI -> Scale(1. / iI);
    }
  }


  // set histogram styles
  const Int_t    cA  = 810;
  const Int_t    cB  = 800;
  const Int_t    cC  = 830;
  const Int_t    cD  = 850;
  const Int_t    cE  = 870;
  const Int_t    cF  = 860;
  const Int_t    cG  = 890;
  const Int_t    cH  = 910;
  const Int_t    cI  = 1;
  const Int_t    mA  = 1;
  const Int_t    mB  = 1;
  const Int_t    mC  = 1;
  const Int_t    mD  = 1;
  const Int_t    mE  = 1;
  const Int_t    mF  = 1;
  const Int_t    mG  = 1;
  const Int_t    mH  = 1;
  const Int_t    mI  = 1;
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
  hInputH -> SetLineColor(cH);
  hInputH -> SetLineColor(cH);
  hInputH -> SetMarkerColor(cH);
  hInputH -> SetMarkerStyle(mH);
  hInputH -> SetTitleFont(txt);
  hInputH -> SetTitle(sTitle.Data());
  hInputH -> SetName(sNameH.Data());
  hInputH -> GetXaxis() -> SetLabelSize(lab);
  hInputH -> GetXaxis() -> CenterTitle(cnt);
  hInputH -> GetXaxis() -> SetTitleFont(txt);
  hInputH -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputH -> GetYaxis() -> SetLabelSize(lab);
  hInputH -> GetYaxis() -> CenterTitle(cnt);
  hInputH -> GetYaxis() -> SetTitleFont(txt);
  hInputH -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputI -> SetLineColor(cI);
  hInputI -> SetLineColor(cI);
  hInputI -> SetMarkerColor(cI);
  hInputI -> SetMarkerStyle(mI);
  hInputI -> SetTitleFont(txt);
  hInputI -> SetTitle(sTitle.Data());
  hInputI -> SetName(sNameI.Data());
  hInputI -> GetXaxis() -> SetLabelSize(lab);
  hInputI -> GetXaxis() -> CenterTitle(cnt);
  hInputI -> GetXaxis() -> SetTitleFont(txt);
  hInputI -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputI -> GetYaxis() -> SetLabelSize(lab);
  hInputI -> GetYaxis() -> CenterTitle(cnt);
  hInputI -> GetYaxis() -> SetTitleFont(txt);
  hInputI -> GetYaxis() -> SetTitle(sTitleY.Data());
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
  lLegendH -> AddEntry(hInputH, sLegendH.Data());
  lLegendH -> AddEntry(hInputI, sLegendI.Data());
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
  hInputA  -> Draw("");
  hInputB  -> Draw("same");
  hInputC  -> Draw("same");
  hInputD  -> Draw("same");
  hInputE  -> Draw("same");
  hInputF  -> Draw("same");
  hInputG  -> Draw("same");
  hInputH  -> Draw("same");
  hInputI  -> Draw("same");
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
  hInputH  -> Write();
  hInputI  -> Write();
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
  fInputH  -> cd();
  fInputH  -> Close();
  fInputI  -> cd();
  fInputI  -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
