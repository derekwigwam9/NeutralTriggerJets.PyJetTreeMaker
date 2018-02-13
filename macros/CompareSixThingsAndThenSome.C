// 'CompareSixThringsAndThenSome.C'
// Derek Anderson
// 05.23.2017
//
// Use this to compare six histograms
// and calculate the ratio of some of
// them.


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
static const TString sOutput("pythiaXembedding.firstEvtsExcluded.pTtrk.d4m6y2017.root");
static const TString sInputA("cut.gam.firstEvtsExcluded.d30m5y2017.root");
static const TString sInputB("cut.gam.firstEvtsExcluded.d30m5y2017.root");
static const TString sInputC("pp200r12.tupleHistograms.d2m6y2017.root");
static const TString sInputD("pp200r12.tupleHistograms.d2m6y2017.root");
static const TString sInputE("pythia.particle.d22m5y2017.root");
static const TString sInputF("../JetData/pp200r9.after.r03rm1chrg.d15m5y2017.root");
static const TString sHistA("FitCut2/hPtTrkG_n2d0");
static const TString sHistB("hPtTrk_NoCuts");
static const TString sHistC("hPtGe");
static const TString sHistD("hPtPy");
static const TString sHistE("QA/hPtTrk");
static const TString sHistF("hTrkPtP");
// histogram parameters
static const TString sNameA("hTrkPt_r12u_cut");
static const TString sNameB("hTrkPt_r12u_uncut");
static const TString sNameC("hTrkPt_r12g");
static const TString sNameD("hTrkPt_r12p");
static const TString sNameE("hTrkPt_Pyth");
static const TString sNameF("hTrkPt_Data");
static const TString sTitle("Primary track p_{T}");
static const TString sTitleX("p_{T}^{trk}");
static const TString sTitleY("(1/N_{trg}) dN_{trk}/dp_{T}^{trk}");
// ratio parameters
static const TString sNameAB("hRatio_e90");
static const TString sNameAC("hRatio_e80");
static const TString sNameAD("hRatio_e70");
static const TString sNameAE("hRatio_e60");
static const TString sTitleR("Primary track p_{T}");
static const TString sTitleRX("p_{T}^{trk}");
static const TString sTitleRY("ratio");
// fit parameters
static const Int_t    nDec(2);
static const Int_t    xFit1(0);
static const Int_t    xFit2(30);
static const Bool_t   doFit(false);
static const Double_t yFitB(0.90);
static const Double_t yFitC(0.80);
static const Double_t yFitD(0.70);
static const Double_t yFitE(0.60);
static const TString  sNameFB("fit_e90");
static const TString  sNameFC("fit_e80");
static const TString  sNameFD("fit_e70");
static const TString  sNameFE("fit_e60");
static const TString  sFunc("[0]*(1-exp(-4*x))");
// legend parameters
static const TString sLegendA("Run 12, MuDst (after cuts)");
static const TString sLegendB("Run 12, MuDst (before cuts)");
static const TString sLegendC("Run 12, Geant");
static const TString sLegendD("Run 12, Pythia");
static const TString sLegendE("Pythia, #pi^{0} trigger");
static const TString sLegendF("Run 9, #pi^{0} trigger");
// label parameters
static const TString sLabel1("pp-collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9, 30) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} #in (0.2, 30) GeV/c, |#eta^{trk}| < 1.0");
static const TString sLabel4("Anti-k_{T}, R = 0.3, charged jet");
// canvas parameters
static const TString sCanvasH("hTrkPt");
static const TString sCanvasR("hRatio");

// rebin histograms if necesary
static const Int_t  rebinA  = 10;
static const Int_t  rebinB  = 10;
static const Int_t  rebinC  = 10;
static const Int_t  rebinD  = 10;
static const Int_t  rebinE  = 1;
static const Int_t  rebinF  = 1;
static const Bool_t doRebin = true;

// scale histograms if necesary
static const Int_t    scaleA1 = 4075061;
static const Int_t    scaleB1 = 4075061;
static const Int_t    scaleC1 = 4424872;
static const Int_t    scaleD1 = 4424872;
static const Int_t    scaleE1 = 1;
static const Int_t    scaleF1 = 1;
static const Double_t scaleA2 = 0.1;
static const Double_t scaleB2 = 0.1;
static const Double_t scaleC2 = 0.1;
static const Double_t scaleD2 = 0.1;
static const Double_t scaleE2 = 1.;
static const Double_t scaleF2 = 1.;
static const Bool_t   doScale = true;



void CompareSixThingsAndThenSome() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  TFile *fInputC = new TFile(sInputC.Data(), "read");
  TFile *fInputD = new TFile(sInputD.Data(), "read");
  TFile *fInputE = new TFile(sInputE.Data(), "read");
  TFile *fInputF = new TFile(sInputF.Data(), "read");
  if (!fOutput || !fInputA || !fInputB || !fInputC || !fInputD || !fInputE || !fInputF) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
    assert(fInputC);
    assert(fInputD);
    assert(fInputE);
    assert(fInputF);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  TH1D *hInputC = (TH1D*) fInputC -> Get(sHistC.Data());
  TH1D *hInputD = (TH1D*) fInputD -> Get(sHistD.Data());
  TH1D *hInputE = (TH1D*) fInputE -> Get(sHistE.Data());
  TH1D *hInputF = (TH1D*) fInputF -> Get(sHistF.Data());
  if (!hInputA || !hInputB || !hInputC || !hInputD || !hInputE || !hInputF) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
    assert(hInputC);
    assert(hInputD);
    assert(hInputE);
    assert(hInputF);
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
  }


  // scale histograms
  if (doScale) {
    hInputA -> Scale(1. / (Double_t) scaleA1);
    hInputA -> Scale(1. / scaleA2);
    hInputB -> Scale(1. / (Double_t) scaleB1);
    hInputB -> Scale(1. / scaleB2);
    hInputC -> Scale(1. / (Double_t) scaleC1);
    hInputC -> Scale(1. / scaleC2);
    hInputD -> Scale(1. / (Double_t) scaleD1);
    hInputD -> Scale(1. / scaleD2);
    hInputE -> Scale(1. / (Double_t) scaleE1);
    hInputE -> Scale(1. / scaleE2);
    hInputF -> Scale(1. / (Double_t) scaleF1);
    hInputF -> Scale(1. / scaleF2);
  }


  // calculate ratios
  const Int_t    nX = hInputA -> GetNbinsX();
  const Double_t x1 = hInputA -> GetBinLowEdge(1);
  const Double_t x2 = hInputA -> GetBinLowEdge(nX + 1);
  TH1D *hRatioAB = new TH1D(sNameAB.Data(), "", nX, x1, x2);
  TH1D *hRatioAC = new TH1D(sNameAC.Data(), "", nX, x1, x2);
  TH1D *hRatioAD = new TH1D(sNameAD.Data(), "", nX, x1, x2);
  TH1D *hRatioAE = new TH1D(sNameAE.Data(), "", nX, x1, x2);
  hRatioAB -> Divide(hInputB, hInputA, 1., 1.);
  hRatioAC -> Divide(hInputC, hInputA, 1., 1.);
  hRatioAD -> Divide(hInputD, hInputA, 1., 1.);
  //hRatioAE -> Divide(hInputE, hInputA, 1., 1.);


  const Int_t    cA  = 810;
  const Int_t    cB  = 830;
  const Int_t    cC  = 850;
  const Int_t    cD  = 870;
  const Int_t    cE  = 890;
  const Int_t    cF  = 1;
  const Int_t    mA  = 1;
  const Int_t    mB  = 1;
  const Int_t    mC  = 1;
  const Int_t    mD  = 1;
  const Int_t    mE  = 1;
  const Int_t    mF  = 1;
  const Int_t    txt = 42;
  const Int_t    cnt = 1;
  const Double_t lab = 0.02;
  // set histogram styles
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
  // set ratio styles
  hRatioAB -> SetLineColor(cB);
  hRatioAB -> SetLineColor(cB);
  hRatioAB -> SetMarkerColor(cB);
  hRatioAB -> SetMarkerStyle(mB);
  hRatioAB -> SetTitleFont(txt);
  hRatioAB -> SetTitle(sTitleR.Data());
  hRatioAB -> SetName(sNameAB.Data());
  hRatioAB -> GetXaxis() -> SetLabelSize(lab);
  hRatioAB -> GetXaxis() -> CenterTitle(cnt);
  hRatioAB -> GetXaxis() -> SetTitleFont(txt);
  hRatioAB -> GetXaxis() -> SetTitle(sTitleRX.Data());
  hRatioAB -> GetYaxis() -> SetLabelSize(lab);
  hRatioAB -> GetYaxis() -> CenterTitle(cnt);
  hRatioAB -> GetYaxis() -> SetTitleFont(txt);
  hRatioAB -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatioAC -> SetLineColor(cC);
  hRatioAC -> SetLineColor(cC);
  hRatioAC -> SetMarkerColor(cC);
  hRatioAC -> SetMarkerStyle(mC);
  hRatioAC -> SetTitleFont(txt);
  hRatioAC -> SetTitle(sTitleR.Data());
  hRatioAC -> SetName(sNameAC.Data());
  hRatioAC -> GetXaxis() -> SetLabelSize(lab);
  hRatioAC -> GetXaxis() -> CenterTitle(cnt);
  hRatioAC -> GetXaxis() -> SetTitleFont(txt);
  hRatioAC -> GetXaxis() -> SetTitle(sTitleRX.Data());
  hRatioAC -> GetYaxis() -> SetLabelSize(lab);
  hRatioAC -> GetYaxis() -> CenterTitle(cnt);
  hRatioAC -> GetYaxis() -> SetTitleFont(txt);
  hRatioAC -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatioAD -> SetLineColor(cD);
  hRatioAD -> SetLineColor(cD);
  hRatioAD -> SetMarkerColor(cD);
  hRatioAD -> SetMarkerStyle(mD);
  hRatioAD -> SetTitleFont(txt);
  hRatioAD -> SetTitle(sTitleR.Data());
  hRatioAD -> SetName(sNameAD.Data());
  hRatioAD -> GetXaxis() -> SetLabelSize(lab);
  hRatioAD -> GetXaxis() -> CenterTitle(cnt);
  hRatioAD -> GetXaxis() -> SetTitleFont(txt);
  hRatioAD -> GetXaxis() -> SetTitle(sTitleRX.Data());
  hRatioAD -> GetYaxis() -> SetLabelSize(lab);
  hRatioAD -> GetYaxis() -> CenterTitle(cnt);
  hRatioAD -> GetYaxis() -> SetTitleFont(txt);
  hRatioAD -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatioAE -> SetLineColor(cE);
  hRatioAE -> SetLineColor(cE);
  hRatioAE -> SetMarkerColor(cE);
  hRatioAE -> SetMarkerStyle(mE);
  hRatioAE -> SetTitleFont(txt);
  hRatioAE -> SetTitle(sTitleR.Data());
  hRatioAE -> SetName(sNameAE.Data());
  hRatioAE -> GetXaxis() -> SetLabelSize(lab);
  hRatioAE -> GetXaxis() -> CenterTitle(cnt);
  hRatioAE -> GetXaxis() -> SetTitleFont(txt);
  hRatioAE -> GetXaxis() -> SetTitle(sTitleRX.Data());
  hRatioAE -> GetYaxis() -> SetLabelSize(lab);
  hRatioAE -> GetYaxis() -> CenterTitle(cnt);
  hRatioAE -> GetYaxis() -> SetTitleFont(txt);
  hRatioAE -> GetYaxis() -> SetTitle(sTitleRY.Data());
  cout << "    Styles set." << endl;


  // fit ratios
  const Int_t cFB = 830;
  const Int_t cFC = 850;
  const Int_t cFD = 870;
  const Int_t cFE = 890;
  const Int_t sF  = 1;
  const Int_t wF  = 1;
  TF1 *fRatioAB = new TF1(sNameFB.Data(), sFunc.Data(), xFit1, xFit2);
  TF1 *fRatioAC = new TF1(sNameFC.Data(), sFunc.Data(), xFit1, xFit2);
  TF1 *fRatioAD = new TF1(sNameFD.Data(), sFunc.Data(), xFit1, xFit2);
  TF1 *fRatioAE = new TF1(sNameFE.Data(), sFunc.Data(), xFit1, xFit2);
  fRatioAB -> SetLineColor(cFB);
  fRatioAB -> SetLineStyle(sF);
  fRatioAB -> SetLineWidth(wF);
  fRatioAB -> SetParameter(0, yFitB);
  fRatioAC -> SetLineColor(cFC);
  fRatioAC -> SetLineStyle(sF);
  fRatioAC -> SetLineWidth(wF);
  fRatioAC -> SetParameter(0, yFitC);
  fRatioAD -> SetLineColor(cFD);
  fRatioAD -> SetLineStyle(sF);
  fRatioAD -> SetLineWidth(wF);
  fRatioAD -> SetParameter(0, yFitD);
  fRatioAE -> SetLineColor(cFE);
  fRatioAE -> SetLineStyle(sF);
  fRatioAE -> SetLineWidth(wF);
  fRatioAE -> SetParameter(0, yFitE);
  if (doFit) {
    hRatioAB -> Fit(fRatioAB, "QR");
    hRatioAC -> Fit(fRatioAC, "QR");
    hRatioAD -> Fit(fRatioAD, "QR");
    //hRatioAE -> Fit(fRatioAE, "QR");
    cout << "    Fits done." << endl;
  }

  TString sLabelAB("");
  TString sLabelAC("");
  TString sLabelAD("");
  TString sLabelAE("");
  if (doFit) {
    const Double_t eAB = fRatioAB -> GetParameter(0);
    const Double_t eAC = fRatioAC -> GetParameter(0);
    const Double_t eAD = fRatioAD -> GetParameter(0);
    const Double_t eAE = fRatioAE -> GetParameter(0);
    TString sRawB("");
    TString sRawC("");
    TString sRawD("");
    TString sRawE("");
    sRawB += eAB;
    sRawC += eAC;
    sRawD += eAD;
    sRawE += eAE;

    const Int_t nRawB = (Int_t) sRawB.First(".");
    const Int_t nRawC = (Int_t) sRawC.First(".");
    const Int_t nRawD = (Int_t) sRawD.First(".");
    const Int_t nRawE = (Int_t) sRawE.First(".");
    const Int_t nTxtB = nRawB + nDec + 1;
    const Int_t nTxtC = nRawC + nDec + 1;
    const Int_t nTxtD = nRawD + nDec + 1;
    const Int_t nTxtE = nRawE + nDec + 1;
    TString sTxtB("");
    TString sTxtC("");
    TString sTxtD("");
    TString sTxtE("");
    sTxtB.Append(sRawB, nTxtB);
    sTxtC.Append(sRawC, nTxtC);
    sTxtD.Append(sRawD, nTxtD);
    sTxtE.Append(sRawE, nTxtE);

    sLabelAB += "#color[";
    sLabelAB += cB;
    sLabelAB += "]{#epsilon_{fit} = ";
    sLabelAB += sTxtB.Data();
    sLabelAB += "}";
    sLabelAC += "#color[";
    sLabelAC += cC;
    sLabelAC += "]{#epsilon_{fit} = ";
    sLabelAC += sTxtC.Data();
    sLabelAC += "}";
    sLabelAD += "#color[";
    sLabelAD += cD;
    sLabelAD += "]{#epsilon_{fit} = ";
    sLabelAD += sTxtD.Data();
    sLabelAD += "}";
    sLabelAE += "#color[";
    sLabelAE += cE;
    sLabelAE += "]{#epsilon_{fit} = ";
    sLabelAE += sTxtE.Data();
    sLabelAE += "}";
    cout << "    Parameters extracted." << endl;
  }

  const Int_t    cL  = 0;
  const Int_t    fL  = 0;
  const Int_t    sL  = 0;
  const Double_t x1L = 0.1;
  const Double_t x2L = 0.3;
  const Double_t y1L = 0.1;
  const Double_t y2L = 0.3;
  TLegend *lLegendH = new TLegend(x1L, y1L, x2L, y2L);
  TLegend *lLegendR = new TLegend(x1L, y1L, x2L, y2L);
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
  lLegendR -> SetFillColor(cL);
  lLegendR -> SetFillStyle(sL);
  lLegendR -> SetLineColor(cL);
  lLegendR -> SetLineStyle(sL);
  lLegendR -> SetTextFont(txt);
  lLegendR -> AddEntry(hInputB, sLegendB.Data());
  lLegendR -> AddEntry(hInputC, sLegendC.Data());
  lLegendR -> AddEntry(hInputD, sLegendD.Data());
  //lLegendR -> AddEntry(hInputE, sLegendE.Data());
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

  const Double_t x1R = 0.5;
  const Double_t x2R = 0.7;
  const Double_t y1R = 0.1;
  const Double_t y2R = 0.3;
  TPaveText *pLabelR = new TPaveText(x1R, y1R, x2R, y2R, "NDC NB");
  pLabelR -> SetFillColor(cL);
  pLabelR -> SetFillStyle(sL);
  pLabelR -> SetLineColor(cL);
  pLabelR -> SetLineStyle(sL);
  pLabelR -> SetTextFont(txt);
  pLabelR -> AddText(sLabelAB.Data());
  pLabelR -> AddText(sLabelAC.Data());
  pLabelR -> AddText(sLabelAD.Data());
  pLabelR -> AddText(sLabelAE.Data());
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
  lLegendH -> Draw();
  pLabelH  -> Draw();
  cPlotH   -> Write();
  cPlotH   -> Close();

  TCanvas *cPlotR = new TCanvas(sCanvasR.Data(), "", wC, hC);
  cPlotR   -> SetGrid(grd, grd);
  cPlotR   -> SetLogy(log);
  hRatioAB -> Draw();
  hRatioAC -> Draw("same");
  hRatioAD -> Draw("same");
  //hRatioAE -> Draw("same");
  lLegendR -> Draw();
  pLabelH  -> Draw();
  if (doFit)
    pLabelR -> Draw();
  cPlotR   -> Write();
  cPlotR   -> Close();
  cout << "    Plot drawn." << endl;


  fOutput  -> cd();
  hInputA  -> Write();
  hInputB  -> Write();
  hInputC  -> Write();
  hInputD  -> Write();
  hInputE  -> Write();
  hInputF  -> Write();
  hRatioAB -> Write();
  hRatioAC -> Write();
  hRatioAD -> Write();
  //hRatioAE -> Write();
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
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
