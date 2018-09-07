// 'MakePrettyPlot.C'
// Derek Anderson
// 09.07.2018
//
// Use this to make a pretty plot comparing
// two jet distributions (and do the direct
// -photon if you want)
//
// NOTE: the gamma-subtraction assumes that
//       the histogram specified by 'sInA'
//       is the gamma-rich distribution.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const UInt_t NPad(2);
static const UInt_t NPlot(2);
static const UInt_t NHist(3);
static const UInt_t NRebin(2);
static const UInt_t NRebinVar(16);
static const Bool_t DoRebin(false);
static const Bool_t DoVariableRebin(false);
static const Bool_t DoGammaRichPlot(false);
static const Bool_t DoGammaSubtraction(false);



void MakePrettyPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distribution..." << endl;

  // io parameters
  const TString sOut("parVsDet.pythia.et920pi0.r05a065rm1chrg.d7m9y2018.root");
  const TString sInA("pp200py.resTestPlotDet.et920pi0.r05a065rm1chrg.d6m9y2018.root");
  const TString sInB("pp200py.resTestPlotPar.et920pi0.r05a065rm1chrg.d6m9y2018.root");
  const TString sHistA("Pi0/hJetPtCorrP");
  const TString sHistB("Pi0/hJetPtCorrP");

  // plot parameters
  const TString sTitle("");
  const TString sNameA("hParticle");
  const TString sNameB("hDetector");
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - #rhoA^{jet} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) dN^{jet}/d(p_{T}^{reco} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleR("ratio");
  const TString sLabelA("detector level");
  const TString sLabelB("particle level");

  // text parameters
  const TString sSys("pythia, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trigger, E_{T}^{trg} #in (9, 20) GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // subtraction parameters
  const UInt_t   fFilGR(3017);
  const UInt_t   fLinGR(1);
  const UInt_t   fWidGR(1);
  const UInt_t   fMarGR(4);
  const UInt_t   fColGR(879);
  const UInt_t   fColP[NHist] = {856, 921, 856};
  const UInt_t   fColG[NHist] = {896, 921, 896};
  const TString  sFileGam("pp200r9.binByBinTest.et911vz55gam.r02a005rm1chrg.p0m3k0n58t4.root");
  const TString  sFilePi0("pp200r9.binByBinTest.et911vz55pi0.r02a005rm1chrg.p0m3k0n58t4.root");
  const TString  sHistGam("hUnfolded");
  const TString  sHistPi0("hUnfolded");
  const TString  sNameGR("hGammaRich");
  const TString  sLabelGR("unfolded [#gamma^{rich}]");
  const Double_t gammaPurity(0.650);

  // misc parameters
  const Double_t plotRange[NPlot]    = {-3., 37.};
  const Double_t varRebin[NRebinVar] = {0., 1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 18., 21., 24., 27., 30.};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInA = new TFile(sInA.Data(), "read");
  TFile *fInB = new TFile(sInB.Data(), "read");
  if (!fOut || !fInA || !fInB) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInA = " << fInA << ", fInB" << fInB
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hHistA = (TH1D*) fInA -> Get(sHistA.Data());
  TH1D *hHistB = (TH1D*) fInB -> Get(sHistB.Data());
  if (!hHistA || !hHistB) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hHistA = " << hHistA << ", hHistB = " << hHistB
         << endl;
    return;
  }
  hHistA -> SetName(sNameA.Data());
  hHistB -> SetName(sNameB.Data());
  cout << "    Grabbed histograms." << endl;


  // do gamma-subtraction (if need be)
  if (DoGammaSubtraction) {
    TFile *fPi0 = new TFile(sFilePi0.Data(), "read");
    if (!fPi0) {
      cerr << "PANIC: couldn't open pi0 file!" << endl;
      return;
    }
    TH1D  *hPi0 = (TH1D*) fPi0    -> Get(sHistPi0.Data());
    TH1D  *hGam = (TH1D*) hHistA -> Clone();
    if (!hPi0 || !hGam) {
      cerr << "PANIC: couldn't grab pi0 or gamma histogram!\n"
           << "       hPi0 = " << hPi0 << ", hGam = " << hGam
           << endl;
      return;
    }
    hPi0   -> Scale(gammaPurity);
    hHistA -> Add(hGam, hPi0, 1., -1.);
    hHistA -> Scale(1. / (1. - gammaPurity));
    fPi0   -> Close();
    cout << "    Did gamma subtraction." << endl;
  }


  // rebin (if need be)
  hHistA -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  hHistB -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  if (DoRebin) {
    hHistA -> Rebin(NRebin);
    hHistB -> Rebin(NRebin);

    const UInt_t nBinA = hHistA -> GetNbinsX();
    const UInt_t nBinB = hHistB -> GetNbinsX();
    for (UInt_t iBinA = 1; iBinA <( nBinA + 1); iBinA++) {
      const Double_t valA = hHistA -> GetBinContent(iBinA);
      const Double_t errA = hHistA -> GetBinError(iBinA);
      const Double_t sizA = hHistA -> GetBinWidth(iBinA);
      hHistA -> SetBinContent(iBinA, valA / sizA);
      hHistA -> SetBinError(iBinA, errA / sizA);
    }
    for (UInt_t iBinB = 1; iBinB < (nBinB + 1); iBinB++) {
      const Double_t valB = hHistB -> GetBinContent(iBinB);
      const Double_t errB = hHistB -> GetBinError(iBinB);
      const Double_t sizB = hHistB -> GetBinWidth(iBinB);
      hHistB -> SetBinContent(iBinB, valB / sizB);
      hHistB -> SetBinError(iBinB, errB / sizB);
    }
    cout << "    Rebinned histograms (uniform)." << endl;
  }
  if (DoVariableRebin) {
    hHistA = (TH1D*) hHistA -> Rebin(NRebinVar - 1, sNameA.Data(), varRebin);
    hHistB = (TH1D*) hHistB -> Rebin(NRebinVar - 1, sNameB.Data(), varRebin);

    const UInt_t nBinA = hHistA -> GetNbinsX();
    const UInt_t nBinB = hHistB -> GetNbinsX();
    for (UInt_t iBinA = 1; iBinA < (nBinA + 1); iBinA++) {
      const Double_t valA = hHistA -> GetBinContent(iBinA);
      const Double_t errA = hHistA -> GetBinError(iBinA);
      const Double_t sizA = hHistA -> GetBinWidth(iBinA);
      hHistA -> SetBinContent(iBinA, valA / sizA);
      hHistA -> SetBinError(iBinA, errA / sizA);
    }
    for (UInt_t iBinB = 1; iBinB < (nBinB + 1); iBinB++) {
      const Double_t valB = hHistB -> GetBinContent(iBinB);
      const Double_t errB = hHistB -> GetBinError(iBinB);
      const Double_t sizB = hHistB -> GetBinWidth(iBinB);
      hHistB -> SetBinContent(iBinB, valB / sizB);
      hHistB -> SetBinError(iBinB, errB / sizB);
    }
    cout << "    Rebinned histograms (variable)." << endl;
  }

  // calculate ratio
  TH1D *hRatio = (TH1D*) hHistA -> Clone();
  hRatio -> SetName("hRatio");

  const UInt_t goodDiv = hRatio  -> Divide(hHistA, hHistB, 1., 1.);
  const UInt_t nBinR   = hHistA -> GetNbinsX();
  if (goodDiv == 0) {
    for (UInt_t iBinR = 0; iBinR < nBinR; iBinR++) {
      const UInt_t   iBinB = hHistB -> FindBin(hHistA -> GetBinCenter(iBinR));
      const Double_t num   = hHistA -> GetBinContent(iBinR);
      const Double_t den   = hHistB -> GetBinContent(iBinB);
      const Double_t nErr  = hHistA -> GetBinError(iBinR);
      const Double_t dErr  = hHistA -> GetBinError(iBinB);
      const Double_t nRel  = nErr / num;
      const Double_t dRel  = dErr / den;
      const Double_t ratio = num / den;
      const Double_t error = ratio * TMath::Sqrt((nRel * nRel) + (dRel * dRel));
      if (den > 0.) {
        hRatio -> SetBinContent(iBinR, ratio);
        hRatio -> SetBinError(iBinR, error);
      }
      else {
        hRatio -> SetBinContent(iBinR, 0.);
        hRatio -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
  }
  hRatio -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  cout << "    Calculated ratio." << endl;


  // set styles
  UInt_t fCol[NHist];
  if (DoGammaSubtraction) {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      fCol[iHist] = fColG[iHist];
    }
  }
  else {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      fCol[iHist] = fColP[iHist];
    }
  }

  const UInt_t  fMar[NHist] = {4, 8, 4};
  const UInt_t  fFil[NHist] = {0, 0, 0};
  const UInt_t  fLin[NHist] = {1, 1, 1};
  const UInt_t  fWid[NHist] = {1, 1, 1};
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hHistA -> SetMarkerColor(fCol[0]);
  hHistA -> SetMarkerStyle(fMar[0]);
  hHistA -> SetFillColor(fCol[0]);
  hHistA -> SetFillStyle(fFil[0]);
  hHistA -> SetLineColor(fCol[0]);
  hHistA -> SetLineStyle(fLin[0]);
  hHistA -> SetLineWidth(fWid[0]);
  hHistA -> SetTitle(sTitle.Data());
  hHistA -> SetTitleFont(fTxt);
  hHistA -> GetXaxis() -> SetTitle(sTitleX.Data());
  hHistA -> GetXaxis() -> SetTitleFont(fTxt);
  hHistA -> GetXaxis() -> SetTitleSize(fTit[1]);
  hHistA -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hHistA -> GetXaxis() -> SetLabelFont(fTxt);
  hHistA -> GetXaxis() -> SetLabelSize(fLab[1]);
  hHistA -> GetXaxis() -> CenterTitle(fCnt);
  hHistA -> GetYaxis() -> SetTitle(sTitleY.Data());
  hHistA -> GetYaxis() -> SetTitleFont(fTxt);
  hHistA -> GetYaxis() -> SetTitleSize(fTit[1]);
  hHistA -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hHistA -> GetYaxis() -> SetLabelFont(fTxt);
  hHistA -> GetYaxis() -> SetLabelSize(fLab[1]);
  hHistA -> GetYaxis() -> CenterTitle(fCnt);
  hHistB -> SetMarkerColor(fCol[1]);
  hHistB -> SetMarkerStyle(fMar[1]);
  hHistB -> SetFillColor(fCol[1]);
  hHistB -> SetFillStyle(fFil[1]);
  hHistB -> SetLineColor(fCol[1]);
  hHistB -> SetLineStyle(fLin[1]);
  hHistB -> SetLineWidth(fWid[1]);
  hHistB -> SetTitle(sTitle.Data());
  hHistB -> SetTitleFont(fTxt);
  hHistB -> GetXaxis() -> SetTitle(sTitleX.Data());
  hHistB -> GetXaxis() -> SetTitleFont(fTxt);
  hHistB -> GetXaxis() -> SetTitleSize(fTit[1]);
  hHistB -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hHistB -> GetXaxis() -> SetLabelFont(fTxt);
  hHistB -> GetXaxis() -> SetLabelSize(fLab[1]);
  hHistB -> GetXaxis() -> CenterTitle(fCnt);
  hHistB -> GetYaxis() -> SetTitle(sTitleY.Data());
  hHistB -> GetYaxis() -> SetTitleFont(fTxt);
  hHistB -> GetYaxis() -> SetTitleSize(fTit[1]);
  hHistB -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hHistB -> GetYaxis() -> SetLabelFont(fTxt);
  hHistB -> GetYaxis() -> SetLabelSize(fLab[1]);
  hHistB -> GetYaxis() -> CenterTitle(fCnt);
  hRatio -> SetMarkerColor(fCol[2]);
  hRatio -> SetMarkerStyle(fMar[2]);
  hRatio -> SetFillColor(fCol[2]);
  hRatio -> SetFillStyle(fFil[2]);
  hRatio -> SetLineColor(fCol[2]);
  hRatio -> SetLineStyle(fLin[2]);
  hRatio -> SetLineWidth(fWid[2]);
  hRatio -> SetTitle(sTitle.Data());
  hRatio -> SetTitleFont(fTxt);
  hRatio -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatio -> GetXaxis() -> SetTitleFont(fTxt);
  hRatio -> GetXaxis() -> SetTitleSize(fTit[0]);
  hRatio -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hRatio -> GetXaxis() -> SetLabelFont(fTxt);
  hRatio -> GetXaxis() -> SetLabelSize(fLab[0]);
  hRatio -> GetXaxis() -> CenterTitle(fCnt);
  hRatio -> GetYaxis() -> SetTitle(sTitleR.Data());
  hRatio -> GetYaxis() -> SetTitleFont(fTxt);
  hRatio -> GetYaxis() -> SetTitleSize(fTit[0]);
  hRatio -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hRatio -> GetYaxis() -> SetLabelFont(fTxt);
  hRatio -> GetYaxis() -> SetLabelSize(fLab[0]);
  hRatio -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;


  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.3};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hHistB, sLabelB.Data());
  leg -> AddEntry(hHistA, sLabelA.Data());
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NPlot * NPlot] = {0.3, 0.1, 0.5, 0.3};
  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(1);
  const UInt_t  fLinLi(2);
  const UInt_t  fWidLi(1);
  const Float_t fLinXY[NPlot * NPlot] = {plotRange[0], 1., plotRange[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;


  // make plot
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.05);
  const Float_t fMarginT1(0.);
  const Float_t fMarginT2(0.05);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.);
  const Float_t fPadXY1[NPlot * NPlot] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NPlot * NPlot] = {0., 0.35, 1., 1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pPad1 = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2 = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot  -> SetGrid(fGrid, fGrid);
  cPlot  -> SetTicks(fTick, fTick);
  cPlot  -> SetBorderMode(fMode);
  cPlot  -> SetBorderSize(fBord);
  pPad1  -> SetGrid(fGrid, fGrid);
  pPad1  -> SetTicks(fTick, fTick);
  pPad1  -> SetLogx(fLogX);
  pPad1  -> SetLogy(fLogY);
  pPad1  -> SetBorderMode(fMode);
  pPad1  -> SetBorderSize(fBord);
  pPad1  -> SetFrameBorderMode(fFrame);
  pPad1  -> SetLeftMargin(fMarginL);
  pPad1  -> SetRightMargin(fMarginR);
  pPad1  -> SetTopMargin(fMarginT1);
  pPad1  -> SetBottomMargin(fMarginB1);
  pPad2  -> SetGrid(fGrid, fGrid);
  pPad2  -> SetTicks(fTick, fTick);
  pPad2  -> SetLogx(fLogX);
  pPad2  -> SetLogy(fLogY);
  pPad2  -> SetBorderMode(fMode);
  pPad2  -> SetBorderSize(fBord);
  pPad2  -> SetFrameBorderMode(fFrame);
  pPad2  -> SetLeftMargin(fMarginL);
  pPad2  -> SetRightMargin(fMarginR);
  pPad2  -> SetTopMargin(fMarginT2);
  pPad2  -> SetBottomMargin(fMarginB2);
  cPlot  -> cd();
  pPad1  -> Draw();
  pPad2  -> Draw();
  pPad1  -> cd();
  hRatio -> Draw("E2");
  line   -> Draw();
  pPad2  -> cd();
  hHistA -> Draw("E2");
  hHistB -> Draw("SAME E2");
  leg    -> Draw();
  txt    -> Draw();
  fOut   -> cd();
  cPlot  -> Write();
  cPlot  -> Close();
  cout << "    Made plot." << endl;


  // close files
  fOut   -> cd();
  hHistA -> Write();
  hHistB -> Write();
  hRatio -> Write();
  fOut   -> Close();
  fInA   -> cd();
  fInA   -> Close();
  fInB   -> cd();
  fInB   -> Close();
  cout << "  Plot made!\n" << endl;

}

// End ------------------------------------------------------------------------
