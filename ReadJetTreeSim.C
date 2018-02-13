// 'ReadJetTreeSim.C'
// Nihar Sahoo, Derek Anderson
// 11.09.2016

#include <vector>
#include <cassert>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TDirectory.h"

using namespace std;


static const Double_t pi = TMath::Pi();

// filepaths
static const TString iPath("pythia.pi0eTtrg9.r05rm1full.d14m11y2017.root");
static const TString oPath("pythia.pi0eTtrg9.r05a065rm1full.plots.d27m11y2017.root");

// parameters
static const Double_t GamTsp = 22;
static const Double_t Pi0Tsp = 111;
static const Double_t MinPt  = 9.;
static const Double_t MaxPt  = 20.;

// jet parameters
static const Double_t Rjet     = 0.5;
static const Double_t MinArea  = 0.65;  // R03: 0.2, R04: 0.5, R05: 0.65, R07: 1.2
static const Double_t MinJetPt = 0.2;
static const Double_t RecoilDf = pi / 4.;



void ReadJetTreeSim(Bool_t inBatchMode=false) {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\nBeginning reading script!" << endl;


  TFile *fIn  = new TFile(iPath, "read");
  TFile *fOut = new TFile(oPath, "recreate");
  if (!fIn) {
    cerr << "PANIC: input file could not be opened!" << endl;
    assert(fIn);
  }

  // get trees
  TTree *JetTree = 0;
  if (fIn)
    fIn -> GetObject("femtoDST", JetTree);
  if (!JetTree) {
    cerr << "PANIC: 'femtoDST' not grabbed!" << endl;
    assert(JetTree);
  }


  // declare Geant event leaves
  cout << "  Setting branch addresses..." << endl;
  Int_t    EventIndex = 0;
  Double_t Refmult    = 0.;
  Double_t TSP        = 0.;
  Double_t TrgEta     = 0.;
  Double_t TrgPhi     = 0.;
  Double_t TrgEt      = 0.;
  Double_t Rho        = 0.;
  Double_t Sigma      = 0.;
  Double_t Vz         = 0.;
  // declare Geant jet leaves
  vector<Double_t> *JetEta    = 0;
  vector<Double_t> *JetPt     = 0;
  vector<Double_t> *JetNCons  = 0;
  vector<Double_t> *JetIndex  = 0;
  vector<Double_t> *JetEta    = 0;
  vector<Double_t> *JetPhi    = 0;
  vector<Double_t> *JetE      = 0;
  vector<Double_t> *JetArea   = 0;
  // declare Geant constituent leaves
  vector<vector<Double_t>> *JetConsPt  = 0;
  vector<vector<Double_t>> *JetConsEta = 0;
  vector<vector<Double_t>> *JetConsPhi = 0;
  vector<vector<Double_t>> *JetConsE   = 0;
  


  // declare Geant branches
  TBranch *bEventIndex = 0;
  TBranch *bRefmult    = 0;
  TBranch *bTsp        = 0;
  TBranch *bTrgEta     = 0;
  TBranch *bTrgPhi     = 0;
  TBranch *bTrgEt      = 0;
  TBranch *bRho        = 0;
  TBranch *bSigma      = 0;
  TBranch *bVz         = 0;
  TBranch *bJetEta     = 0;
  TBranch *bJetPt      = 0;
  TBranch *bJetNCons   = 0;
  TBranch *bJetIndex   = 0;
  TBranch *bJetEta     = 0;
  TBranch *bJetPhi     = 0;
  TBranch *bJetE       = 0;
  TBranch *bJetArea    = 0;
  TBranch *bJetConsPt  = 0;
  TBranch *bJetConsEta = 0;
  TBranch *bJetConsPhi = 0;
  TBranch *bJetConsE   = 0;


  // set Geant branches
  JetTree -> SetBranchAddress("eventIndex", &EventIndex, &bEventIndex);
  JetTree -> SetBranchAddress("Refmult", &Refmult, &bRefmult);
  JetTree -> SetBranchAddress("TSP", &TSP, &bTsp);
  JetTree -> SetBranchAddress("TrgEta", &TrgEta, &bTrgEta);
  JetTree -> SetBranchAddress("TrgPhi", &TrgPhi, &bTrgPhi);
  JetTree -> SetBranchAddress("TrgEt", &TrgEt, &bTrgEt);
  JetTree -> SetBranchAddress("Rho", &Rho, &bRho);
  JetTree -> SetBranchAddress("Sigma", &Sigma, &bSigma);
  JetTree -> SetBranchAddress("Vz", &Vz,&bVz);
  JetTree -> SetBranchAddress("JetIndex", &JetIndex, &bJetIndex);
  JetTree -> SetBranchAddress("JetEta", &JetEta, &bJetEta);
  JetTree -> SetBranchAddress("JetPt", &JetPt, &bJetPt);
  JetTree -> SetBranchAddress("JetNCons", &JetNCons, &bJetNCons);
  JetTree -> SetBranchAddress("JetEta", &JetEta, &bJetEta);
  JetTree -> SetBranchAddress("JetPhi",&JetPhi, &bJetPhi); 
  JetTree -> SetBranchAddress("JetE", &JetE, &bJetE); 
  JetTree -> SetBranchAddress("JetArea",&JetArea, &bJetArea);
  JetTree -> SetBranchAddress("JetConsPt", &JetConsPt, &bJetConsPt);
  JetTree -> SetBranchAddress("JetConsEta", &JetConsEta, &bJetConsEta);
  JetTree -> SetBranchAddress("JetConsPhi", &JetConsPhi, &bJetConsPhi);
  JetTree -> SetBranchAddress("JetConsE", &JetConsE, &bJetConsE);



  // declare histograms [0 for pi0 trigger, 1 for gamma trigger]
  cout << "  Creating histograms..." << endl;
  TH1D *hRefmult[2];
  TH1D *hNumJets[2];
  TH1D *hRho[2];
  TH1D *hSigma[2];
  TH1D *hJetArea[2];
  TH1D *hJetEta[2];
  TH1D *hAllPtRaw[2];
  TH1D *hJetPtRaw[2];
  TH1D *hAllDeltaPhi[2];
  TH1D *hBinDeltaPhi[2][3];
  TH1D *hJetDeltaPhi[2];
  TH1D *hJetPtCorr[2];
  TH2D *hJetPtVsDf[2];
  TH2D *hJetPtVsA[2];

  const Int_t    nM  = 200;
  const Int_t    nN  = 50;
  const Int_t    nP  = 50;
  const Int_t    nA  = 500;
  const Int_t    nH  = 100;
  const Int_t    nDf = 360;
  const Int_t    nPt = 120;
  const Double_t m1  = 0.;
  const Double_t m2  = 200.;
  const Double_t n1  = 0.;
  const Double_t n2  = 50.;
  const Double_t p1  = 0.;
  const Double_t p2  = 5.;
  const Double_t a1  = 0.;
  const Double_t a2  = 5.;
  const Double_t h1  = -5.;
  const Double_t h2  = 5.;
  const Double_t dF1 = -1.* (pi / 2.);
  const Double_t dF2 = 3.* (pi / 2.);
  const Double_t pT1 = -20.;
  const Double_t pT2 = 100.;
  hRefmult[0]        = new TH1D("hRefmultP", "Refmult, #pi^{0} trigger; refmult; counts", nM, m1, m2);
  hRefmult[1]        = new TH1D("hRefmultG", "Refmult, #gamma^{rich} trigger; refmult; counts", nM, m1, m2);
  hNumJets[0]        = new TH1D("hNumJetsP", "No. of jets, #pi^{0} trigger; N_{jet}; counts", nN, n1, n2);
  hNumJets[1]        = new TH1D("hNumJetsG", "No. of jets, #gamma^{rich} trigger; N_{jet}; counts", nN, n1, n2);
  hRho[0]            = new TH1D("hRhoP", "Rho, #pi^{0} trigger; #rho; counts", nP, p1, p2);
  hRho[1]            = new TH1D("hRhoG", "Rho, #gamma^{rich} trigger; #rho; counts", nP, p1, p2);
  hSigma[0]          = new TH1D("hSigmaP", "Sigma, #pi^{0} trigger; #sigma; counts", nP, p1, p2);
  hSigma[1]          = new TH1D("hSigmaG", "Sigma, #gamma^{rich} trigger; #sigma; counts)", nP, p1, p2);
  hJetArea[0]        = new TH1D("hJetAreaP", "Recoil jet area, #pi^{0} trigger; A_{jet}; (1/N_{trg})dN_{jet}/dA_{jet}", nA, a1, a2);
  hJetArea[1]        = new TH1D("hJetAreaG", "Recoil jet area, #gamma^{rich} trigger; A_{jet}; (1/N_{trg})dN_{jet}/dA_{jet}", nA, a1, a2);
  hJetEta[0]         = new TH1D("hJetEtaP", "Recoil jet #eta, #pi^{0} trigger; #eta; (1/N_{trg})dN_{jet}/d#eta", nH, h1, h2);
  hJetEta[1]         = new TH1D("hJetEtaG", "Recoil jet #eta, #gamma^{rich} trigger; #eta; (1/N_{trg})dN_{jet}/d#eta", nH, h1, h2);
  hAllDeltaPhi[0]    = new TH1D("hAllDeltaPhiP", "All jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hAllDeltaPhi[1]    = new TH1D("hAllDeltaPhiG", "All jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[0][0] = new TH1D("hBinDeltaPhiP_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[0][1] = new TH1D("hBinDeltaPhiP_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[0][2] = new TH1D("hBinDeltaPhiP_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[1][0] = new TH1D("hBinDeltaPhiG_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[1][1] = new TH1D("hBinDeltaPhiG_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[1][2] = new TH1D("hBinDeltaPhiG_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hJetDeltaPhi[0]    = new TH1D("hJetDeltaPhiP", "Recoil jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hJetDeltaPhi[1]    = new TH1D("hJetDeltaPhiG", "Recoil jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hAllPtRaw[0]       = new TH1D("hAllPtRawP", "All jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
  hAllPtRaw[1]       = new TH1D("hAllPtRawG", "All jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
  hJetPtRaw[0]       = new TH1D("hJetPtRawP", "Recoil jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
  hJetPtRaw[1]       = new TH1D("hJetPtRawG", "Recoil jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
  hJetPtCorr[0]      = new TH1D("hJetPtCorrP", "Recoil jet p_{T}^{corr}, #pi^{0} trigger; p_{T}^{corr}; (1/N_{trg})dN_{jet}/dp_{T}^{corr}", nPt, pT1, pT2);
  hJetPtCorr[1]      = new TH1D("hJetPtCorrG", "Recoil jet p_{T}^{corr}, #gamma^{rich} trigger; p_{T}^{corr}; (1/N_{trg})dN_{jet}/dp_{T}^{corr}", nPt, pT1, pT2);
  hJetPtVsDf[0]      = new TH2D("hJetPtVsDfP", "Jet p_{T}^{raw} vs. #Delta#varphi, #pi^{0} trigger; #Delta#varphi; p_{T}^{raw}", nDf, dF1, dF2, nPt, pT1, pT2);
  hJetPtVsDf[1]      = new TH2D("hJetPtVsDfG", "Jet p_{T}^{raw} vs. #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; p_{T}^{raw}", nDf, dF1, dF2, nPt, pT1, pT2);
  hJetPtVsA[0]       = new TH2D("hJetPtVsAp", "Jet p_{T}^{raw} vs. Area, #pi^{0} trigger; A_{jet}; p_{T}^{raw}", nA, a1, a2, nPt, pT1, pT2);
  hJetPtVsA[1]       = new TH2D("hJetPtVsAg", "Jet p_{T}^{raw} vs. Area, #gamma^{rich} trigger; A_{jet}; p_{T}^{raw}", nA, a1, a2, nPt, pT1, pT2);
  // errors
  for (Int_t i = 0; i < 2; i++) {
    hRefmult[i]        -> Sumw2();
    hNumJets[i]        -> Sumw2();
    hRho[i]            -> Sumw2();
    hSigma[i]          -> Sumw2();
    hJetArea[i]        -> Sumw2();
    hJetEta[i]         -> Sumw2();
    hAllDeltaPhi[i]    -> Sumw2();
    hBinDeltaPhi[i][0] -> Sumw2();
    hBinDeltaPhi[i][1] -> Sumw2();
    hBinDeltaPhi[i][2] -> Sumw2();
    hJetDeltaPhi[i]    -> Sumw2();
    hAllPtRaw[i]       -> Sumw2();
    hJetPtRaw[i]       -> Sumw2();
    hJetPtCorr[i]      -> Sumw2();
    hJetPtVsDf[i]      -> Sumw2();
    hJetPtVsA[i]       -> Sumw2();
  }



  cout << "  Beginning event loop..." << endl;
  Int_t nEvts = (Int_t) JetTree -> GetEntries();

  // event loop
  Int_t nTrgPi0  = 0;
  Int_t nTrgGam  = 0;
  Int_t nBytes   = 0;
  Int_t breakVal = 0;
  for (Int_t i = 0; i < nEvts; i++) {

    Int_t bytes = JetTree -> GetEntry(i);
    if (bytes < 0) {
      cerr << "ERROR: problem with event " << i << "...\n" << endl;
      breakVal = 1;
      break;
    }

    nBytes += bytes;
    if (!inBatchMode) {
      cout << "    Processing event " << i+1 << "/" << nEvts << "...\r" << flush;
      if (i+1 == nEvts) cout << endl;
    }
    else
      cout << "    Processing event " << i+1 << "/" << nEvts << "..." << endl;


    // determine triggers [0 for pi0, 1 for gamma-rich]
    Int_t  tspFlag = -1;
    if (TSP == 111)
      tspFlag = 0;
    else if (TSP == 22)
      tspFlag = 1;


    // fill event histograms
    if (tspFlag == 0) {
      hRefmult[0] -> Fill(Refmult);
      hRho[0]     -> Fill(Rho);
      hSigma[0]   -> Fill(Sigma);
      ++nTrgPi0;
    }
    else if (tspFlag == 1) {
      hRefmult[1] -> Fill(Refmult);
      hRho[1]     -> Fill(Rho);
      hSigma[1]   -> Fill(Sigma);
      ++nTrgGam;
    }


    // jet loop
    Int_t nJets = (Int_t) JetEta -> size();
    for (Int_t j = 0; j < nJets; j++) {

      Double_t h     = JetEta  -> at(j);
      Double_t a     = JetArea -> at(j);
      Double_t pT    = JetPt   -> at(j);
      Double_t f     = JetPhi  -> at(j);
      Double_t pTc   = pT - (Rho * a);
      Double_t dF    = f - TrgPhi;
      if (dF < -1.*pi/2.) dF += 2.*pi;
      if (dF > 3.*pi/2.)  dF -= 2.*pi;
      Double_t dFre  = f - TrgPhi;
      if (dFre < 0.)    dFre += 2.*pi;
      if (dFre > 2.*pi) dFre -= 2.*pi;
      Double_t dFcut = abs(dFre - pi);

      // fill histograms for all jets
      Bool_t isInBin1 = ((pT > 0.2) && (pT < 1.));
      Bool_t isInBin2 = ((pT > 1.) && (pT < 5.));
      Bool_t isInBin3 = (pT > 5.);
      Bool_t isRecoil = (dFcut < RecoilDf);
      if (tspFlag == 0) {
        hAllPtRaw[0]    -> Fill(pT);
        hAllDeltaPhi[0] -> Fill(dF);
        hJetPtVsDf[0]   -> Fill(dF, pT);
        hJetPtVsA[0]    -> Fill(a, pT);
        if (isInBin1)
          hBinDeltaPhi[0][0] -> Fill(dF);
        if (isInBin2)
          hBinDeltaPhi[0][1] -> Fill(dF);
        if (isInBin3)
          hBinDeltaPhi[0][2] -> Fill(dF);
      }
      if (tspFlag == 1) {
        hAllPtRaw[1]    -> Fill(pT);
        hAllDeltaPhi[1] -> Fill(dF);
        hJetPtVsDf[1]   -> Fill(dF, pT);
        hJetPtVsA[1]    -> Fill(a, pT);
        if (isInBin1)
          hBinDeltaPhi[1][0] -> Fill(dF);
        if (isInBin2)
          hBinDeltaPhi[1][1] -> Fill(dF);
        if (isInBin3)
          hBinDeltaPhi[1][2] -> Fill(dF);
      }


      // jet cuts
      if (pT < MinJetPt)
        continue;
      if (a < MinArea)
        continue;
      if (!isRecoil)
        continue;

      // fill histograms
      if (tspFlag == 0) {
        hJetArea[0]     -> Fill(a);
        hJetEta[0]      -> Fill(h);
        hJetDeltaPhi[0] -> Fill(dF);
        hJetPtRaw[0]    -> Fill(pT);
        hJetPtCorr[0]   -> Fill(pTc);
      }
      else if (tspFlag == 1) {
        hJetArea[1]     -> Fill(a);
        hJetEta[1]      -> Fill(h);
        hJetDeltaPhi[1] -> Fill(dF);
        hJetPtRaw[1]    -> Fill(pT);
        hJetPtCorr[1]   -> Fill(pTc);
      }

    }  // end jet loop

    if (tspFlag == 0)
      hNumJets[0] -> Fill(nJets);
    else if (tspFlag == 1)
      hNumJets[1] -> Fill(nJets);

  }  // end event loop

  if (breakVal == 1) {
    cerr << "ERROR: Occured during event loop!\n"
         << "       Aborting program!"
         << endl;
    cassert(0);
  }
  else
    cout << "  Event loop finished!" << endl;



  // normalize histograms
  cout << "  Normalizing histograms..." << endl;
  const Double_t mBin  = (m2 - m1) / nM;
  const Double_t nBin  = (n2 - n1) / nN;
  const Double_t pBin  = (p2 - p1) / nP;
  const Double_t aBin  = (a2 - a1) / nA;
  const Double_t hBin  = (h2 - h1) / nH;
  const Double_t dFbin = (dF2 - dF1) / nDf;
  const Double_t pTbin = (pT2 - pT1) / nPt;
  const Double_t hNorm = 2. * (1. - Rjet);
  for (Int_t i = 0; i < 2; i++) {
    Double_t norm = 1.;
/*
    if (i == 0)
      norm = nTrgPi0;
    else if (i == 1)
      norm = nTrgGam;
*/
    hRefmult[i]        -> Scale(1. / mBin);
    hNumJets[i]        -> Scale(1. / nBin);
    hRho[i]            -> Scale(1. / pBin);
    hSigma[i]          -> Scale(1. / pBin);
    hJetArea[i]        -> Scale(1. / norm);
    hJetArea[i]        -> Scale(1. / hNorm);
    hJetArea[i]        -> Scale(1. / aBin);
    hJetEta[i]         -> Scale(1. / norm);
    hJetEta[i]         -> Scale(1. / hNorm);
    hJetEta[i]         -> Scale(1. / hBin);
    hAllDeltaPhi[i]    -> Scale(1. / norm);
    hAllDeltaPhi[i]    -> Scale(1. / hNorm);
    hAllDeltaPhi[i]    -> Scale(1. / dFbin);
    hBinDeltaPhi[i][0] -> Scale(1. / norm);
    hBinDeltaPhi[i][0] -> Scale(1. / hNorm);
    hBinDeltaPhi[i][0] -> Scale(1. / dFbin);
    hBinDeltaPhi[i][1] -> Scale(1. / norm);
    hBinDeltaPhi[i][1] -> Scale(1. / hNorm);
    hBinDeltaPhi[i][1] -> Scale(1. / dFbin);
    hBinDeltaPhi[i][2] -> Scale(1. / norm);
    hBinDeltaPhi[i][2] -> Scale(1. / hNorm);
    hBinDeltaPhi[i][2] -> Scale(1. / dFbin);
    hJetDeltaPhi[i]    -> Scale(1. / norm);
    hJetDeltaPhi[i]    -> Scale(1. / hNorm);
    hJetDeltaPhi[i]    -> Scale(1. / dFbin);
    hAllPtRaw[i]       -> Scale(1. / norm);
    hAllPtRaw[i]       -> Scale(1. / hNorm);
    hAllPtRaw[i]       -> Scale(1. / pTbin);
    hJetPtRaw[i]       -> Scale(1. / norm);
    hJetPtRaw[i]       -> Scale(1. / hNorm);
    hJetPtRaw[i]       -> Scale(1. / pTbin);
    hJetPtCorr[i]      -> Scale(1. / norm);
    hJetPtCorr[i]      -> Scale(1. / hNorm);
    hJetPtCorr[i]      -> Scale(1. / pTbin);
    hJetPtVsDf[i]      -> Scale(1. / norm);
    hJetPtVsDf[i]      -> Scale(1. / hNorm);
    hJetPtVsDf[i]      -> Scale(1. / dFbin);
    hJetPtVsDf[i]      -> Scale(1. / pTbin);
    hJetPtVsA[i]       -> Scale(1. / norm);
    hJetPtVsA[i]       -> Scale(1. / hNorm);
    hJetPtVsA[i]       -> Scale(1. / aBin);
    hJetPtVsA[i]       -> Scale(1. / pTbin);
  }


  // create directories
  TDirectory *dPi0 = (TDirectory*) fOut -> mkdir("Pi0");
  TDirectory *dGam = (TDirectory*) fOut -> mkdir("Gam");

  // write and close output
  fOut -> cd();
  for (Int_t i = 0; i < 2; i++) {
    if (i == 0)
      dPi0 -> cd();
    else if (i == 1)
      dGam -> cd();
    hRefmult[i]        -> Write();
    hNumJets[i]        -> Write();
    hRho[i]            -> Write();
    hSigma[i]          -> Write();
    hJetArea[i]        -> Write();
    hJetEta[i]         -> Write();
    hAllDeltaPhi[i]    -> Write();
    hBinDeltaPhi[i][0] -> Write();
    hBinDeltaPhi[i][1] -> Write();
    hBinDeltaPhi[i][2] -> Write();
    hJetDeltaPhi[i]    -> Write();
    hAllPtRaw[i]       -> Write();
    hJetPtRaw[i]       -> Write();
    hJetPtCorr[i]      -> Write();
    hJetPtVsDf[i]      -> Write();
    hJetPtVsA[i]       -> Write();
    fOut -> cd();
  }
  fOut -> Close();

  // close input
  fIn -> cd();
  fIn -> Close();


  cout << "Reading script finished!\n" << endl;

}

// End ------------------------------------------------------------------------