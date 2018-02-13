// 'MatchJets.C'
// Nihar Sahoo, Derek Anderson
// 11.15.2016
//
// This version was made to compare charged jets to neutral jets 
//   [Derek, 01.28.2017]

#include <vector>
#include <cassert>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TDirectory.h"

using namespace std;


static const Double_t pi = TMath::Pi();

// filepaths
static const TString cName("pythiaCheckG.r03a02rm1chrg.d27m1y2017.root");
static const TString fName("pythiaCheckG.r03a02rm1full.d27m1y2017.root");
static const TString oName("debug.root");

// parameters
static const Double_t MinJetPt = 0.2;
static const Double_t MinArea  = 0.2;  // R03: 0.2, R04: 0.5, R05: 0.65, R07: 1.2
static const Double_t Rcut    = 0.3;   // Rcut = Rjet
static const Double_t Qmin    = 0.15;  // fraction of jet pT must be above Qmin
static const Double_t Qmax    = 1.5;   // fraction of jet pT must be below Qmax
static const Double_t HardCut = 17.;   // jets w/ pT>HardCut are considered 'hard'
static const Double_t Ccut    = 0.;    // pTchrg must be above this
static const Double_t Fcut    = 0.;    // pTfull must be above this



void MatchJets(const TString cPath=cName, const TString fPath=fName, const TString oPath=oName, Bool_t inBatchMode=false) {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;

  cout << "\nBeginning match script!" << endl;


  TFile *fChrg = new TFile(cPath, "read");
  TFile *fFull = new TFile(fPath, "read");
  TFile *fOut  = new TFile(oPath, "recreate");
  if (!fChrg) {
    cerr << "PANIC: Charged input file could not be opened!" << endl;
    assert(fChrg);
  }
  if (!fFull) {
    cerr << "PANIC: Full input file could not be opened!" << endl;
    assert(fFull);
  }

  // get trees
  TTree *tChrg;
  TTree *tFull;
  if (fChrg && fFull) {
    fChrg -> GetObject("femtoDST", tChrg);
    fFull -> GetObject("femtoDST", tFull);
  }


  // declare charged event leaves
  cout << "  Setting branch addresses..." << endl;
  Int_t    cEventIndex = 0;
  Int_t    cNJets      = 0;
  Double_t cRefmult    = 0.;
  Double_t cTSP        = 0.;
  Double_t cTrgEta     = 0.;
  Double_t cTrgPhi     = 0.;
  Double_t cTrgEt      = 0.;
  Double_t cRho        = 0.;
  Double_t cSigma      = 0.;
  // declare charged jet leaves
  vector<Double_t> *cJetEta    = 0;
  vector<Double_t> *cJetPt     = 0;
  vector<Double_t> *cJetNCons  = 0;
  vector<Double_t> *cJetIndex  = 0;
  vector<Double_t> *cJetEta    = 0;
  vector<Double_t> *cJetPhi    = 0;
  vector<Double_t> *cJetE      = 0;
  vector<Double_t> *cJetArea   = 0;
  // declare charged constituent leaves
  vector<vector<Double_t>> *cJetConsPt  = 0;
  vector<vector<Double_t>> *cJetConsEta = 0;
  vector<vector<Double_t>> *cJetConsPhi = 0;
  vector<vector<Double_t>> *cJetConsE   = 0;
  
  // declare full event leaves
  Int_t    fEventIndex = 0;
  Int_t    fNJets      = 0;
  Double_t fRefmult    = 0.;
  Double_t fTSP        = 0.;
  Double_t fTrgEta     = 0.;
  Double_t fTrgPhi     = 0.;
  Double_t fTrgEt      = 0.;
  Double_t fRho        = 0.;
  Double_t fSigma      = 0.;
  // declare full jet leaves
  vector<Double_t> *fJetEta    = 0;
  vector<Double_t> *fJetPt     = 0;
  vector<Double_t> *fJetNCons  = 0;
  vector<Double_t> *fJetIndex  = 0;
  vector<Double_t> *fJetEta    = 0;
  vector<Double_t> *fJetPhi    = 0;
  vector<Double_t> *fJetE      = 0;
  vector<Double_t> *fJetArea   = 0;
  // declare full constituent leaves
  vector<vector<Double_t>> *fJetConsPt  = 0;
  vector<vector<Double_t>> *fJetConsEta = 0;
  vector<vector<Double_t>> *fJetConsPhi = 0;
  vector<vector<Double_t>> *fJetConsE   = 0;


  // declare charged branches
  TBranch *bEventIndexC = 0;
  TBranch *bNJetsC      = 0;
  TBranch *bRefmultC    = 0;
  TBranch *bTspC        = 0;
  TBranch *bTrgEtaC     = 0;
  TBranch *bTrgPhiC     = 0;
  TBranch *bTrgEtC      = 0;
  TBranch *bRhoC        = 0;
  TBranch *bSigmaC      = 0;
  TBranch *bJetEtaC     = 0;
  TBranch *bJetPtC      = 0;
  TBranch *bJetNConsC   = 0;
  TBranch *bJetIndexC   = 0;
  TBranch *bJetEtaC     = 0;
  TBranch *bJetPhiC     = 0;
  TBranch *bJetEC       = 0;
  TBranch *bJetAreaC    = 0;
  TBranch *bJetConsPtC  = 0;
  TBranch *bJetConsEtaC = 0;
  TBranch *bJetConsPhiC = 0;
  TBranch *bJetConsEC   = 0;

  // declare full branches
  TBranch *bEventIndexF = 0;
  TBranch *bNJetsF      = 0;
  TBranch *bRefmultF    = 0;
  TBranch *bTspF        = 0;
  TBranch *bTrgEtaF     = 0;
  TBranch *bTrgPhiF     = 0;
  TBranch *bTrgEtF      = 0;
  TBranch *bRhoF        = 0;
  TBranch *bSigmaF      = 0;
  TBranch *bJetEtaF     = 0;
  TBranch *bJetPtF      = 0;
  TBranch *bJetNConsF   = 0;
  TBranch *bJetIndexF   = 0;
  TBranch *bJetEtaF     = 0;
  TBranch *bJetPhiF     = 0;
  TBranch *bJetEF       = 0;
  TBranch *bJetAreaF    = 0;
  TBranch *bJetConsPtF  = 0;
  TBranch *bJetConsEtaF = 0;
  TBranch *bJetConsPhiF = 0;
  TBranch *bJetConsEF   = 0;


  // set charged branches
  tChrg -> SetBranchAddress("EventIndex", &cEventIndex, &bEventIndexC);
  tChrg -> SetBranchAddress("Refmult", &cRefmult, &bRefmultC);
  tChrg -> SetBranchAddress("NJets", &cNJets, &bNJetsC);
  tChrg -> SetBranchAddress("TSP", &cTSP, &bTspC);
  tChrg -> SetBranchAddress("TrgEta", &cTrgEta, &bTrgEtaC);
  tChrg -> SetBranchAddress("TrgPhi", &cTrgPhi, &bTrgPhiC);
  tChrg -> SetBranchAddress("TrgEt", &cTrgEt, &bTrgEtC);
  tChrg -> SetBranchAddress("Rho", &cRho, &bRhoC);
  tChrg -> SetBranchAddress("Sigma", &cSigma, &bSigmaC);
  tChrg -> SetBranchAddress("JetIndex", &cJetIndex, &bJetIndexC);
  tChrg -> SetBranchAddress("JetEta", &cJetEta, &bJetEtaC);
  tChrg -> SetBranchAddress("JetPt", &cJetPt, &bJetPtC);
  tChrg -> SetBranchAddress("JetNCons", &cJetNCons, &bJetNConsC);
  tChrg -> SetBranchAddress("JetEta", &cJetEta, &bJetEtaC);
  tChrg -> SetBranchAddress("JetPhi",&cJetPhi, &bJetPhiC); 
  tChrg -> SetBranchAddress("JetE", &cJetE, &bJetEC); 
  tChrg -> SetBranchAddress("JetArea",&cJetArea, &bJetAreaC);
  tChrg -> SetBranchAddress("JetConsPt", &cJetConsPt, &bJetConsPtC);
  tChrg -> SetBranchAddress("JetConsEta", &cJetConsEta, &bJetConsEtaC);
  tChrg -> SetBranchAddress("JetConsPhi", &cJetConsPhi, &bJetConsPhiC);
  tChrg -> SetBranchAddress("JetConsE", &cJetConsE, &bJetConsEC);

  // set full branches
  tFull -> SetBranchAddress("EventIndex", &fEventIndex, &bEventIndexF);
  tFull -> SetBranchAddress("Refmult", &fRefmult, &bRefmultF);
  tFull -> SetBranchAddress("NJets", &fNJets, &bNJetsF);
  tFull -> SetBranchAddress("TSP", &fTSP, &bTspF);
  tFull -> SetBranchAddress("TrgEta", &fTrgEta, &bTrgEtaF);
  tFull -> SetBranchAddress("TrgPhi", &fTrgPhi, &bTrgPhiF);
  tFull -> SetBranchAddress("TrgEt", &fTrgEt, &bTrgEtF);
  tFull -> SetBranchAddress("Rho", &fRho, &bRhoF);
  tFull -> SetBranchAddress("Sigma", &fSigma, &bSigmaF);
  tFull -> SetBranchAddress("JetIndex", &fJetIndex, &bJetIndexF);
  tFull -> SetBranchAddress("JetEta", &fJetEta, &bJetEtaF);
  tFull -> SetBranchAddress("JetPt", &fJetPt, &bJetPtF);
  tFull -> SetBranchAddress("JetNCons", &fJetNCons, &bJetNConsF);
  tFull -> SetBranchAddress("JetEta", &fJetEta, &bJetEtaF);
  tFull -> SetBranchAddress("JetPhi",&fJetPhi, &bJetPhiF); 
  tFull -> SetBranchAddress("JetE", &fJetE, &bJetEF); 
  tFull -> SetBranchAddress("JetArea",&fJetArea, &bJetAreaF);
  tFull -> SetBranchAddress("JetConsPt", &fJetConsPt, &bJetConsPtF);
  tFull -> SetBranchAddress("JetConsEta", &fJetConsEta, &bJetConsEtaF);
  tFull -> SetBranchAddress("JetConsPhi", &fJetConsPhi, &bJetConsPhiF);
  tFull -> SetBranchAddress("JetConsE", &fJetConsE, &bJetConsEF);



  cout << "  Creating histograms..." << endl;

  const Int_t nJetTypes   = 7;
  const Int_t nMatchTypes = 5;
  TH1D *hEfficiencyA;
  TH1D *hEfficiencyPt;
  TH2D *hResponseA;
  TH2D *hResponseAn;
  TH2D *hResponsePt;
  TH2D *hResponsePtn;
  TH2D *hResponsePtc;
  TH2D *hResponsePtcn;
  // event histograms
  TH1D *hRefmultC;
  TH1D *hRefmultF;
  TH1D *hNumJetsC;
  TH1D *hNumJetsF;
  TH1D *hNumToMatch;
  TH1D *hNumMatched;
  TH1D *hNumHardC;
  TH1D *hNumHardF;
  TH1D *hChrgArea;
  TH1D *hMatchArea;
  TH1D *hChrgPtCorr;
  TH1D *hMatchPtCorr;
  // jet histograms  [0='C',1='F',2='A',3='M',4='J',5='Y',6='N']
  TH1D *hJetArea[nJetTypes];
  TH1D *hJetEta[nJetTypes];
  TH1D *hJetPhi[nJetTypes];
  TH1D *hJetPt[nJetTypes];
  TH1D *hJetPtCorr[nJetTypes];
  TH2D *hJetPhiVsEta[nJetTypes];
  // matching histograms [0='F',1='A',2='M',3='J',4='Y']
  TH1D *hJetQt[nMatchTypes];
  TH1D *hJetDr[nMatchTypes];
  TH1D *hJetS[nMatchTypes];
  TH1D *hJetDp[nMatchTypes];
  TH2D *hJetQtVsDr[nMatchTypes];
  TH2D *hJetSvsDr[nMatchTypes];

  // histogram parameters
  const Int_t    nM = 200;
  const Int_t    nN = 100;
  const Int_t    nA = 500;
  const Int_t    nH = 1000;
  const Int_t    nF = 720;
  const Int_t    nP = 1100;
  const Int_t    nQ = 600;
  const Int_t    nR = 600;
  const Int_t    nS = 600;
  const Int_t    nD = 1000;
  const Double_t m1 = 0.;
  const Double_t m2 = 200.;
  const Double_t n1 = 0.;
  const Double_t n2 = 100;
  const Double_t a1 = 0.;
  const Double_t a2 = 5.;
  const Double_t h1 = -5.;
  const Double_t h2 = 5.;
  const Double_t f1 = -2.*pi;
  const Double_t f2 = 2.*pi;
  const Double_t p1 = -10.;
  const Double_t p2 = 100.;
  const Double_t q1 = 0.;
  const Double_t q2 = 3.;
  const Double_t r1 = 0.;
  const Double_t r2 = 3.;
  const Double_t s1 = 0.;
  const Double_t s2 = 3.;
  const Double_t d1 = -50.;
  const Double_t d2 = 50.;
  hEfficiencyA    = new TH1D("hEfficiencyA", "Efficiency, #epsilon(A_{jet}) = N_{match}(A_{jet})/N_{charge}(A_{jet})", nA, a1, a2);
  hEfficiencyPt   = new TH1D("hEfficiencyPt", "Efficiency, #epsilon(p_{T}^{jet}) = N_{match}(p_{T}^{jet})/N_{charge}(p_{T}^{jet})", nP, p1, p2);
  hResponseA      = new TH2D("hResponseA", "Response matrix, jet area; match; charge", nA, a1, a2, nA, a1, a2);
  hResponseAn     = new TH2D("hResponseAn", "Response matrix, jet area (normalized); match; charge", nA, a1, a2, nA, a1, a2);
  hResponsePt     = new TH2D("hResponsePt", "Response matrix, jet p_{T}; match; charge", nP, p1, p2, nP, p1, p2);
  hResponsePtN    = new TH2D("hResponsePtN", "Response matrix, jet p_{T} (normalized); match; charge", nP, p1, p2, nP, p1, p2);
  hResponsePtc    = new TH2D("hResponsePtc", "Response matrix, jet p_{T}^{corr}; match; charge", nP, p1, p2, nP, p1, p2);
  hResponsePtcN   = new TH2D("hResponsePtcN", "Response matrix, jet p_{T}^{corr} (normalized); match; charge", nP, p1, p2, nP, p1, p2);
  // event histograms
  hRefmultC       = new TH1D("hRefmultC", "Charged Refmult", nM, m1, m2);
  hRefmultF       = new TH1D("hRefmultF", "Full Refmult", nM, m1, m2);
  hNumJetsC       = new TH1D("hNumJetsC", "no. of jets, Charged", nN, n1, n2);
  hNumJetsF       = new TH1D("hNumJetsF", "No. of jets, full", nN, n1, n2);
  hNumToMatch     = new TH1D("hNumToMatch", "No. of jets to match (ie. no. of jets w/ pT > pTcut)", nN, n1, n2);
  hNumMatched     = new TH1D("hNumMatched", "No. of jets matched", nN, n1, n2);
  hNumHardC       = new TH1D("hNumHardC", "No. of jets w/ p_{T} above a threshold, charged", nN, n1, n2);
  hNumHardF       = new TH1D("hNumHardF", "No. of jets w/ p_{T} above a threshold, full", nN, n1, n2);
  hChrgArea       = new TH1D("hChrgArea", "Total no. of jets to match per A_{jet} bin (for efficiency)", nA, a1, a2);
  hMatchArea      = new TH1D("hMatchArea", "Total no. of jets matched per A_{jet} bin (for efficiency)", nA, a1, a2);
  hChrgPtCorr     = new TH1D("hChrgPtCorr", "Total no. of jets to match per p_{T}^{corr} bin (for efficiency)", nP, p1, p2);
  hMatchPtCorr    = new TH1D("hMatchPtCorr", "Total no. of jets matched per p_{T}^{corr} bin (for efficiency)", nP, p1, p2);
  // charge jets
  hJetArea[0]     = new TH1D("hJetAreaC", "Jet area, charged", nA, a1, a2);
  hJetEta[0]      = new TH1D("hJetEtaC", "Jet eta, charged", nH, h1, h2);
  hJetPhi[0]      = new TH1D("hJetPhiC", "Jet phi, charged", nF, f1, f2);
  hJetPt[0]       = new TH1D("hJetPtC", "Jet p_{T}, charged", nP, p1, p2);
  hJetPtCorr[0]   = new TH1D("hJetPtCorrC", "Jet p_{T}^{corr}, charged", nP, p1, p2);
  hJetPhiVsEta[0] = new TH2D("hJetPhiVsEtaC", "Jet #varphi vs. #eta, charged", nH, h1, h2, nF, f1, f2);
  // full jets
  hJetArea[1]     = new TH1D("hJetAreaF", "Jet area, full", nA, a1, a2);
  hJetEta[1]      = new TH1D("hJetEtaF", "Jet eta, full", nH, h1, h2);
  hJetPhi[1]      = new TH1D("hJetPhiF", "Jet phi, full", nF, f1, f2);
  hJetPt[1]       = new TH1D("hJetPtF", "Jet p_{T}, full", nP, p1, p2);
  hJetPtCorr[1]   = new TH1D("hJetPtCorrF", "Jet p_{T}^{corr}, full", nP, p1, p2);
  hJetPhiVsEta[1] = new TH2D("hJetPhiVsEtaF", "Jet #varphi vs. #eta, full", nH, h1, h2, nF, f1, f2);
  hJetQt[0]       = new TH1D("hJetQtF", "Jet q_{T}, full (normalization different!)", nQ, q1, q2);
  hJetDr[0]       = new TH1D("hJetDrF", "Jet #Deltar, full (normalization different!)", nR, r1, r2);
  hJetS[0]        = new TH1D("hJetSf", "Jet s=A_{full}/A_{charge}, full (normalization different!)", nS, s1, s2);
  hJetDp[0]       = new TH1D("hJetDpF", "Jet #Deltap_{T}=p_{T}^{full}-p_{T}^{charge}, full (normalization different!)", nD, d1, d2);
  hJetQtVsDr[0]   = new TH2D("hJetQtVsDrF", "Jet q_{T} vs. #Deltar, full (normalization different!); #Deltar; q_{T}", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[0]    = new TH2D("hJetSvsDrF", "Jet s vs. #Deltar, full (normalization different!); #Deltar; s", nR, r1, r2, nS, s1, s2);
  // candidate matches
  hJetArea[2]     = new TH1D("hJetAreaA", "Jet area, candidates", nA, a1, a2);
  hJetEta[2]      = new TH1D("hJetEtaA", "Jet eta, candidates", nH, h1, h2);
  hJetPhi[2]      = new TH1D("hJetPhiA", "Jet phi, candidates", nF, f1, f2);
  hJetPt[2]       = new TH1D("hJetPtA", "Jet p_{T}, candidates", nP, p1, p2);
  hJetPtCorr[2]   = new TH1D("hJetPtCorrA", "Jet p_{T}^{corr}, candidates", nP, p1, p2);
  hJetPhiVsEta[2] = new TH2D("hJetPhiVsEtaA", "Jet #varphi vs. #eta, candidates", nH, h1, h2, nF, f1, f2);
  hJetQt[1]       = new TH1D("hJetQtA", "Jet q_{T}, candidates", nQ, q1, q2);
  hJetDr[1]       = new TH1D("hJetDrA", "Jet #Deltar, candidates", nR, r1, r2);
  hJetS[1]        = new TH1D("hJetSa", "Jet s=A_{cand.}/A_{charge}, candidates", nS, s1, s2);
  hJetDp[1]       = new TH1D("hJetDpA", "Jet #Deltap_{T}=p_{T}^{cand.}-p_{T}^{charge}, candidates (normalization different!)", nD, d1, d2);
  hJetQtVsDr[1]   = new TH2D("hJetQtVsDrA", "Jet q_{T} vs. #Deltar, candidates; #Deltar; q_{T}", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[1]    = new TH2D("hJetSvsDrA", "Jet s vs. #Deltar, candidates; #Deltar; s", nR, r1, r2, nS, s1, s2);
  // matches
  hJetArea[3]     = new TH1D("hJetAreaM", "Jet area, matches", nA, a1, a2);
  hJetEta[3]      = new TH1D("hJetEtaM", "Jet eta, matches", nH, h1, h2);
  hJetPhi[3]      = new TH1D("hJetPhiM", "Jet phi, matches", nF, f1, f2);
  hJetPt[3]       = new TH1D("hJetPtM", "Jet p_{T}, matches", nP, p1, p2);
  hJetPtCorr[3]   = new TH1D("hJetPtCorrM", "Jet p_{T}^{corr}, matches", nP, p1, p2);
  hJetPhiVsEta[3] = new TH2D("hJetPhiVsEtaM", "Jet #varphi vs. #eta, matches", nH, h1, h2, nF, f1, f2);
  hJetQt[2]       = new TH1D("hJetQtM", "Jet q_{T}, matches", nQ, q1, q2);
  hJetDr[2]       = new TH1D("hJetDrM", "Jet #Deltar, matches", nR, r1, r2);
  hJetS[2]        = new TH1D("hJetSm", "Jet s=A_{match}/A_{charge}, matches", nS, s1, s2);
  hJetDp[2]       = new TH1D("hJetDpM", "Jet #Deltap_{T}=p_{T}^{match}-p_{T}^{charge}, matches (normalization different!)", nD, d1, d2);
  hJetQtVsDr[2]   = new TH2D("hJetQtVsDrM", "Jet q_{T} vs. #Deltar, matches", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[2]    = new TH2D("hJetSvsDrM", "Jet s vs. #Deltar, matches", nR, r1, r2, nS, s1, s2);
  // junk (charged jets that weren't matched)
  hJetArea[4]     = new TH1D("hJetAreaJ", "Jet area, junk", nA, a1, a2);
  hJetEta[4]      = new TH1D("hJetEtaJ", "Jet eta, junk", nH, h1, h2);
  hJetPhi[4]      = new TH1D("hJetPhiJ", "Jet phi, junk", nF, f1, f2);
  hJetPt[4]       = new TH1D("hJetPtJ", "Jet p_{T}, junk", nP, p1, p2);
  hJetPtCorr[4]   = new TH1D("hJetPtCorrJ", "Jet p_{T}^{corr}, junk", nP, p1, p2);
  hJetPhiVsEta[4] = new TH2D("hJetPhiVsEtaJ", "Jet #varphi vs. #eta, junk", nH, h1, h2, nF, f1, f2);
  hJetQt[3]       = new TH1D("hJetQtJ", "Jet q_{T}, junk (normalization different!)", nQ, q1, q2);
  hJetDr[3]       = new TH1D("hJetDrJ", "Jet #Deltar, junk (normalization different!)", nR, r1, r2);
  hJetS[3]        = new TH1D("hJetSj", "Jet s=A_{junk}/A_{charge}, junk (normalization different!)", nS, s1, s2);
  hJetDp[3]       = new TH1D("hJetDpJ", "Jet #Deltap_{T}=p_{T}^{junk}-p_{T}^{charge}, junk (normalization different!)", nD, d1, d2);
  hJetQtVsDr[3]   = new TH2D("hJetQtVsDrJ", "Jet q_{T} vs. #Deltar, junk (normalization different!)", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[3]    = new TH2D("hJetSvsDrJ", "Jet s vs. #Deltar, junk (normalization different!)", nR, r1, r2, nS, s1, s2);
  // mystery (jets w/ dR > Rjet and |qT-1|<.1)
  hJetArea[5]     = new TH1D("hJetAreaY", "Jet area, mystery", nA, a1, a2);
  hJetEta[5]      = new TH1D("hJetEtaY", "Jet eta, mystery", nH, h1, h2);
  hJetPhi[5]      = new TH1D("hJetPhiY", "Jet phi, mystery", nF, f1, f2);
  hJetPt[5]       = new TH1D("hJetPtY", "Jet p_{T}, mystery", nP, p1, p2);
  hJetPtCorr[5]   = new TH1D("hJetPtCorrY", "Jet p_{T}^{corr}, mystery", nP, p1, p2);
  hJetPhiVsEta[5] = new TH2D("hJetPhiVsEtaY", "Jet #varphi vs. #eta, mystery", nH, h1, h2, nF, f1, f2);
  hJetQt[4]       = new TH1D("hJetQtY", "Jet q_{T}, mystery", nQ, q1, q2);
  hJetDr[4]       = new TH1D("hJetDrY", "Jet #Deltar, mystery", nR, r1, r2);
  hJetS[4]        = new TH1D("hJetSy", "Jet s=A_{?}/A_{charge}, mystery", nS, s1, s2);
  hJetDp[4]       = new TH1D("hJetDpY", "Jet #Deltap_{T}=p_{T}^{?}-p_{T}^{charge}, mystery (normalization different!)", nD, d1, d2);
  hJetQtVsDr[4]   = new TH2D("hJetQtVsDrY", "Jet q_{T} vs. #Deltar, mystery", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[4]    = new TH2D("hJetSvsDrY", "Jet s vs. #Deltar, mystery", nR, r1, r2, nS, s1, s2);
  // not matches (charged jets that weren't matched)
  hJetArea[6]     = new TH1D("hJetAreaN", "Jet area, (charge) not matches", nA, a1, a2);
  hJetEta[6]      = new TH1D("hJetEtaN", "Jet eta, (charge) not matches", nH, h1, h2);
  hJetPhi[6]      = new TH1D("hJetPhiN", "Jet phi, (charge) not matches", nF, f1, f2);
  hJetPt[6]       = new TH1D("hJetPtN", "Jet p_{T}, (charge) not matches", nP, p1, p2);
  hJetPtCorr[6]   = new TH1D("hJetPtCorrN", "Jet p_{T}^{corr}, (charge) not matches", nP, p1, p2);
  hJetPhiVsEta[6] = new TH2D("hJetPhiVsEtaN", "Jet #varphi vs. #eta, (charge) not matches", nH, h1, h2, nF, f1, f2);

  // errors
  hEfficiencyA   -> Sumw2();
  hEfficiencyPt  -> Sumw2();
  hResponseA     -> Sumw2();
  hResponseAn    -> Sumw2();
  hResponsePt    -> Sumw2();
  hResponsePtN   -> Sumw2();
  hResponsePtc   -> Sumw2();
  hResponsePtcN  -> Sumw2();
  hRefmultC      -> Sumw2();
  hRefmultF      -> Sumw2();
  hNumJetsC      -> Sumw2();
  hNumJetsF      -> Sumw2();
  hNumToMatch    -> Sumw2();
  hNumMatched    -> Sumw2();
  hNumHardC      -> Sumw2();
  hNumHardF      -> Sumw2();
  hChrgArea      -> Sumw2();
  hMatchArea     -> Sumw2();
  hChrgPtCorr    -> Sumw2();
  hMatchPtCorr   -> Sumw2();
  for (Int_t i = 0; i < nJetTypes; i++) {
    hJetArea[i]   -> Sumw2();
    hJetEta[i]    -> Sumw2();
    hJetPhi[i]    -> Sumw2();
    hJetPt[i]     -> Sumw2();
    hJetPtCorr[i] -> Sumw2();
  }
  for (Int_t i = 0; i < nMatchTypes; i++) {
    hJetQt[i] -> Sumw2();
    hJetDr[i] -> Sumw2();
    hJetS[i]  -> Sumw2();
    hJetDp[i] -> Sumw2();
  }


  // make sure both trees have the same no. of events
  Int_t nEvts = 0;
  Int_t cEvts = (Int_t) tChrg -> GetEntries();
  Int_t fEvts = (Int_t) tFull -> GetEntries();
  if (cEvts != fEvts) {
    cerr << "PANIC: Charged and full files do NOT have the same no. of events!\n"
         << "       nChargedEvt = " << cEvts << ", nFullEvt = " << fEvts
         << endl;
    assert(cEvts == fEvts);
  }
  else {
    nEvts = cEvts;
    cout << "  Beginning event loop..." << endl;
  }


  // vector for matching
  vector<Int_t> matchIndices;

  // event loop
  Int_t nBytesC  = 0;
  Int_t nBytesF  = 0;
  Int_t breakVal = 0;
  for (Int_t i = 0; i < nEvts; i++) {

    Int_t cBytes = tChrg -> GetEntry(i);
    Int_t fBytes = tFull -> GetEntry(i);
    if (cBytes < 0) {
      cerr << "ERROR: problem with charged event " << i << "...\n" << endl;
      breakVal = 1;
      break;
    }
    if (fBytes < 0) {
      cerr << "ERROR: problem with full event " << i << "..." << endl;
      breakVal = 1;
      break;
    }

    // should be same event
    Bool_t isSameEvent = (cEventIndex == fEventIndex);
    if (!isSameEvent) {
      cerr << "PANIC: event index is NOT the same! Stopped at i = " << i << "\n"
           << "       ChargedEvt = " << cEventIndex << ", FullEvt = " << fEventIndex << "\n"
           << endl;
      assert(isSameEvent);
    }

    nBytesC += cBytes;
    nBytesF += fBytes;
    if (!inBatchMode) {
      cout << "    Processing event " << i+1 << "/" << nEvts << "...\r" << flush;
      if (i+1 == nEvts) cout << endl;
    }
    else 
      cout << "    Processing event " << i+1 << "/" << nEvts << "..." << endl;


    // charged jet loop
    Int_t nHardC   = 0;
    Int_t nToMatch = 0;
    Int_t nMatched = 0;
    Int_t nCjets   = (Int_t) cJetEta -> size();
    Int_t nFjets   = (Int_t) fJetEta -> size();
    for (Int_t j = 0; j < nCjets; j++) {

      Double_t cA   = cJetArea -> at(j);
      Double_t cF   = cJetPhi  -> at(j);
      Double_t cH   = cJetEta  -> at(j);
      Double_t cPt  = cJetPt   -> at(j);
      Double_t cPtc = cPt - (cRho * cA);
      if (cPt > HardCut)
        ++nHardC;

      if (cPt < MinJetPt)
        continue;
      if (cA < MinArea)
        continue;

      // fill charged histograms
      hJetArea[0]     -> Fill(cA);
      hJetEta[0]      -> Fill(cH);
      hJetPhi[0]      -> Fill(cF);
      hJetPt[0]       -> Fill(cPt);
      hJetPtCorr[0]   -> Fill(cPtc);
      hJetPhiVsEta[0] -> Fill(cH, cF);


      if (cPt < Ccut)
        continue;
      else {
        hChrgArea   -> Fill(cA);
        hChrgPtCorr -> Fill(cPtc);
        ++nToMatch;
      }

      // match jets ['b' for best]
      Int_t    bIndex = 0.;
      Double_t bH     = 0.;
      Double_t bPt    = 0.;
      Double_t bPtc   = 0.;
      Double_t bF     = 0.;
      Double_t bA     = 0.;
      Double_t bQt    = 0.;
      Double_t bS     = 0.;
      Double_t bDp    = 0.;
      Double_t bDq    = 0.;
      Double_t bDr    = 999.;

      // full jet loop
      Bool_t isMatched = false;
      for (Int_t k = 0; k < nFjets; k++) {

        Double_t fA   = fJetArea -> at(k);
        Double_t fF   = fJetPhi  -> at(k);
        Double_t fH   = fJetEta  -> at(k);
        Double_t fPt  = fJetPt   -> at(k);
        Double_t fPtc = fPt - (fRho * fA);

        Bool_t isInAcceptance = true;
        if ((fPt < MinJetPt) || (fA < MinArea))
          isInAcceptance = false;


        // match jets
        Double_t qT = fPt / cPt;
        Double_t s  = fA / cA;
        Double_t dP = fPt - cPt;
        Double_t dQ = dP / cPt;
        Double_t dH = fH - cH;
        Double_t dF = fF - cF;
        Double_t dR = sqrt(dH*dH + dF*dF);

        if (fPt < Fcut)
          continue;

        // fill full histograms
        hJetQt[0]     -> Fill(qT);
        hJetDr[0]     -> Fill(dR);
        hJetS[0]      -> Fill(s);
        hJetDp[0]     -> Fill(dP);
        hJetQtVsDr[0] -> Fill(dR, qT);
        hJetSvsDr[0]  -> Fill(dR, s);


        Bool_t isBetter = false;
        Bool_t isInRcut = (dR < Rcut);
        //Bool_t isInQcut = ((qT > Qmin) && (qT < Qmax));
        if (isInRcut && isInAcceptance) {
          isMatched = true;
          isBetter  = ((dR < bDr) && (qT > bQt));

          // fill candidate histograms
          hJetArea[2]     -> Fill(fA);
          hJetEta[2]      -> Fill(fH);
          hJetPhi[2]      -> Fill(fF);
          hJetPt[2]       -> Fill(fPt);
          hJetPtCorr[2]   -> Fill(fPtc);
          hJetPhiVsEta[2] -> Fill(fH, fF);
          hJetQt[1]       -> Fill(qT);
          hJetDr[1]       -> Fill(dR);
          hJetS[1]        -> Fill(s);
          hJetDp[1]       -> Fill(dP);
          hJetQtVsDr[1]   -> Fill(dR, qT);
          hJetSvsDr[1]    -> Fill(dR, s);
        }
        else {
          // fill junk histograms
          hJetQt[3]     -> Fill(qT);
          hJetDr[3]     -> Fill(dR);
          hJetS[3]      -> Fill(s);
          hJetDp[3]     -> Fill(dP);
          hJetQtVsDr[3] -> Fill(dR, qT);
          hJetSvsDr[3]  -> Fill(dR, s);
        }

        // fill mystery histograms
        Double_t qCut      = TMath::Abs(qT - 1);
        Bool_t   isNearOne = (qCut < 0.1);
        if (!isInRcut && isNearOne) {
          hJetArea[5]     -> Fill(fA);
          hJetEta[5]      -> Fill(fH);
          hJetPhi[5]      -> Fill(fF);
          hJetPt[5]       -> Fill(fPt);
          hJetPtCorr[5]   -> Fill(fPtc);
          hJetPhiVsEta[5] -> Fill(fH, fF);
          hJetQt[4]       -> Fill(qT);
          hJetDr[4]       -> Fill(dR);
          hJetS[4]        -> Fill(s);
          hJetDp[4]       -> Fill(dP);
          hJetQtVsDr[4]   -> Fill(dR, qT);
          hJetSvsDr[4]    -> Fill(dR, s);
        }

        // check if candidate is best match
        if (isMatched && isBetter) {
          bIndex = k;
          bA     = fA;
          bH     = fH;
          bF     = fF;
          bPt    = fPt;
          bPtc   = fPtc;
          bQt    = qT;
          bS     = s;
          bDp    = dP;
          bDq    = dQ;
          bDr    = dR;
        }

      }  // end full jet loop


      // fill match histograms
      if (isMatched) {
        hResponseA      -> Fill(bA, cA);
        hResponseAn     -> Fill(bA, cA);
        hResponsePt     -> Fill(bPt, cPt);
        hResponsePtN    -> Fill(bPt, cPt);
        hResponsePtc    -> Fill(bPtc, cPtc);
        hResponsePtcN   -> Fill(bPtc, cPtc);
        hMatchArea      -> Fill(bA);
        hMatchPtCorr    -> Fill(bPtc);
        hJetArea[3]     -> Fill(bA);
        hJetEta[3]      -> Fill(bH);
        hJetPhi[3]      -> Fill(bF);
        hJetPt[3]       -> Fill(bPt);
        hJetPtCorr[3]   -> Fill(bPtc);
        hJetPhiVsEta[3] -> Fill(bH, bF);
        hJetQt[2]       -> Fill(bQt);
        hJetDr[2]       -> Fill(bDr);
        hJetS[2]        -> Fill(bS);
        hJetDp[2]       -> Fill(bDp);
        hJetQtVsDr[2]   -> Fill(bDr, bQt);
        hJetSvsDr[2]    -> Fill(bDr, bS);
        matchIndices.push_back(bIndex);
        ++nMatched;
      }
      else {
        hJetArea[6]     -> Fill(cA);
        hJetEta[6]      -> Fill(cH);
        hJetPhi[6]      -> Fill(cF);
        hJetPt[6]       -> Fill(cPt);
        hJetPtCorr[6]   -> Fill(cPtc);
        hJetPhiVsEta[6] -> Fill(cH, cF);
      }

    }  // end charged jet loop

    Int_t matchSize = (Int_t) matchIndices.size();
    if (nMatched != matchSize) {
      cerr << "ERROR: matchIndices did something weird in event " << i << "..." << endl;
      breakVal = 1;
      break;
    }


    // full jet loop
    Int_t nHardF = 0;
    for (Int_t j = 0; j < nFjets; j++) {

      Double_t fA   = fJetArea   -> at(j);
      Double_t fF   = fJetPhi    -> at(j);
      Double_t fH   = fJetEta    -> at(j);
      Double_t fPt  = fJetPt     -> at(j);
      Double_t fPtc = fPt - (fRho * fA);
      if (fPt > HardCut)
        ++nHardF;

      if (fPt < MinJetPt)
        continue;
      if (fA < MinArea)
        continue;

      // fill fullhistograms
      hJetArea[1]     -> Fill(fA);
      hJetEta[1]      -> Fill(fH);
      hJetPhi[1]      -> Fill(fF);
      hJetPt[1]       -> Fill(fPt);
      hJetPtCorr[1]   -> Fill(fPtc);
      hJetPhiVsEta[1] -> Fill(fH, fF);


      // check if full jet matches charged jet
      Bool_t isMatch = false;
      for (Int_t k = 0; k < nMatched; k++) {
        Int_t m = matchIndices.at(k);
        if (j == m) {
          isMatch = true;
          break;
        }
      }

      // fill junk histograms
      if (!isMatch) {
        hJetArea[4]     -> Fill(fA);
        hJetEta[4]      -> Fill(fH);
        hJetPhi[4]      -> Fill(fF);
        hJetPt[4]       -> Fill(fPt);
        hJetPtCorr[4]   -> Fill(fPtc);
        hJetPhiVsEta[4] -> Fill(fH, fF);
      }

    }  // end full jet loop


    // fill event histograms
    hRefmultC   -> Fill(cRefmult);
    hRefmultF   -> Fill(fRefmult);
    hNumJetsC   -> Fill(cNJets);
    hNumJetsF   -> Fill(fNJets);
    hNumToMatch -> Fill(nToMatch);
    hNumMatched -> Fill(nMatched);
    hNumHardC   -> Fill(nHardC);
    hNumHardF   -> Fill(nHardF);


    matchIndices.clear();

  }  // end event loop


  if (breakVal == 1) {
    cerr << "ERROR: Occured during event loop!\n"
         << "       Aborting program!"
         << endl;
    cassert(0);
  }
  else
    cout << "  Event loop finished!" << endl;


  // calculate efficiency
  cout << "  Calculating efficiency..." << endl;
  TH1D *hChargeA  = (TH1D*) hJetArea[0] -> Clone();
  TH1D *hMatchA   = (TH1D*) hJetArea[3] -> Clone();
  TH1D *hChargePt = (TH1D*) hJetPt[0]   -> Clone();
  TH1D *hMatchPt  = (TH1D*) hJetPt[3]   -> Clone();
  hEfficiencyA  -> Divide(hMatchA, hChargeA, 1., 1.);
  hEfficiencyPt -> Divide(hMatchPt, hChargePt, 1., 1.);


  cout << "  Normalizing histograms..." << endl;

  // bin widths
  const Double_t mBin   = (m2 - m1) / nM;
  const Double_t nBin   = (n2 - n1) / nN;
  const Double_t aBin   = (a2 - a1) / nA;
  const Double_t hBin   = (h2 - h1) / nH;
  const Double_t fBin   = (f2 - f1) / nF;
  const Double_t pBin   = (p2 - p1) / nP;
  const Double_t qBin   = (q2 - q1) / nQ;
  const Double_t rBin   = (r2 - r1) / nR;
  const Double_t sBin   = (s2 - s1) / nS;
  const Double_t dBin   = (d2 - d1) / nD;
  // overall normalizations
  const Double_t mNorm  = nEvts * mBin;
  const Double_t nNorm  = nEvts * nBin;
  const Double_t aNorm  = nEvts * aBin;
  const Double_t hNorm  = nEvts * hBin;
  const Double_t fNorm  = nEvts * fBin;
  const Double_t pNorm  = nEvts * pBin;
  const Double_t qNorm  = nEvts * qBin;
  const Double_t rNorm  = nEvts * rBin;
  const Double_t sNorm  = nEvts * sBin;
  const Double_t dNorm  = nEvts * dBin;
  const Double_t hfNorm = hNorm * fBin;
  // normalize histograms
  hRefmultC   -> Scale(1. / mNorm);
  hRefmultF   -> Scale(1. / mNorm);
  hNumJetsC   -> Scale(1. / nNorm);
  hNumJetsF   -> Scale(1. / nNorm);
  hNumToMatch -> Scale(1. / nNorm);
  hNumMatched -> Scale(1. / nNorm);
  hNumHardC   -> Scale(1. / nNorm);
  hNumHardF   -> Scale(1. / nNorm);
  for (Int_t i = 0; i < nJetTypes; i++) {
    hJetArea[i]     -> Scale(1. / aNorm);
    hJetEta[i]      -> Scale(1. / hNorm);
    hJetPhi[i]      -> Scale(1. / fNorm);
    hJetPt[i]       -> Scale(1. / pNorm);
    hJetPtCorr[i]   -> Scale(1. / pNorm);
    hJetPhiVsEta[i] -> Scale(1. / hfNorm);
  }

  // normalize response matrices
  const Int_t nAbinsX = hResponseAn   -> GetNbinsX();
  const Int_t nAbinsY = hResponseAn   -> GetNbinsY();
  const Int_t nPbinsX = hResponsePtN  -> GetNbinsX();
  const Int_t nPbinsY = hResponsePtN  -> GetNbinsY();
  const Int_t nCbinsX = hResponsePtcN -> GetNbinsX();
  const Int_t nCbinsY = hResponsePtcN -> GetNbinsY();
  for (Int_t i = 1; i < nAbinsY+1; i++) {
    const Double_t aNorm = hResponseAn -> Integral(1, nAbinsX, i, i);
    if (aNorm == 0.) continue;

    for (Int_t j = 1; j < nAbinsX+1; j++) {
      const Double_t oldCnt = hResponseAn -> GetBinContent(i, j);
      const Double_t oldErr = hResponseAn -> GetBinError(i, j);
      const Double_t newCnt = oldCnt / aNorm;
      const Double_t newErr = oldErr / aNorm;
      hResponseAn -> SetBinContent(i, j, newCnt);
      hResponseAn -> SetBinError(i, j, newErr);
    }
  }
  for (Int_t i = 1; i < nPbinsY+1; i++) {
    const Double_t pNorm = hResponsePtN -> Integral(1, nPbinsX, i, i);
    if (pNorm == 0.) continue;

    for (Int_t j = 1; j < nPbinsX+1; j++) {
      const Double_t oldCnt = hResponsePtN -> GetBinContent(i, j);
      const Double_t oldErr = hResponsePtN -> GetBinError(i, j);
      const Double_t newCnt = oldCnt / pNorm;
      const Double_t newErr = oldErr / pNorm;
      hResponsePtN -> SetBinContent(i, j, newCnt);
      hResponsePtN -> SetBinError(i, j, newErr);
    }
  }
  for (Int_t i = 1; i < nCbinsY+1; i++) {
    const Double_t cNorm = hResponsePtcN -> Integral(1, nCbinsX, i, i);
    if (cNorm == 0.) continue;

    for (Int_t j = 1; j < nCbinsX+1; j++) {
      const Double_t oldCnt = hResponsePtcN -> GetBinContent(i, j);
      const Double_t oldErr = hResponsePtcN -> GetBinError(i, j);
      const Double_t newCnt = oldCnt / cNorm;
      const Double_t newErr = oldErr / cNorm;
      hResponsePtcN -> SetBinContent(i, j, newCnt);
      hResponsePtcN -> SetBinError(i, j, newErr);
    }
  }


  // create directory structure
  const Int_t nDir = nJetTypes + 1;
  TDirectory *dir[nDir];
  dir[7] = (TDirectory*) fOut -> mkdir("EventInfo");
  dir[0] = (TDirectory*) fOut -> mkdir("ChrgJets");
  dir[1] = (TDirectory*) fOut -> mkdir("FullJets");
  dir[2] = (TDirectory*) fOut -> mkdir("Candidates");
  dir[3] = (TDirectory*) fOut -> mkdir("MatchJets");
  dir[4] = (TDirectory*) fOut -> mkdir("JunkJets");
  dir[5] = (TDirectory*) fOut -> mkdir("Mystery");
  dir[6] = (TDirectory*) fOut -> mkdir("NotMatches");

  // write and close output
  fOut          -> cd();
  hEfficiencyA  -> Write();
  hEfficiencyPt -> Write();
  hResponseA    -> Write();
  hResponseAn   -> Write();
  hResponsePt   -> Write();
  hResponsePtN  -> Write();
  hResponsePtc  -> Write();
  hResponsePtcN -> Write();
  dir[7]        -> cd();
  hRefmultC     -> Write();
  hRefmultF     -> Write();
  hNumJetsC     -> Write();
  hNumJetsF     -> Write();
  hNumToMatch   -> Write();
  hNumMatched   -> Write();
  hNumHardC     -> Write();
  hNumHardF     -> Write();
  hChrgArea     -> Write();
  hMatchArea    -> Write();
  hChrgPtCorr   -> Write();
  hMatchPtCorr  -> Write();
  for (Int_t i = 0; i < nJetTypes; i++) {
    dir[i]          -> cd();
    hJetArea[i]     -> Write();
    hJetEta[i]      -> Write();
    hJetPhi[i]      -> Write();
    hJetPt[i]       -> Write();
    hJetPtCorr[i]   -> Write();
    hJetPhiVsEta[i] -> Write();
  }
  for (Int_t i = 0; i < nMatchTypes; i++) {
    Int_t iDir = i + 1;
    dir[iDir]     -> cd();
    hJetQt[i]     -> Write();
    hJetDr[i]     -> Write();
    hJetS[i]      -> Write();
    hJetDp[i]     -> Write();
    hJetQtVsDr[i] -> Write();
    hJetSvsDr[i]  -> Write();
  }
  fOut -> cd();
  fOut -> Close();

  // close input
  fChrg -> cd();
  fChrg -> Close();
  fFull -> cd();
  fFull -> Close();


  cout << "Matching script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
