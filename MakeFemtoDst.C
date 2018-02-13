// 'MakeFemtoDst.C'
// Derek Anderson
// 08.05.2016
//
// This macro produces a compact tree of jets (the "FemtoDST") via the
// 'StFemtoDstMaker' class.  Use the switches 'particle' and 'detector'
// to turn on/off production of particle-level and detector-level jets
// respectively.
//
// Last updated: 05.19.2017

#include <TSystem>
#include <iostream>
#include "TString.h"

using namespace std;

class StFemtoDstMaker;


// i/o parameters
const TString ipFile("./input/Pythia23p.gMerged.root");
const TString idFile("./input/Pythia23d.gMerged.root");
// trigger parameters
const Int_t    tID     = 22;
const Double_t eTmin   = 9.;
const Double_t eTmax   = 20.;
const Double_t hTrkMax = 1.;
const Double_t hTrgMax = 0.9;
// jet parameters
const Int_t    nRM   = 1;
const Double_t rJet  = 0.3;
const Double_t aMin  = 0.2;
const Double_t pTmin = 0.2;
const Double_t pTmax = 20.;
const Double_t qMin  = 0.15;



void MakeFemtoDst(const Int_t nTrgs=-1, const Int_t StartEvt=0, const Int_t StopEvt=-1, const Bool_t particle, const Bool_t detector, const Int_t type=0, const TString oFile="oops.root") {

  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjet.so");
  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjettools.so");
  gSystem -> Load("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Ana_nutralTgr_Jet/JetMacro/.sl64_gcc482/lib/libStFemtoDstMaker.so");

  // lower verbosity
  gErrorIgnoreLevel = kError;


  // create particle-level jets
  if (particle) {
    StFemtoDstMaker p(ipFile, oFile, 0, type);
    p.Init(nRM, tID, rJet, aMin, pTmin, pTmax, qMin, eTmin, eTmax, hTrkMax, hTrgMax);
    p.Make(nTrgs, StartEvt, StopEvt);
    p.Finish();
  }

  // create detector-level jets
  if (detector) {
    StFemtoDstMaker d(idFile, oFile, 1, type);
    d.Init(nRM, tID, rJet, aMin, pTmin, pTmax, qMin, eTmin, eTmax, hTrkMax, hTrgMax);
    d.Make(nTrgs, StartEvt, StopEvt);
    d.Finish();
  }

}

// End ------------------------------------------------------------------------
