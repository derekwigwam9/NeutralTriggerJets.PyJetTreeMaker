// 'MergeData.C'
// Derek Anderson
// 10.24.2016
//
// This macro merges numerous TTree's into a TChain.

#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;


void MergeData() {

  const Int_t   nFiles = 3675;
  const TString oName("py23pi01par.merged.root");
  const TString tName("ParTree");
  const TString cName("ParTree");
  const TString fList("pythia23pi0.set1.list");
  const TString fBad("pythia23pi0.bad1.list");

  cout << "\n  Merging files..." << endl;


  TFile  *oFile = new TFile(oName, "recreate");
  TChain *chain = new TChain(tName, cName);

  ifstream files(fList.Data());
  ofstream bad(fBad.Data());

  Int_t  add;
  UInt_t iFile(0);
  string file;
  while (files) {
    files >> file;
    TString name(file);
    add = chain -> Add(name.Data(), 0);
    if (add == 0) {
      cout << "    File '" << name.Data() << "' no good!" << endl;
      bad  << name.Data();
      bad  << endl;
    }
    else
      cout << "    File '" << name.Data() << "' merged." << endl;

    iFile++;
    if (iFile == nFiles) break;
  }

  cout << "  " << iFile << " files merged!\n" << endl;


  // draw a few plots to check
  TCanvas *cTotalMult = new TCanvas("cTotalMult", "Total multiplicity", 200, 10, 700, 500);
  cTotalMult -> SetGrid(0, 0);
  chain      -> Draw("Events_refmult");
  cTotalMult -> Write();
  cTotalMult -> Close();

  TCanvas *cTrackPt = new TCanvas("cTrackPt", "Track pT", 200, 10, 700, 500);
  cTrackPt -> SetGrid(0, 0);
  cTrackPt -> SetLogy(1);
  chain    -> Draw("pTracks_pT");
  cTrackPt -> Write();
  cTrackPt -> Close(); 


  oFile -> cd();
  chain -> Write();
  oFile -> Close();

}
