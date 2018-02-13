// 'MakeRhoComparison.C'
// Derek Anderson
//
// This plots the rho distributions of several sets of events against
// each other.

#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;



void MakeRhoComparison() {

  cout << "\nPlotting Rho distributions..." << endl;


  TFile *oFile  = new TFile("RhoPlots.root", "recreate");
  TFile *pFile1 = new TFile("Pythia17p.r05a065rm1.pP.Sep27.root");
  TFile *pFile2 = new TFile("Pythia17p.r05a065rm2.pP.Sep27.root");
  TFile *pFile3 = new TFile("Pythia17p.r05a065rm3.pP.Sep27.root");
  TFile *nFile1 = new TFile("Pythia17p.r05a065rm1.p0.Sep27.root");
  TFile *nFile2 = new TFile("Pythia17p.r05a065rm2.p0.Sep27.root");
  TFile *nFile3 = new TFile("Pythia17p.r05a065rm3.p0.Sep27.root");

  cout << "  Files loaded..." << endl;


  oFile -> cd();
  TH1D *hRhoP[3];
  TH1D *hRhoN[3];

  pFile1 -> cd();
  hRhoP[0] = (TH1D*) pFile1 -> Get("QA/hRho");
  pFile2 -> cd();
  hRhoP[1] = (TH1D*) pFile2 -> Get("QA/hRho");
  pFile3 -> cd();
  hRhoP[2] = (TH1D*) pFile3 -> Get("QA/hRho");
  nFile1 -> cd();
  hRhoN[0] = (TH1D*) nFile1 -> Get("QA/hRho");
  nFile2 -> cd();
  hRhoN[1] = (TH1D*) nFile2 -> Get("QA/hRho");
  nFile3 -> cd();
  hRhoN[2] = (TH1D*) nFile3 -> Get("QA/hRho");

  cout << "  Histograms grabbed:" << endl;
  for (Int_t i = 0; i < 3; i++) {
    cout << "    hRhoP[" << i << "] = " << hRhoP[i]
         << ", hRhoN[" << i << "] = " << hRhoN[i]
         << endl;
  }


  hRhoP[0] -> SetMarkerColor(kRed);
  hRhoP[0] -> SetLineColor(kRed);
  hRhoP[1] -> SetMarkerColor(kGreen);
  hRhoP[1] -> SetLineColor(kGreen);
  hRhoP[2] -> SetMarkerColor(kBlue);
  hRhoP[2] -> SetLineColor(kBlue);

}

// End ------------------------------------------------------------------------
