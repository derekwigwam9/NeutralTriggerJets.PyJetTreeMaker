
#include "./Nihar_functions.h"
const Int_t NbinPt = 100;//100
const Int_t NbinPt_min = -20;
const Int_t NbinPt_max = 100; 
const Int_t Max_Norm_xBin = 0;
const Int_t Min_Norm_xBin = -10;

const Int_t JetR= 3;
TString JetRadiusL = "R= 0.3";

const Int_t XNbinPt= 75;
const Int_t   XNbinPt_min = -5;
const Int_t   XNbinPt_max= 43;

void PlotPythiaRatios()
{
 

 
  TFile *fOut    = new TFile("pythiaRatios.r03a02rm1.d28m1y2017.root", "recreate");
  TFile *fFuncDC = new TFile("pythiaCheckD.r03a02rm1chrg.d28m1y2017.root");
  TFile *fFuncDF = new TFile("pythiaCheckD.r03a02rm1full.d28m1y2017.root");
  TFile *fFuncRC = new TFile("pythiaCheckR.p40.r03a02rm1chrg.d28m1y2017.root");
  TFile *fFuncRF = new TFile("pythiaCheckR.p40.r03a02rm1full.d28m1y2017.root");
  TH1D  *fPythDC = (TH1D*) fFuncDC -> Get("h_Ratio");
  TH1D  *fPythDF = (TH1D*) fFuncDF -> Get("h_Ratio");
  TH1D  *fPythRC = (TH1D*) fFuncRC -> Get("h_Ratio");
  TH1D  *fPythRF = (TH1D*) fFuncRF -> Get("h_Ratio");
  fPythDC -> SetNameTitle("fDirectCharged", "");
  fPythDF -> SetNameTitle("fDirectFull", "");
  fPythRC -> SetNameTitle("fRichChargedPP", "");
  fPythRF -> SetNameTitle("fRichFullPP", "");

  const Bool_t addAuAu = false;
  TH1D *fPythAC;
  TH1D *fPythAF;
  if (addAuAu) {
    TFile *fFuncAC = new TFile("rich/pythiaRichPlots.AA.r05a065rm1chrg.d26m1y2017.root");
    TFile *fFuncAF = new TFile("rich/pythiaRichPlots.AA.r05a065rm1full.d26m1y2017.root");
    fPythAC = (TH1D*) fFuncAC -> Get("h_Ratio");
    fPythAF = (TH1D*) fFuncAF -> Get("h_Ratio");
    fPythAC -> SetNameTitle("fRichChargedAA", "");
    fPythAF -> SetNameTitle("fRichFullAA", "");
  }

#if 0
  h_Datajet_JetPtCorr_wobgsub_g = (TH1F*)*f_g->Get("h_JetPtCorr_gamma_Recoil");
  h_Datajet_JetPtCorr_wobgsub_pi0 = (TH1F*)*f_pi->Get("h_JetPtCorr_pi0_Recoil");
  TH1F *h_pi0 =h_Datajet_JetPtCorr_wobgsub_pi0;
  TH1F *h_g =h_Datajet_JetPtCorr_wobgsub_g;
  TString nameStatus = "W/ Bg. Subtracted";
  Double_t Y_max= 2;
  
#endif


#if 1
 
  TString nameStatus = "";
  Double_t Y_max= 1.52;

#endif

   //_________PLOTTINE PART_______
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(1);
  gStyle->SetOptLogz(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetFrameLineWidth(2); 
  gStyle->SetLineColor(2);

  TCanvas *c = new TCanvas("c", "canvas",0,45,821,460);
  gStyle->SetOptStat(0);
  c->Range(-37.61978,-0.003397649,57.08655,0.01106914);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetGrid(0, 0);  // [Derek]
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.1860465);
  c->SetRightMargin(0.02203182);
  c->SetTopMargin(0.073903);
  c->SetBottomMargin(0.2286374);
  c->SetFrameLineWidth(2);
  c->SetFrameBorderMode(0);
  c->SetFrameLineWidth(2);
  c->SetFrameBorderMode(0);

  TH2F *h5 = new TH2F("h5","",XNbinPt,XNbinPt_min,XNbinPt_max,1000,0.2,8.7);
  h5->SetStats(0);

  Int_t ci;      // for color index setting
  TColor *color; // for color definition with alpha
  ci = TColor::GetColor("#000099");
  h5->SetLineColor(ci);
  h5->GetXaxis()->SetTitle("p_{T,jet}^{par} [GeV/c]");
  h5->GetXaxis()->CenterTitle(true);
  h5->GetXaxis()->SetLabelFont(42);
  h5->GetXaxis()->SetLabelOffset(0.001);
  h5->GetXaxis()->SetLabelSize(0.08);
  h5->GetXaxis()->SetTitleSize(0.08);
  h5->GetXaxis()->SetTickLength(0.06);
  h5->GetXaxis()->SetTitleOffset(1.14);
  h5->GetXaxis()->SetTitleFont(42);
  h5->GetYaxis()->SetTitle("#pi^{0}_{rich}/#gamma_{rich,dir}");
  h5->GetYaxis()->CenterTitle(true);
  h5->GetYaxis()->SetNdivisions(505);
  h5->GetYaxis()->SetLabelFont(42);
  h5->GetYaxis()->SetLabelSize(0.08);
  h5->GetYaxis()->SetTitleSize(0.1);
  h5->GetYaxis()->SetTickLength(0.04);
  h5->GetYaxis()->SetTitleOffset(0.94);
  h5->GetYaxis()->SetTitleFont(42);
  h5->GetYaxis()->SetNdivisions(8, 5, 0);
  h5->GetZaxis()->SetLabelFont(42);
  h5->GetZaxis()->SetLabelSize(0.035);
  h5->GetZaxis()->SetTitleSize(0.035);
  h5->GetZaxis()->SetTitleFont(42);
  h5->Draw("");

  Double_t lcolor;
  Double_t Mstyle;
  Double_t Msize;

  // set styles
  fPythDC -> SetLineStyle(1);
  fPythDC -> SetLineWidth(2);
  fPythDC -> SetFillStyle(3001);
  fPythDC -> SetMarkerStyle(1);
  fPythDF -> SetLineStyle(1);
  fPythDF -> SetLineWidth(2);
  fPythDF -> SetFillStyle(3001);
  fPythDF -> SetMarkerStyle(1);
  if (addAuAu) {
    fPythAC -> SetLineStyle(1);
    fPythAC -> SetLineWidth(2);
    fPythAC -> SetFillStyle(3001);
    fPythAC -> SetMarkerStyle(1);
    fPythAF -> SetLineStyle(1);
    fPythAF -> SetLineWidth(2);
    fPythAF -> SetFillStyle(3001);
    fPythAF -> SetMarkerStyle(1);
  }
  fPythRC -> SetLineStyle(1);
  fPythRC -> SetLineWidth(2);
  fPythRC -> SetFillStyle(3001);
  fPythRC -> SetMarkerStyle(1); 
  fPythRF -> SetLineStyle(1);
  fPythRF -> SetLineWidth(2);
  fPythRF -> SetFillStyle(3001);
  fPythRF -> SetMarkerStyle(1);
  // set colors
  fPythDC -> SetLineColor(kRed);
  fPythDC -> SetFillColor(kRed);
  fPythDC -> SetMarkerColor(kRed);
  fPythDF -> SetLineColor(kGreen);
  fPythDF -> SetFillColor(kGreen);
  fPythDF -> SetMarkerColor(kGreen);
  if (addAuAu) {
    fPythDF -> SetLineColor(kOrange);
    fPythDF -> SetFillColor(kOrange);
    fPythDF -> SetMarkerColor(kOrange);
    fPythAC -> SetLineColor(kYellow);
    fPythAC -> SetFillColor(kYellow);
    fPythAC -> SetMarkerColor(kYellow);
    fPythAF -> SetLineColor(kGreen);
    fPythAF -> SetFillColor(kGreen);
    fPythAF -> SetMarkerColor(kGreen);
  }
  fPythRC -> SetLineColor(kBlue);
  fPythRC -> SetFillColor(kBlue);
  fPythRC -> SetMarkerColor(kBlue);
  fPythRF -> SetLineColor(kViolet);
  fPythRF -> SetFillColor(kViolet);
  fPythRF -> SetMarkerColor(kViolet);
  // draw ratios
  fPythDC -> Draw("same");
  fPythDF -> Draw("same");
  if (addAuAu) {
    fPythAC -> Draw("same");
    fPythAF -> Draw("same");
  }
  fPythRC -> Draw("same");
  fPythRF -> Draw("same");

  c -> Update();
  TLine *l = new TLine(XNbinPt_min,1,XNbinPt_max,1);
  l->SetLineStyle(2);
  l->SetLineColor(2);
  l->Draw();


  TLegend *leg2;
  TString JetRecoLegnd;
  JetRecoLegnd += "anti-k_{T} algo.,";
  JetRecoLegnd += JetRadiusL; 
  Draw_Legend(leg2, 0.54,0.73,0.66,0.97,"Pythia: p+p 200 GeV",JetRecoLegnd.Data(),"9 < p_{T}^{trig} < 30 GeV/c",nameStatus.Data(),"");
  
  leg = new TLegend(0.64,0.60,0.86,0.70);
  leg -> SetFillColor(0);
  leg -> SetTextSize(0.03);  
  leg -> SetLineColor(0);
  leg -> AddEntry(fPythDC,"#gamma_{dir} [charged]", "l");
  leg -> AddEntry(fPythDF,"#gamma_{dir} [full]", "l");
  if (addAuAu) {
    leg -> AddEntry(fPythAC, "#gamma_{rich}(70\045) [charged]", "l");
    leg -> AddEntry(fPythAF, "#gamma_{rich}(70\045) [full]", "l");
  }
  leg -> AddEntry(fPythRC,"#gamma_{rich}(40\045) [charged]", "l");
  leg -> AddEntry(fPythRF,"#gamma_{rich}(40\045) [full]", "l");
  leg -> Draw();


  // save output
  fOut    -> cd();
  fPythDC -> Write();
  fPythDF -> Write();
  fPythRC -> Write();
  fPythRF -> Write();
  if (addAuAu) {
    fPythAC -> Write();
    fPythAF -> Write();
  }
  c       -> Write();
  c       -> Close();
  fOut    -> Close();

}
