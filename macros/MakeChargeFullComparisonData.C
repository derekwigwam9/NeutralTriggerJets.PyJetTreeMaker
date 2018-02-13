
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

void MakeChargeFullComparisonData()
{
 


  // co-opted old notation:
  //   g  ==> charge
  //   pi ==> full 
  TFile *fOut = new TFile("dataCompP.r03a02rm1.d27m1y2017.root", "recreate");
  f_g= new TFile("data.r03a02rm1chrg.plots.d26m1y2017.root");
  f_pi= new TFile("data.r03a02rm1full.plots.d26m1y2017.root");

  // For comparing against data
  const Bool_t doComparison = false;
  TF1 *fPyth;
  if (doComparison) {
    TFile *fFunc = new TFile("pythiaPlots.r05a065rm1.full.root");
    fPyth = (TF1*) fFunc -> Get("f1");
    fPyth -> SetName("fPyth");
  }

  f_pi.ls();
  f_g.ls();
  
  
#if 0
  h_Datajet_JetPtCorr_wobgsub_g = (TH1F*)*f_g->Get("h_JetPtCorr_gamma_Recoil");
  h_Datajet_JetPtCorr_wobgsub_pi0 = (TH1F*)*f_pi->Get("h_JetPtCorr_pi0_Recoil");
  TH1F *h_pi0 =h_Datajet_JetPtCorr_wobgsub_pi0;
  TH1F *h_g =h_Datajet_JetPtCorr_wobgsub_g;
  TString nameStatus = "W/ Bg. Subtracted";
  Double_t Y_max= 2;
  
#endif


#if 1
  
   
  h_Datajet_JetPtCorr_g = (TH1F*)*f_g->Get("Pi0/hJetPtCorrP");
  h_Datajet_JetPtCorr_pi0 = (TH1F*)*f_pi->Get("Pi0/hJetPtCorrP");

  TH1F *h_pi0 =h_Datajet_JetPtCorr_pi0;
  TH1F *h_g =h_Datajet_JetPtCorr_g;
  TString nameStatus = "#pi^{0} trigger";
  Double_t Y_max= 1.52;

  //h_pi0 -> Rebin(2);
  h_pi0 -> SetNameTitle("hFull", "Full jets");
  //h_g   -> Rebin(2);
  h_g   -> SetNameTitle("hCharge", "Charge jets");


#endif


  
  //  h_Datajet_JetPtCorr_g->Draw("P");
  //  h_Datajet_JetPtCorr_pi0->Draw("PSAME");
  TH1D *h_Ratio = new TH1D("h_Ratio","",NbinPt,NbinPt_min,NbinPt_max);
  Int_t nbin_g=  h_g->GetNbinsX();
  Int_t nbin_pi=  h_pi0->GetNbinsX();
  if( nbin_g != nbin_pi){ cout<< "Bin Mismatch..... CHECK!!!!!"<<endl; break;}
  
  Int_t ratio_pnt=0;
  for(int i=1; i <=nbin_g ; i++)
    {

      double g_bincnter = h_g->GetBinCenter(i);
      double y_g = h_g->GetBinContent(i);
      double y_g_err = h_g->GetBinError(i);
      //      cout<<" 1...   "<<h_normalizedErr_full->GetBinContent(i)<<"  "<<h_normalizedErr_full->GetBinError(i)<<endl;
      
      double y_pi = h_pi0->GetBinContent(i);
      double y_pi_err = h_pi0->GetBinError(i);
      if(!(y_pi > 0 && y_g > 0 ))continue;

      double ratio = y_pi/y_g ;
      double ratio_err = sqrt(pow(ratio,2)*(pow((y_pi_err/y_pi),2) + pow((y_g_err/y_g),2)));

      cout<<i<<"  "<<h_g->GetBinCenter(i)<< "  ratio= "<<ratio<<"   err= "<<ratio_err<<endl;

      int     bin_me_max = h_Ratio->FindBin(g_bincnter);
      h_Ratio->SetBinContent(bin_me_max,ratio);
      h_Ratio->SetBinError(bin_me_max,ratio_err);

      ratio_pnt++;
    }
    h_Ratio->SetEntries(ratio_pnt);


    //_________Difference
    TH1D *h_Diff = new TH1D("h_Diff","",NbinPt,NbinPt_min,NbinPt_max);
    Int_t diff_pnt=0;
    for(int i=0; i <nbin_g ; i++)
    {
      
      double y_g = h_g->GetBinContent(i);
      double y_g_err = h_g->GetBinError(i);
      //      cout<<" 1...   "<<h_normalizedErr_full->GetBinContent(i)<<"  "<<h_normalizedErr_full->GetBinError(i)<<endl;
      
      double y_pi = h_pi0->GetBinContent(i);
      double y_pi_err = h_pi0->GetBinError(i);
      if(!(y_pi > 0 && y_g > 0 ))continue;

      double diff = y_g - y_pi ;
      double diff_err = sqrt(pow(y_pi_err,2) + pow(y_g_err,2));

      //      cout<<i<< "  Diff= "<<diff<<"   err= "<<diff_err<<endl;
      h_Diff->SetBinContent(i,diff);
      h_Diff->SetBinError(i,diff_err);

      diff_pnt++;
    }
    h_Diff->SetEntries(ratio_pnt);


   //_________PLOTTINE PART_______
   //  gROOT->Reset();
  //  gStyle->SetPadGridX(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(1);
  gStyle->SetOptLogz(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetPadGridY(1);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadTopMargin(0.05);
  //    gStyle->SetPadBorderSize(2);
  gStyle->SetFrameLineWidth(2); 
  gStyle->SetLineColor(2);
  
  //_________
  //   TCanvas *c = new TCanvas("c", "canvas",285,45,656,743);
  TCanvas *c = new TCanvas("c", "canvas",700,800);
  //   c->Range(-7.241379,-2.096007,5.249042,2.699528);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetGrid(0, 0);
   c->SetLogy();
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetLeftMargin(0.1794479);
   c->SetRightMargin(0.0993865);
   c->SetTopMargin(0.01536313);
   c->SetBottomMargin(0.3743017);
   c->SetFrameBorderMode(0);
   c->SetFrameBorderMode(0);
   
   //   TH1F *h1 = new TH1F("h1","",100,-5,5);
   TH2D *h1 = new TH2D("h1","",XNbinPt,XNbinPt_min,XNbinPt_max,1000,0.000004,Y_max);
   //____X-axis
   h1->GetXaxis()->SetLabelFont(42);
   h1->GetXaxis()->SetLabelSize(0.0);
   h1->GetXaxis()->SetTitleSize(0.033);
   h1->GetXaxis()->SetTitleFont(42);
   //____Y-axis
   h1->GetYaxis()->SetTitle("1/N_{trig} dN_{jets}/(dp_{T,jet}^{reco} d#eta_{jet})");
   h1->GetYaxis()->CenterTitle(true);
   h1->GetYaxis()->SetLabelFont(42);
   h1->GetYaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetTitleSize(0.05);
   h1->GetYaxis()->SetTitleOffset(1.63);
   h1->GetYaxis()->SetTitleFont(42);
   h1->GetZaxis()->SetLabelFont(42);
   h1->GetZaxis()->SetLabelSize(0.035);
   h1->GetZaxis()->SetTitleSize(0.035);
   h1->GetZaxis()->SetTitleFont(42);
   h1->Draw("");
   
   Double_t lcolor;
   Double_t Mstyle;
   Double_t Msize;
    
   //   Draw_hist(h_g, "PSAME", kOrange+4,kOrange+4, 20, 1.1 );
   Draw_hist(h_g, "PSAME", 4, 4, 29, 1.1 );
   //   Draw_hist(h_pi0, "PSAME", kGreen+4,kGreen+4, 29, 1.6 );
   Draw_hist(h_pi0, "PSAME", 4, 4, 30, 1.1 );
  
   

    TLegend *leg2;
    
    TString JetRecoLegnd;
    JetRecoLegnd += "anti-k_{T} algo.,";
    JetRecoLegnd += JetRadiusL;
    
    Draw_Legend(leg2, 0.54,0.73,0.66,0.97,"p+p 200 GeV",JetRecoLegnd.Data(),"9 < p_{T}^{trig} < 30 GeV/c",nameStatus.Data(),"#color[2]{STAR Preliminary}"  );
    
    leg = new TLegend(0.64,0.60,0.86,0.70);
    leg->SetFillColor(0);   leg->SetTextSize(0.03);  
    leg->SetLineColor(0);
    leg->AddEntry(h_g,"Charged jets","p");
    leg->AddEntry(h_pi0,"Full jets","p");
 
    leg->Draw();

      
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "pad2",0.02,0.007,0.917,0.377);
   pad2->Draw();
   pad2->cd();
   pad2->SetGrid(0, 0);  // [Derek]
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   //   pad2->SetLogy();
   pad2->SetLeftMargin(0.1784615);
   pad2->SetRightMargin(0.01692308);
   pad2->SetTopMargin(0.007434944);
   pad2->SetBottomMargin(0.3745725);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderMode(0);
   
   TH2F *h3 = new TH2F("h3","",XNbinPt,XNbinPt_min,XNbinPt_max,1000,0.2,8.7);
   //_____X-axis
   h3->GetXaxis()->SetTitle("p_{T,jet}^{reco} [GeV/c]");
   h3->GetXaxis()->CenterTitle(true);
   h3->GetXaxis()->SetLabelFont(42);
   h3->GetXaxis()->SetLabelOffset(0.001);
   h3->GetXaxis()->SetLabelSize(0.14);
   h3->GetXaxis()->SetTitleSize(0.14);
   h3->GetXaxis()->SetTitleFont(42);
   h3->GetXaxis()->SetTickLength(0.06);
   //_____X-axis
   h3->GetYaxis()->SetTickLength(0.04);
   h3->GetYaxis()->SetTitle("full/charged");
   h3->GetYaxis()->CenterTitle(true);
   h3->GetYaxis()->SetLabelFont(42);
   h3->GetYaxis()->SetLabelSize(0.12);
   h3->GetYaxis()->SetTitleSize(0.12);
   h3->GetYaxis()->SetTitleOffset(0.71);
   h3->GetYaxis()->SetTitleFont(42);
   h3->GetYaxis()->SetNdivisions(8, 5, 0);
   h3->GetZaxis()->SetLabelFont(42);
   h3->GetZaxis()->SetLabelSize(0.035);
   h3->GetZaxis()->SetTitleSize(0.035);
   h3->GetZaxis()->SetTitleFont(42);
   h3->Draw("");


   TLine *l = new TLine(XNbinPt_min,1,XNbinPt_max,1);
   l->SetLineStyle(2);
   l->SetLineColor(2);
   l->Draw();
   Draw_hist(h_Ratio, "PESAME", 1,1, 20, 1.2  );

   TF1 *f1 = new TF1("f1", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", XNbinPt_min+5, 25);
   f1 -> SetLineColor(kViolet);
   h_Ratio->Fit("f1", "R", "SAME");
   h_Ratio->SetStats(false);
   if (doComparison) {
     fPyth -> SetLineColor(kGreen);
     fPyth -> Draw("same");
   }

   cout<<"Paramters for fit function: "<<" p0: "<<f1->GetParameter(0)<<
     "  p1: "<<f1->GetParameter(1)<<
     "  p2: "<<f1->GetParameter(2)<<
     "  p3: "<<f1->GetParameter(3)<<endl;
     
   
   leg = new TLegend(0.64,0.60,0.86,0.70);
    leg->SetFillColor(0);   leg->SetTextSize(0.1);  
    leg->SetLineColor(0);
    leg->AddEntry(f1,"Pol. fit","l");
    if (doComparison) leg->AddEntry(fPyth,"PYTHIA", "l");
 
    leg->Draw();


   //___________Difference Plots
   //_________

  TCanvas *c1 = new TCanvas("c1", "canvas",0,45,821,460);
   gStyle->SetOptStat(0);
   c1->Range(-37.61978,-0.003397649,57.08655,0.01106914);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGrid(0, 0);  // [Derek]
   //   c1->SetLogz();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.1860465);
   c1->SetRightMargin(0.02203182);
   c1->SetTopMargin(0.073903);
   c1->SetBottomMargin(0.2286374);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);

   //   TH2F *h5 = new TH2F("h5","",75,-20,55,1000,-9e-08,0.01);
   //   TH2F *h5 = new TH2F("h5","",75,-20,55,1000,-9e-08,0.1);
   TH2F *h5 = new TH2F("h5","",75,-20,55,1000,-0.1,0.1);
   h5->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   h5->SetLineColor(ci);
   h5->GetXaxis()->SetTitle("p_{T,jet}^{reco} [GeV/c]");
   h5->GetXaxis()->CenterTitle(true);
   h5->GetXaxis()->SetLabelFont(42);
   h5->GetXaxis()->SetLabelOffset(0.001);
   h5->GetXaxis()->SetLabelSize(0.08);
   h5->GetXaxis()->SetTitleSize(0.08);
   h5->GetXaxis()->SetTickLength(0.06);
   h5->GetXaxis()->SetTitleOffset(1.14);
   h5->GetXaxis()->SetTitleFont(42);
   h5->GetYaxis()->SetTitle("charged - full");
   h5->GetYaxis()->CenterTitle(true);
   h5->GetYaxis()->SetNdivisions(505);
   h5->GetYaxis()->SetLabelFont(42);
   h5->GetYaxis()->SetLabelSize(0.08);
   h5->GetYaxis()->SetTitleSize(0.1);
   h5->GetYaxis()->SetTickLength(0.04);
   h5->GetYaxis()->SetTitleOffset(0.94);
   h5->GetYaxis()->SetTitleFont(42);
   h5->GetZaxis()->SetLabelFont(42);
   h5->GetZaxis()->SetLabelSize(0.035);
   h5->GetZaxis()->SetTitleSize(0.035);
   h5->GetZaxis()->SetTitleFont(42);
   h5->Draw("");
   TLine *line = new TLine(-20,1,55,1);
   line->SetLineStyle(2);
   line->Draw();

   TLine *l = new TLine(XNbinPt_min,1,XNbinPt_max,1);
   l->SetLineStyle(2);
   l->SetLineColor(1);
   l->Draw();
   Draw_hist(h_Diff, "PESAME", kOrange+4,kOrange+4, 20, 1.2  );

  // save output
  fOut    -> cd();
  h_pi0   -> Write();
  h_g     -> Write();
  h_Ratio -> Write();
  h_Diff  -> Write();
  f1      -> Write();
  if (doComparison)
    fPyth -> Write();
  c       -> Write();
  c       -> Close();
  c1      -> Write();
  c1      -> Close();
  fOut    -> Close();

}
