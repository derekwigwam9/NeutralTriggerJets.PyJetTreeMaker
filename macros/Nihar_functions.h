
///________
TH1 Draw_hist(TH1 *gr, Option_t *option, Double_t mColor, Double_t mlineColor, Double_t mStyle, Double_t mSize )
{
  gr->SetFillColor(mColor);
  gr->SetFillStyle(mStyle);
  gr->SetMarkerColor(mColor);
  gr->SetMarkerStyle(mStyle);
  gr->SetMarkerSize(mSize);
  gr->SetLineColor(mlineColor);
  //  gr->Draw("HIST");

  //  return gr->Draw("HIST");
  return gr->Draw(option);

}

TH1 Draw_hist(TH1 *gr, Option_t *option, Double_t mColor, Double_t mlineColor, Double_t mStyle, Double_t mSize )
{
  gr->SetFillColor(mColor);
  gr->SetFillStyle(mStyle);
  gr->SetMarkerColor(mColor);
  gr->SetMarkerStyle(mStyle);
  gr->SetMarkerSize(mSize);
  gr->SetLineColor(mlineColor);
  //  gr->Draw("HIST");

  //  return gr->Draw("HIST");
  return gr->Draw(option);

}

TH1 Draw_hist(TGraphErrors *gr, Option_t *option, Double_t mColor, Double_t mlineColor, Double_t mStyle, Double_t mSize )
{
  gr->SetFillColor(mColor);
  gr->SetFillStyle(mStyle);
  gr->SetMarkerColor(mColor);
  gr->SetMarkerStyle(mStyle);
  gr->SetMarkerSize(mSize);
  gr->SetLineColor(mlineColor);
  //  gr->Draw("HIST");

  //  return gr->Draw("HIST");
  return gr->Draw(option);

   

  
}

//______

TH1 Draw_hist(TH1*h, Option_t *option, Double_t mColor, Double_t mlineColor, Double_t mStyle, Double_t mSize )
{
  h->SetFillColor(mColor);
  h->SetFillStyle(mStyle);
  h->SetMarkerColor(mColor);
  h->SetMarkerStyle(mStyle);
  h->SetMarkerSize(mSize);
  h->SetLineColor(mlineColor);
  //  h->Draw("HIST");

  //  return h->Draw("HIST");
  return h->Draw(option);

}

///_______
 
//__________________________
TLegend Draw_Legend(TLegend *leg, double x1, double y1, double x2, double y2, const char *t1, const char *t2, const char *t3, const char *t4,const char *t5)
{
  leg = new TLegend(x1,y1,x2,y2);
  leg->AddEntry((TObject*)0, t1, "");
  leg->AddEntry((TObject*)0, t2, "");
  leg->AddEntry((TObject*)0, t3, "");
  leg->AddEntry((TObject*)0, t4, "");
  leg->AddEntry((TObject*)0, t5, "");
  
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);  
  leg->SetFillColor(0);  
  leg->SetLineColor(0);
  //  leg->AddEntry(h_Datajet_JetPtCorr_gamma,"#pi^{0}-h^{#pm}","p");  
  //  leg->AddEntry(h_normalized_full,"#gamma_{dir}-h^{#pm}","p");
  leg->Draw();

}

TLegend Draw_Legend_Marker(TLegend *leg, double x1, double y1, double x2, double y2, TH1* h1, TH1* h2, const char *t1, const char *t2, Option_t *op1, Option_t *op2)
{
  leg = new TLegend(x1,y1,x2,y2);
  //  leg->AddEntry((TObject*)0, t1, "");
  //  leg->AddEntry((TObject*)0, t2, "");
  //  leg->AddEntry((TObject*)0, t3, "");
  //  leg->AddEntry((TObject*)0, t4, "");
  
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);  
  leg->SetFillColor(0);  
  leg->SetLineColor(0);
  leg->AddEntry(h1,t1, op1);  
  leg->AddEntry(h2,t2,op2);
  leg->Draw();

}


//____________________
TH1D* DirectPhotonCorrectedHist(TH1* hg, TH1* hpi, Double_t R,Double_t R_err,Int_t Nbin, Int_t Nbin_min, Int_t Nbin_max)
{
  // This method used in gamma-hadron analysis
  // for detail you can look at analysis note http://www.star.bnl.gov/protected/jetcorr/nihar/PaperProposal/AnalysisNote_v5.pdf
  // Equation no. 5 in the analysis note 
  
  
  Int_t Nbin_hg = hg->GetNbinsX();
  Int_t Nbin_hpi = hpi->GetNbinsX();

  if(Nbin_hg != Nbin_hpi)
    {
      cout<<"CHECK: # of bin not same ( "<<" h_g: "<< Nbin_hg <<" , h_pi: "<< Nbin_hpi<< ")" <<endl;
      return kFALSE;;
    }


  TH1D *hCorrDirPho = new TH1D("hCorrDirPho","",Nbin_hg,Nbin_min,Nbin_max);hCorrDirPho ->Sumw2();
  Int_t Sub_pnt =0;
  for(int i=1; i <= Nbin_hg; i++)
    {

      double xg =  hg->GetBinCenter(i);
      double xpi =  hpi->GetBinCenter(i);
      
      double yg =  hg->GetBinContent(i);
      double yg_err =  hg->GetBinError(i);
      double ypi =  hpi->GetBinContent(i);
      double ypi_err =  hpi->GetBinError(i);

      //      double sub_y = yg - ypi;
      double dirPho_yield_bin = (yg -  R*ypi)/(1-R);
      double dirPho_yield_bin_err = sqrt(
					       pow(1./(1.-R),2)*( 
				      pow(yg_err,2)+
				      pow(R,2)*pow(ypi_err,2)+
				      pow((yg-ypi)/(1-R),2)*pow(R_err,2)
								  )
					 );

      /*      cout<<"yg: "<<yg<<" yg_err: "<<yg_err
	<<"ypi: "<<ypi<<" ypi_err: "<<ypi_err
	  <<" dirPho y: "<<dirPho_yield_bin<<" dirPho_yield err: "<<dirPho_yield_bin_err<<endl;
      */
      
      hCorrDirPho->SetBinContent(i,dirPho_yield_bin);
      hCorrDirPho->SetBinError(i,dirPho_yield_bin_err);

      Sub_pnt++;
    }
    hCorrDirPho->SetEntries(Sub_pnt);


    //    hSub->Draw(option);
    return hCorrDirPho;
}



//____________________
//void Subtract_2Hist(TH1* h1, TH1* h2, Int_t Nbin, Int_t Nbin_min, Int_t Nbin_max,Option_t *option )
//void Subtract_2Hist(TH1* h1, TH1* h2, Int_t Nbin, Int_t Nbin_min, Int_t Nbin_max,TH1D &hSub )
TH1D* Subtract_2Hist(TH1* h1, TH1* h2, Int_t Nbin, Int_t Nbin_min, Int_t Nbin_max )
{

  // h1: RE   and h2: ME
  
  Int_t Nbin_h1 = h1->GetNbinsX();
  Int_t Nbin_h2 = h2->GetNbinsX();

  if(Nbin_h1 != Nbin_h2)
    {
      cout<<"CHECK: # of bin not same ( "<<" h_RE: "<< Nbin_h1<<" , h_ME: "<< Nbin_h2<< ")" <<endl;
      return kFALSE;;
    }

  //  TH1D *hSub = new TH1D("hSub","",Nbin,Nbin_min,Nbin_max);hSub->Sumw2();
  TH1D *hSub1 = new TH1D("hSub1","",Nbin,Nbin_min,Nbin_max);hSub->Sumw2();
  Int_t Sub_pnt =0;
  for(int i=1; i <= Nbin_h1; i++)
    {

      double x1 =  h1->GetBinCenter(i);
      double x2 =  h2->GetBinCenter(i);
      
      double y1 =  h1->GetBinContent(i);
      double y1_err =  h1->GetBinError(i);
      double y2 =  h2->GetBinContent(i);
      double y2_err =  h2->GetBinError(i);

      double sub_y = y1 - y2;

      //      cout<<i<<" bin cntr=  "<<x1<<"  "<<x2<<"  "<<y1<<" - "<<y2<<" = "<<sub_y<<endl;//"  "<<y2<<"+-"<<y2_err<<"  "<<y1<<"+-"<<y1_err<<endl;
      
      hSub1->SetBinContent(i,sub_y);
      hSub1->SetBinError(i,y1_err);

      Sub_pnt++;
    }
    hSub1->SetEntries(Sub_pnt);


    //    hSub->Draw(option);
    return hSub1;
}


//________
Double_t LevyFitFunc_pT(Double_t* x_val, Double_t* par)
{
    // One over pT term is removed -> original pT distribution
    // Fit function for d2N/(2pi dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = pT*B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------


/*
//----------------------------------------------------------------------------------------
Int_t Hist_interpolate_and_error(TH1D* hist, Double_t x, Double_t &Int_val, Double_t &Int_err)
{
    // Linear interpolation of a histogram
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[2]   = {0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetEntries() == 0) // no data poins to interpolate
    {
        return_val = -1;
        Int_val = 0;
        Int_err = 0;
        //cout << "No entries in histogram" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t binx = 1; binx < hist->GetNbinsX(); binx++)
        {
            Double_t bin_error_val   = hist->GetBinError(binx);
            Double_t bin_x_pos       = hist->GetBinCenter(binx);
            if(bin_error_val != 0)
            {
                err_counter++;
                bin_entries[1] = hist->GetBinContent(binx);
                bin_error[1]   = hist->GetBinError(binx);
                bin_x_val[1]   = hist->GetBinCenter(binx);
                if(bin_x_pos >= x)
                {
                    binx_high = binx;
                    flag_max  = 1;
                    break;
                }
                else flag_max = 0;
            }
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val = bin_entries[1];
            Int_err = bin_error[1];
            return return_val;
        }
        for(Int_t binx_low = binx_high; binx_low > 0; binx_low--)
        {
            bin_entries[0] = hist->GetBinContent(binx_low);
            bin_error[0]   = hist->GetBinError(binx_low);
            bin_x_val[0]   = hist->GetBinCenter(binx_low);
            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        if(bin_error[0] != 0 && bin_error[1] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;
                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[1];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;
                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val = bin_entries[0];
                Int_err = bin_error[0];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------

*/
