#include <Riostream.h>
#include <string>
#include "TROOT.h"
#include <vector>
#include <sstream>
#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"


/*--------------------------------------------------------------------*/
void makeNicePlotStyle(TH1* plot)
/*--------------------------------------------------------------------*/
{ 
  plot->GetXaxis()->CenterTitle(true);
  plot->GetYaxis()->CenterTitle(true);
  plot->GetXaxis()->SetTitleFont(42); 
  plot->GetYaxis()->SetTitleFont(42);  
  plot->GetXaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetXaxis()->SetTitleOffset(0.9);
  plot->GetYaxis()->SetTitleOffset(1.6);
  plot->GetXaxis()->SetLabelFont(42);
  plot->GetYaxis()->SetLabelFont(42);
  plot->GetYaxis()->SetLabelSize(.05);
  plot->GetXaxis()->SetLabelSize(.05);
}

//*************************************************************
double landaufun(double *x,double *par){
//*************************************************************
   double t;
   t=par[0]*TMath::Landau(x[0],par[1],par[2],0); //p1 is MP,p2 is sigma
   return t;
}

//*************************************************************
Double_t langaufun(Double_t *x, Double_t *par)  
//*************************************************************
{

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;     // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  return (par[2] * step * sum * invsq2pi / par[3]);
}


//*************************************************************
TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
//*************************************************************
{
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function              
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[4]  reasonable start values for the fit
  //   parlimitslo[4]  lower parameter limits
  //   parlimitshi[4]  upper parameter limits
  //   fitparams[4]    returns the final fit parameters
  //   fiterrors[4]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf
  
  Int_t i;
  Char_t FunName[100];
  
  std::cout<<"within langaufit"<<std::endl;

  sprintf(FunName,"Fitfcn_%s",his->GetName());
  
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;
  
  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");
  
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  
  his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
  
  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  
  return (ffit);              // return fit function
  
}

//*************************************************************
void AnalyzeChargeDistibutions(TString file) 
//*************************************************************
{

  TFile *f = TFile::Open(file);
  if (f->IsZombie()) {
    //something very wrong, cannot use this file, exit
    std::cout<< "puppa" << std::endl;
    return;
  }

  TString dets[10] = {"TIB_L1","TIB_L2","TIB_L3","TIB_L4","TOB_L1","TOB_L2","TOB_L3","TOB_L4","TOB_L5","TOB_L6"};
  TH1F *histosToFit[10];

  TCanvas* dummyC = new TCanvas("Canv","Canv",800,600);
  dummyC->cd()->SetLeftMargin(0.17);
  dummyC->cd()->SetRightMargin(0.07);

  dummyC->Print("fits.pdf[");

  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  //fr[0]=0.01*hSNR->GetMean();
  //fr[1]=10.0*hSNR->GetMean();
  
  fr[0] = 200.;
  fr[1] = 500.;
  
  // limits LOW
  pllo[0]=2.;          // Width (scale) parameter of Landau density			 
  pllo[1]=210; 	       // Most Probable (MP, location) parameter of Landau density	 
  pllo[2]=100.0;      // Total area (integral -inf to inf, normalization constant)	 
  pllo[3]=2.;	       // Width (sigma) of convoluted Gaussian function              
  
  // limits HIGH       
  plhi[0]=1000.0;      // Width (scale) parameter of Landau density		      	  
  plhi[1]=500.0;       // Most Probable (MP, location) parameter of Landau density   
  plhi[2]=99999999.0;  // Total area (integral -inf to inf, normalization constant)  
  plhi[3]=1000.0;      // Width (sigma) of convoluted Gaussian function             
  
  // starting values
  sv[0]=10.0;         // Width (scale) parameter of Landau density		  
  sv[1]=300.0; 	       // Most Probable (MP, location) parameter of Landau density	  
  sv[2]=25000.0;       // Total area (integral -inf to inf, normalization constant)  
  sv[3]=10.0;	       // Width (sigma) of convoluted Gaussian function              
      
  for(Int_t i=0;i<10;i++){
    std::cout<<"Analyzing: "<<dets[i].Data()<<std::endl;
    histosToFit[i] = (TH1F*)f->Get(Form("SiStripGainsValidator/ClusterCharge/clusterChargeOverPathNewG2_%s",dets[i].Data()));
    makeNicePlotStyle(histosToFit[i]);

    // Fitting SNR histo
    printf("Fitting...\n");
    
    // make sure there is at least 1000k tracks in the plot
    if(histosToFit[i]->GetEntries()<5000){
      std::cout<<"not enough clusters to make sensible plot"<<std::endl;
      return;
    }

    Double_t chisqr;
    Int_t    ndf;

    TF1 *fitsnr = langaufit(histosToFit[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  
    //Double_t SNRPeak, SNRFWHM;
    //langaupro(fp,SNRPeak,SNRFWHM);
    
    printf("Fitting done\nPlotting results...\n");
    
    fitsnr->SetLineColor(kBlue);
    fitsnr->SetLineWidth(2);
    
    // Global style settings
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    gStyle->SetLabelSize(0.03,"x");
    gStyle->SetLabelSize(0.03,"y");
    
    histosToFit[i]->GetXaxis()->SetRange(0.,1000.);
    histosToFit[i]->SetMarkerStyle(20);
    histosToFit[i]->SetMarkerSize(1);
 
    histosToFit[i]->Draw();
    fitsnr->Draw("lsame");

    dummyC->Update();
    
    TLine *l = new TLine(300.,dummyC->GetUymin(),300.,dummyC->GetUymax());
    l->SetLineColor(kCyan);
    l->SetLineWidth(2);
    l->SetLineStyle(8);
    l->Draw("same");

    dummyC->Print("fits.pdf");
    delete fitsnr;
    delete l;
  }
  /*
  TPaveText *pt = new TPaveText(0.4,0.75,0.6,0.87,"NDC");
  pt->SetFillColor(10);
  pt->SetTextColor(1);
  pt->SetTextFont(42);
  pt->SetTextAlign(11);
  pt->Draw("same");
  */

  dummyC->Print("fits.pdf]");

}
