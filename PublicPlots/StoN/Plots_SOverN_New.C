//#include tdrStyle.C
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TGAxis.h"
#include "TROOT.h"
#include "CMS_lumi.C"

Double_t langaufun(Double_t *x, Double_t *par) {
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
  Double_t invsq2pi = 0.3989422804014;  // (2 pi)^(-1/2)
  Double_t mpshift = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;  // number of convolution steps
  Double_t sc = 5.0;    // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp - xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for (i = 1.0; i <= np / 2; i++) {
    xx = xlow + (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);

    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

void settdrStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle", "Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600);  //Height of canvas
  tdrStyle->SetCanvasDefW(600);  //Width of canvas
  tdrStyle->SetCanvasDefX(0);    //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  //tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  //tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0);  // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.06);
  tdrStyle->SetPadBottomMargin(0.18);
  tdrStyle->SetPadLeftMargin(0.20);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.06, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20., 20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}

TH1 *GetSOverNHisto(TDirectory *current_sourcedir) {
  TIter nextkey(current_sourcedir->GetListOfKeys());
  TKey *key, *oldkey = 0;
  TH1 *h = 0;
  while ((key = (TKey *)nextkey())) {
    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(), key->GetName()))
      continue;

    current_sourcedir->cd();
    TObject *obj = key->ReadObj();
    TString obj_name(obj->GetName());
    //cout<<obj->GetName()<<endl;

    if (obj->IsA()->InheritsFrom(TH1::Class())) {
      // descendant of TH1 -> merge it
      if (!obj_name.Contains("Summary_ClusterStoNCorr_OnTrack__") &&
          !obj_name.Contains("Summary_ClusterStoNCorr__OnTrack__"))
        continue;
      cout << "Getting histogram " << obj_name << endl;
      h = (TH1 *)obj;
    }
  }

  //return (TH1*)h->Clone("new");
  return h;
}

void DrawSOverNHisto(TCanvas *&c, TH1 *h, int run, string subdet, bool thin = true) {
  c->cd();
  TGaxis::SetMaxDigits(4);  // Pour TEC
  h->SetTitle(0);
  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);
  h->SetMarkerSize(0.75);
  h->SetMarkerStyle(20);
  h->SetLineWidth(2);
  h->GetXaxis()->SetRangeUser(0, 100);
  h->GetYaxis()->SetTitle("SiStrip Clusters");
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.3);
  h->GetXaxis()->SetTitle("Signal-To-Noise ratio");
  h->GetXaxis()->SetTitleOffset(1.);
  h->Draw("E1");

  //TF1 *fit = new TF1("mylandau", "langaus", 15.5, 100);
  TF1 *fit = new TF1("mylandau", langaufun, 15.5, 100, 4);
  fit->SetParameter(0, 2.);
  fit->SetParameter(1, 18);
  fit->SetParameter(2, h->Integral());
  fit->SetParameter(3, 0.5);
  if (subdet == "TOB")
    fit->SetRange(17.5, 40);
  if (subdet == "TEC")
    fit->SetRange(12.5, 40);
  if (subdet == "TIB")
    fit->SetRange(10.5, 30);
  if (subdet == "TID")
    fit->SetRange(10.5, 30);
  //else fit->SetRange(12.5,40);
  h->Fit("mylandau", "R");
  double mpv = fit->GetParameter(1);
  double mpv_err = fit->GetParError(1);
  fit->SetLineWidth(2.);
  //fit->SetLineStyle(9);
  fit->SetLineColor(kRed);
  fit->Draw("same");

  TPaveText *text = new TPaveText(0.45, 0.55, 0.70, 0.75, "NDC");
  //text->AddText("CMS Preliminary");
  if (run == 251883)
    text->AddText("2015 Data, 50ns");
  if (run == 258749)
    text->AddText("2015 Data, 25ns");
  if (run == 260627)
    text->AddText("2015 Data, 25ns");
  if (run == 262235)
    text->AddText("2015 Data, sqrt(s)=5TeV");
  if (run == 273450)
    text->AddText("2016 Data");
  if (run == 278770)
    text->AddText("2016 Data - old APV setting");
  if (run == 278803)
    text->AddText("2016 Data - new APV setting");
  if (run == 283407)
    text->AddText("2016 Data");
  //text->AddText("\n");
  text->AddText("#LTinst. lumi.#GT: 1.5#times10^{34} cm^{-2}s^{-1}");
  //text->AddText(Form("%s MPV: %.2f #pm %.3f", subdet.c_str(), mpv, mpv_err));
  if (subdet == "TEC") {
    if (thin)
      text->AddText(Form("%s thin MPV: %.1f", subdet.c_str(), mpv));
    else
      text->AddText(Form("%s thick MPV: %.1f", subdet.c_str(), mpv));
  } else
    text->AddText(Form("%s MPV: %.1f", subdet.c_str(), mpv));
  text->SetTextFont(42);
  text->SetTextAlign(12);
  text->SetTextSize(0.04);
  //text->SetMargin(0.01);
  //text->SetY1NDC(0.5);
  //text->SetY2NDC(0.75);
  text->SetBorderSize(0);
  text->SetShadowColor(10);
  text->SetFillColor(10);
  text->Draw();
}

void Plot(TCanvas *c, TFile *f, int run, string subdet = "TIB", bool thin = true) {
  TH1 *h = 0;

  writeExtraText = true;       // if extra text
  extraText = "Preliminary";   // default extra text is "Preliminary"
  lumi_8TeV = "19.1 fb^{-1}";  // default is "19.7 fb^{-1}"
  lumi_7TeV = "4.9 fb^{-1}";   // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  // For TEC differentiate sensors thicknesses and merge sides
  string sides[2];
  sides[0] = "MINUS";
  sides[1] = "PLUS";

  TH1 *hTECthin = 0;
  TH1 *hTECthick = 0;

  if (subdet == "TEC") {
    // Merge TEC+ and TEC-
    for (unsigned iside = 0; iside < 2; iside++) {
      //TEC thin sensors
      for (unsigned iring = 1; iring < 8; iring++) {
        f->cd(
            Form("DQMData/Run %i/SiStrip/Run summary/MechanicalView/TEC/%s/ring_%i", run, sides[iside].c_str(), iring));

        TDirectory *current_sourcedir = gDirectory;

        TH1 *h = GetSOverNHisto(current_sourcedir);
        if (!h)
          continue;

        if (iring < 5) {
          if (!hTECthin)
            hTECthin = h;
          else
            hTECthin->Add(h);
        } else {
          if (!hTECthick)
            hTECthick = h;
          else
            hTECthick->Add(h);
        }
      }
    }
  }
  // For the other subdets (TIB, TOB, TID)
  else {
    f->cd(Form("DQMData/Run %i/SiStrip/Run summary/MechanicalView/%s", run, subdet.c_str()));
    if (subdet == "TID")
      gDirectory->cd("MINUS");

    TDirectory *current_sourcedir = gDirectory;
    h = GetSOverNHisto(current_sourcedir);

    // Merge both sides for TID
    if (subdet == "TID") {
      f->cd(Form("DQMData/Run %i/SiStrip/Run summary/MechanicalView/%s/PLUS", run, subdet.c_str()));
      current_sourcedir = gDirectory;
      TH1 *h2 = GetSOverNHisto(current_sourcedir);

      if (h2)
        h->Add(h2);
    }
  }

  if (!h && !hTECthin && !hTECthick)
    return;
  if (thin && !h && hTECthin)
    h = hTECthin;
  if (!thin && !h && hTECthick)
    h = hTECthick;

  // Draw histos
  if (h)
    DrawSOverNHisto(c, h, run, subdet, thin);
  c->Modified();
  c->Update();
}

float *PlotPerLayer(TCanvas *c, TFile *f, int run, string subdet = "TIB") {
  int nlayer = 4;
  if (subdet == "TIB")
    nlayer = 4;
  if (subdet == "TOB")
    nlayer = 6;
  float *MPVs = new float[10];

  for (int ilayer = 1; ilayer < nlayer + 1; ilayer++) {
    f->cd(Form("DQMData/Run %i/SiStrip/Run summary/MechanicalView/%s/layer_%i", run, subdet.c_str(), ilayer));
    TDirectory *current_sourcedir = gDirectory;
    TH1 *h = GetSOverNHisto(current_sourcedir);
    if (!h)
      continue;

    c->cd();
    h->SetTitle(0);
    h->GetXaxis()->SetRangeUser(0, 100);
    h->GetYaxis()->SetTitle("Entries");
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetTitle("Signal-To-Noise");
    h->GetXaxis()->SetTitleOffset(1.);
    if (ilayer == 1)
      h->Draw("E1");
    else
      h->Draw("same");

    //TF1 *fit = new TF1("mylandau", "langaus", 15.5, 100);
    TF1 *fit = new TF1("mylandau", langaufun, 15.5, 100, 4);
    fit->SetParameter(0, 2.);
    fit->SetParameter(1, 18);
    fit->SetParameter(2, h->Integral());
    fit->SetParameter(3, 0.5);
    //if(subdet=="TOB") fit->SetRange(10.5,100);
    if (subdet == "TOB")
      fit->SetRange(15.5, 100);
    else
      fit->SetRange(15.5, 100.);
    h->Fit("mylandau", "R");
    double mpv = fit->GetParameter(1);
    double mpv_err = fit->GetParError(1);
    fit->Draw("same");

    TPaveText *text = new TPaveText(0.5, 0.63, 0.82, 0.8, "NDC");
    text->AddText("CMS Preliminary");
    //if(run==251883) text->AddText("2015 Data, 50ns");
    //if(run==258749) text->AddText("2015 Data, 25ns");
    //if(run==260627) text->AddText("2015 Data, 25ns");
    //text->AddText(Form("%s MPV: %.2f #pm %.3f", subdet.c_str(), mpv, mpv_err));
    text->AddText(Form("%s MPV: %.1f", subdet.c_str(), mpv));
    text->SetTextAlign(12);
    //text->SetMargin(0.01);
    //text->SetY1NDC(0.5);
    //text->SetY2NDC(0.75);
    text->SetBorderSize(0);
    text->SetShadowColor(10);
    text->SetFillColor(10);
    text->Draw();
    c->Modified();
    c->Update();

    cout << "layer: " << ilayer << "  mpv: " << mpv << endl;
    MPVs[ilayer - 1] = mpv;
  }
  return MPVs;
}

void Plots_SOverN_New(string file = "DQM_V0001_R000278803__ZeroBias__Run2016F-PromptReco-v1__DQMIO.root") {
  // R000251883__ZeroBias__Run2015B-PromptReco-v1
  // R000251883__Global__CMSSW_X_Y_Z__RECO
  // R000258749__ZeroBias__Run2015D-PromptReco-v4
  // R000258749__Global__CMSSW_X_Y_Z__RECO
  // R000257682__ZeroBias__Run2015D-PromptReco-v3 // PIX BAD
  // R000260627__ZeroBias__Run2015D-PromptReco-v4
  // R000262235__ZeroBias__Run2015E-PromptReco-v1
  // R000273450__ZeroBias__Run2016B-PromptReco-v2
  // R000278770__ZeroBias__Run2016F-PromptReco-v1
  // R000278803__ZeroBias__Run2016F-PromptReco-v1
  //gROOT->LoadMacro("tdrStyle.C");
  settdrStyle();

  // Could be get from the file name
  int run = 283407;          // for choosing plots legends
  string suffix = "_2016H";  // for plots file name

  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas *c = new TCanvas();
  TFile *f = new TFile(file.c_str(), "read");
  if (!f)
    return;
  Plot(c, f, run, "TIB");
  c->Modified();
  c->Update();
  getchar();
  CMS_lumi(c, 0, 33);
  c->Print(Form("TIB_SOverN_run%i%s.png", run, suffix.c_str()));
  c->Print(Form("TIB_SOverN_run%i%s.pdf", run, suffix.c_str()));

  Plot(c, f, run, "TOB");
  c->Modified();
  c->Update();
  getchar();
  CMS_lumi(c, 0, 33);
  c->Print(Form("TOB_SOverN_run%i%s.png", run, suffix.c_str()));
  c->Print(Form("TOB_SOverN_run%i%s.pdf", run, suffix.c_str()));

  Plot(c, f, run, "TID");
  c->Modified();
  c->Update();
  getchar();
  CMS_lumi(c, 0, 33);
  c->Print(Form("TID_SOverN_run%i%s.png", run, suffix.c_str()));
  c->Print(Form("TID_SOverN_run%i%s.pdf", run, suffix.c_str()));

  Plot(c, f, run, "TEC", true);
  c->Modified();
  c->Update();
  getchar();
  CMS_lumi(c, 0, 33);
  c->Print(Form("TECthin_SOverN_run%i%s.png", run, suffix.c_str()));
  c->Print(Form("TECthin_SOverN_run%i%s.pdf", run, suffix.c_str()));

  Plot(c, f, run, "TEC", false);
  c->Modified();
  c->Update();
  getchar();
  CMS_lumi(c, 0, 33);
  c->Print(Form("TECthick_SOverN_run%i%s.png", run, suffix.c_str()));
  c->Print(Form("TECthick_SOverN_run%i%s.pdf", run, suffix.c_str()));

  // Plots per layer for the barrel

  float *MPV_TIB;
  MPV_TIB = PlotPerLayer(c, f, run, "TIB");
  c->Modified();
  c->Update();
  getchar();

  float *MPV_TOB;
  MPV_TOB = PlotPerLayer(c, f, run, "TOB");
  c->Modified();
  c->Update();
  getchar();

  c->cd();
  TGraph *glayer = new TGraph();
  for (int ilayer = 0; ilayer < 4; ilayer++)
    glayer->SetPoint(ilayer, ilayer + 1, MPV_TIB[ilayer]);
  for (int ilayer = 0; ilayer < 6; ilayer++)
    glayer->SetPoint(ilayer + 4, ilayer + 5, MPV_TOB[ilayer]);

  TH1 *h = glayer->GetHistogram();
  h->SetTitle(0);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitle("S/N");
  TAxis *ax = h->GetXaxis();
  //ax->SetBinLabel(5, "TIB");
  int ibin = ax->FindBin(1);
  ax->SetBinLabel(ibin, "TIB1");
  ibin = ax->FindBin(2);
  ax->SetBinLabel(ibin, "TIB2");
  ibin = ax->FindBin(3);
  ax->SetBinLabel(ibin, "TIB3");
  ibin = ax->FindBin(4);
  ax->SetBinLabel(ibin, "TIB4");
  //ax->SetBinLabel(42, "TOB");
  ibin = ax->FindBin(5);
  ax->SetBinLabel(ibin, "TOB1");
  ibin = ax->FindBin(6);
  ax->SetBinLabel(ibin, "TOB2");
  ibin = ax->FindBin(7);
  ax->SetBinLabel(ibin, "TOB3");
  ibin = ax->FindBin(8);
  ax->SetBinLabel(ibin, "TOB4");
  ibin = ax->FindBin(9);
  ax->SetBinLabel(ibin, "TOB5");
  ibin = ax->FindBin(10);
  ax->SetBinLabel(ibin, "TOB6");
  ax->LabelsOption("u");
  glayer->Draw("AP");
  c->Modified();
  c->Update();
  getchar();

  c->Print(Form("SoverN_vslayer_run%i%s.png", run, suffix.c_str()));
  c->Print(Form("SoverN_vslayer_run%i%s.pdf", run, suffix.c_str()));
}
