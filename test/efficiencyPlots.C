#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "tdrstyle.C"


using namespace std;


void efficiency_curves_grooming(const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
				const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);

  // signal files
  TFile *file_S          = new TFile("output_files/WW500DefaultJetMass.root");
  TFile *file_S_Filtered = new TFile("output_files/WW500FilteredJetMass.root");
  TFile *file_S_Pruned   = new TFile("output_files/WW500PrunedJetMass.root");
  TFile *file_S_Trimmed  = new TFile("output_files/WW500TrimmedJetMass.root");

  // background files
  TFile *file_B          = new TFile("output_files/QCDPythiaDefaultJetMass.root");
  TFile *file_B_Filtered = new TFile("output_files/QCDPythiaFilteredJetMass.root");
  TFile *file_B_Pruned   = new TFile("output_files/QCDPythiaPrunedJetMass.root");
  TFile *file_B_Trimmed  = new TFile("output_files/QCDPythiaTrimmedJetMass.root");

  // signal histograms
  TH2D *h2_nPV_JetMass_S          = (TH2D*)file_S->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Filtered = (TH2D*)file_S_Filtered->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Pruned   = (TH2D*)file_S_Pruned->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Trimmed  = (TH2D*)file_S_Trimmed->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S          = (TH2D*)file_S->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Filtered = (TH2D*)file_S_Filtered->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Pruned   = (TH2D*)file_S_Pruned->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Trimmed  = (TH2D*)file_S_Trimmed->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B          = (TH2D*)file_B->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Filtered = (TH2D*)file_B_Filtered->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Pruned   = (TH2D*)file_B_Pruned->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Trimmed  = (TH2D*)file_B_Trimmed->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B          = (TH2D*)file_B->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Filtered = (TH2D*)file_B_Filtered->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Pruned   = (TH2D*)file_B_Pruned->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Trimmed  = (TH2D*)file_B_Trimmed->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S          = h2_nPV_JetMass_S->Integral(fPVLow,fPVHigh,0,201);
  double denom_S_Filtered = h2_nPV_JetMass_S_Filtered->Integral(fPVLow,fPVHigh,0,201);
  double denom_S_Pruned   = h2_nPV_JetMass_S_Pruned->Integral(fPVLow,fPVHigh,0,201);
  double denom_S_Trimmed  = h2_nPV_JetMass_S_Trimmed->Integral(fPVLow,fPVHigh,0,201);

  // background denominator counts
  double denom_B = h2_nPV_JetMass_B->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_Filtered = h2_nPV_JetMass_B_Filtered->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_Pruned   = h2_nPV_JetMass_B_Pruned->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_Trimmed  = h2_nPV_JetMass_B_Trimmed->Integral(fPVLow,fPVHigh,0,201);

  // Default jets
  TGraph *g_eff = new TGraph(21);
  g_eff->SetName("g_eff");
  g_eff->SetLineColor(1);
  g_eff->SetLineWidth(2);
  g_eff->SetLineStyle(1);
  g_eff->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff->SetPoint(i,(num_B/denom_B),(num_S/denom_S));
  }

  // Filtered jets
  TGraph *g_eff_Filtered = new TGraph(21);
  g_eff_Filtered->SetName("g_eff_Filtered");
  g_eff_Filtered->SetLineColor(kBlue);
  g_eff_Filtered->SetLineWidth(2);
  g_eff_Filtered->SetLineStyle(7);
  g_eff_Filtered->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Filtered->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Filtered->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Filtered->SetPoint(i,(num_B/denom_B_Filtered),(num_S/denom_S_Filtered));
  }

  // Pruned jets
  TGraph *g_eff_Pruned = new TGraph(21);
  g_eff_Pruned->SetName("g_eff_Pruned");
  g_eff_Pruned->SetLineColor(kRed);
  g_eff_Pruned->SetLineWidth(2);
  g_eff_Pruned->SetLineStyle(3);
  g_eff_Pruned->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Pruned->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Pruned->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Pruned->SetPoint(i,(num_B/denom_B_Pruned),(num_S/denom_S_Pruned));
  }

  // Trimmed jets
  TGraph *g_eff_Trimmed = new TGraph(21);
  g_eff_Trimmed->SetName("g_eff_Trimmed");
  g_eff_Trimmed->SetLineColor(kGreen+2);
  g_eff_Trimmed->SetLineWidth(2);
  g_eff_Trimmed->SetLineStyle(8);
  g_eff_Trimmed->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Trimmed->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Trimmed->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Trimmed->SetPoint(i,(num_B/denom_B_Trimmed),(num_S/denom_S_Trimmed));
  }
  
  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff->Draw("L");
  g_eff_Filtered->Draw("L");
  g_eff_Pruned->Draw("L");
  g_eff_Trimmed->Draw("L");

  TLegend *legend = new TLegend(.55,.25,.85,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_Trimmed, "Trimmed jet mass","l");
  legend->AddEntry(g_eff_Filtered, "Filtered jet mass","l");
  legend->AddEntry(g_eff_Pruned, "Pruned jet mass","l");
  legend->AddEntry(g_eff, "Default jet mass","l");
  legend->Draw();
  
  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.14,0.97, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());
}


void efficiency_curves_jetpt(const int fPVLow, const int fPVHigh, const string& fTitle,
		             const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);

  // signal files
  TFile *file_S      = new TFile("output_files/WW500TrimmedJetMass.root");
  TFile *file_S_incl = new TFile("output_files/WWTrimmedJetMass.root");

  // background files
  TFile *file_B  = new TFile("output_files/QCDPythiaTrimmedJetMass.root");

  // signal histograms
  TH2D *h2_nPV_JetMass_S_Pt300to500 = (TH2D*)file_S_incl->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt300to500");
  TH2D *h2_nPV_JetMass_S_Pt500to700 = (TH2D*)file_S->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt500to700");
  TH2D *h2_nPV_JetMass_S_Pt700to900 = (TH2D*)file_S->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt700to900");
  TH2D *h2_nPV_JetMass_S_Pt900toInf = (TH2D*)file_S->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt900toInf");

  TH2D *h2_nPV_tau2tau1_S_Pt300to500 = (TH2D*)file_S_incl->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt300to500");
  TH2D *h2_nPV_tau2tau1_S_Pt500to700 = (TH2D*)file_S->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt500to700");
  TH2D *h2_nPV_tau2tau1_S_Pt700to900 = (TH2D*)file_S->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt700to900");
  TH2D *h2_nPV_tau2tau1_S_Pt900toInf = (TH2D*)file_S->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt900toInf");

  // background histograms
  TH2D *h2_nPV_JetMass_B_Pt300to500 = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt300to500");
  TH2D *h2_nPV_JetMass_B_Pt500to700 = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt500to700");
  TH2D *h2_nPV_JetMass_B_Pt700to900 = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt700to900");
  TH2D *h2_nPV_JetMass_B_Pt900toInf = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_JetMass_Pt900toInf");

  TH2D *h2_nPV_tau2tau1_B_Pt300to500 = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt300to500");
  TH2D *h2_nPV_tau2tau1_B_Pt500to700 = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt500to700");
  TH2D *h2_nPV_tau2tau1_B_Pt700to900 = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt700to900");
  TH2D *h2_nPV_tau2tau1_B_Pt900toInf = (TH2D*)file_B->Get("rutgersJetAnalyzer/h2_nPV_tau2tau1_Pt900toInf");

  // signal denominator counts
  double denom_S_Pt300to500 = h2_nPV_JetMass_S_Pt300to500->Integral(fPVLow,fPVHigh,0,201);
  double denom_S_Pt500to700 = h2_nPV_JetMass_S_Pt500to700->Integral(fPVLow,fPVHigh,0,201);
  double denom_S_Pt700to900 = h2_nPV_JetMass_S_Pt700to900->Integral(fPVLow,fPVHigh,0,201);
  double denom_S_Pt900toInf = h2_nPV_JetMass_S_Pt900toInf->Integral(fPVLow,fPVHigh,0,201);

  // background denominator counts
  double denom_B_Pt300to500 = h2_nPV_JetMass_B_Pt300to500->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_Pt500to700 = h2_nPV_JetMass_B_Pt500to700->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_Pt700to900 = h2_nPV_JetMass_B_Pt700to900->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_Pt900toInf = h2_nPV_JetMass_B_Pt900toInf->Integral(fPVLow,fPVHigh,0,201);

  // Pt300to500 jets
  TGraph *g_eff_Pt300to500 = new TGraph(21);
  g_eff_Pt300to500->SetName("g_eff_Pt300to500");
  g_eff_Pt300to500->SetLineColor(kGreen+2);
  g_eff_Pt300to500->SetLineWidth(2);
  g_eff_Pt300to500->SetLineStyle(1);
  g_eff_Pt300to500->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Pt300to500->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Pt300to500->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Pt300to500->SetPoint(i,(num_B/denom_B_Pt300to500),(num_S/denom_S_Pt300to500));
  }

  // Pt500to700 jets
  TGraph *g_eff_Pt500to700 = new TGraph(21);
  g_eff_Pt500to700->SetName("g_eff_Pt500to700");
  g_eff_Pt500to700->SetLineColor(kGreen+2);
  g_eff_Pt500to700->SetLineWidth(2);
  g_eff_Pt500to700->SetLineStyle(7);
  g_eff_Pt500to700->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Pt500to700->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Pt500to700->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Pt500to700->SetPoint(i,(num_B/denom_B_Pt500to700),(num_S/denom_S_Pt500to700));
  }

  // Pt700to900 jets
  TGraph *g_eff_Pt700to900 = new TGraph(21);
  g_eff_Pt700to900->SetName("g_eff_Pt700to900");
  g_eff_Pt700to900->SetLineColor(kGreen+2);
  g_eff_Pt700to900->SetLineWidth(2);
  g_eff_Pt700to900->SetLineStyle(3);
  g_eff_Pt700to900->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Pt700to900->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Pt700to900->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Pt700to900->SetPoint(i,(num_B/denom_B_Pt700to900),(num_S/denom_S_Pt700to900));
  }

  // Pt900toInf jets
  TGraph *g_eff_Pt900toInf = new TGraph(21);
  g_eff_Pt900toInf->SetName("g_eff_Pt900toInf");
  g_eff_Pt900toInf->SetLineColor(kGreen+2);
  g_eff_Pt900toInf->SetLineWidth(2);
  g_eff_Pt900toInf->SetLineStyle(8);
  g_eff_Pt900toInf->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_Pt900toInf->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_Pt900toInf->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_Pt900toInf->SetPoint(i,(num_B/denom_B_Pt900toInf),(num_S/denom_S_Pt900toInf));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_Pt300to500->Draw("L");
  g_eff_Pt500to700->Draw("L");
  g_eff_Pt700to900->Draw("L");
  g_eff_Pt900toInf->Draw("L");

  TLegend *legend = new TLegend(.55,.25,.85,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_Pt300to500, "300<p_{T}<500 GeV","l");
  legend->AddEntry(g_eff_Pt500to700, "500<p_{T}<700 GeV","l");
  legend->AddEntry(g_eff_Pt700to900, "700<p_{T}<900 GeV","l");
  legend->AddEntry(g_eff_Pt900toInf, "p_{T}>900 GeV","l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.14,0.97, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());
}


void efficiency_curves_pileup(const string& fPtRange, const string& fTitle,
		              const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);

  // signal files
  TFile *file_S      = new TFile("output_files/WW500TrimmedJetMass.root");

  // background files
  TFile *file_B  = new TFile("output_files/QCDPythiaTrimmedJetMass.root");

  // signal histograms
  TH2D *h2_nPV_JetMass_S = (TH2D*)file_S->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S = (TH2D*)file_S->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B = (TH2D*)file_B->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B = (TH2D*)file_B->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S_nPV0to10 = h2_nPV_JetMass_S->Integral(1,11,0,201);
  double denom_S_nPV10to20 = h2_nPV_JetMass_S->Integral(11,21,0,201);
  double denom_S_nPV20to30 = h2_nPV_JetMass_S->Integral(21,31,0,201);
  double denom_S_nPV30toInf = h2_nPV_JetMass_S->Integral(31,52,0,201);

  // background denominator counts
  double denom_B_nPV0to10 = h2_nPV_JetMass_B->Integral(1,11,0,201);
  double denom_B_nPV10to20 = h2_nPV_JetMass_B->Integral(11,21,0,201);
  double denom_B_nPV20to30 = h2_nPV_JetMass_B->Integral(21,31,0,201);
  double denom_B_nPV30toInf = h2_nPV_JetMass_B->Integral(31,52,0,201);

  // nPV0to10 jets
  TGraph *g_eff_nPV0to10 = new TGraph(21);
  g_eff_nPV0to10->SetName("g_eff_nPV0to10");
  g_eff_nPV0to10->SetLineColor(kGreen+2);
  g_eff_nPV0to10->SetLineWidth(2);
  g_eff_nPV0to10->SetLineStyle(1);
  g_eff_nPV0to10->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(1,11,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(1,11,0,101-(1+i*5));

    g_eff_nPV0to10->SetPoint(i,(num_B/denom_B_nPV0to10),(num_S/denom_S_nPV0to10));
  }

  // nPV10to20 jets
  TGraph *g_eff_nPV10to20 = new TGraph(21);
  g_eff_nPV10to20->SetName("g_eff_nPV10to20");
  g_eff_nPV10to20->SetLineColor(kGreen+2);
  g_eff_nPV10to20->SetLineWidth(2);
  g_eff_nPV10to20->SetLineStyle(7);
  g_eff_nPV10to20->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(11,21,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(11,21,0,101-(1+i*5));

    g_eff_nPV10to20->SetPoint(i,(num_B/denom_B_nPV10to20),(num_S/denom_S_nPV10to20));
  }

  // nPV20to30 jets
  TGraph *g_eff_nPV20to30 = new TGraph(21);
  g_eff_nPV20to30->SetName("g_eff_nPV20to30");
  g_eff_nPV20to30->SetLineColor(kGreen+2);
  g_eff_nPV20to30->SetLineWidth(2);
  g_eff_nPV20to30->SetLineStyle(3);
  g_eff_nPV20to30->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(21,31,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(21,31,0,101-(1+i*5));

    g_eff_nPV20to30->SetPoint(i,(num_B/denom_B_nPV20to30),(num_S/denom_S_nPV20to30));
  }

  // nPV30toInf jets
  TGraph *g_eff_nPV30toInf = new TGraph(21);
  g_eff_nPV30toInf->SetName("g_eff_nPV30toInf");
  g_eff_nPV30toInf->SetLineColor(kGreen+2);
  g_eff_nPV30toInf->SetLineWidth(2);
  g_eff_nPV30toInf->SetLineStyle(8);
  g_eff_nPV30toInf->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(31,52,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(31,52,0,101-(1+i*5));

    g_eff_nPV30toInf->SetPoint(i,(num_B/denom_B_nPV30toInf),(num_S/denom_S_nPV30toInf));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_nPV0to10->Draw("L");
  g_eff_nPV10to20->Draw("L");
  g_eff_nPV20to30->Draw("L");
  g_eff_nPV30toInf->Draw("L");

  TLegend *legend = new TLegend(.55,.25,.85,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_nPV0to10, "0#leqnPV#leq10","l");
  legend->AddEntry(g_eff_nPV10to20, "10#leqnPV#leq20","l");
  legend->AddEntry(g_eff_nPV20to30, "20#leqnPV#leq30","l");
  legend->AddEntry(g_eff_nPV30toInf, "nPV#geq30","l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.14,0.97, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());
}


void efficiency_curves_bkgmc(const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
			     const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);

  // signal files
  TFile *file_S  = new TFile("output_files/WW500TrimmedJetMass.root");

  // background files
  TFile *file_B_1 = new TFile("output_files/QCDPythiaTrimmedJetMass.root");
  TFile *file_B_2 = new TFile("output_files/QCDHerwigTrimmedJetMass.root");

  // signal histograms
  TH2D *h2_nPV_JetMass_S = (TH2D*)file_S->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S = (TH2D*)file_S->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B_1 = (TH2D*)file_B_1->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_2 = (TH2D*)file_B_2->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B_1 = (TH2D*)file_B_1->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_2 = (TH2D*)file_B_2->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S = h2_nPV_JetMass_S->Integral(fPVLow,fPVHigh,0,201);

  // background denominator counts
  double denom_B_1 = h2_nPV_JetMass_B_1->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_2 = h2_nPV_JetMass_B_2->Integral(fPVLow,fPVHigh,0,201);

  // Pyhtia QCD
  TGraph *g_eff_1 = new TGraph(21);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_1->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_1->SetPoint(i,(num_B/denom_B_1),(num_S/denom_S));
  }

  // Herwig QCD
  TGraph *g_eff_2 = new TGraph(21);
  g_eff_2->SetName("g_eff_2");
  g_eff_2->SetLineColor(kRed);
  g_eff_2->SetLineWidth(2);
  g_eff_2->SetLineStyle(2);
  g_eff_2->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_2->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_2->SetPoint(i,(num_B/denom_B_2),(num_S/denom_S));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_1->Draw("L");
  g_eff_2->Draw("L");

  TLegend *legend = new TLegend(.55,.25,.85,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, "Pythia QCD","l");
  legend->AddEntry(g_eff_2, "Herwig QCD","l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.14,0.97, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());
}


void efficiency_curves_nsjgroomed(const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
			          const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);

  // signal files
  TFile *file_S_1  = new TFile("output_files/WW500TrimmedJetMass.root");
  TFile *file_S_2  = new TFile("output_files/WW500TrimmedJets.root");

  // background files
  TFile *file_B_1 = new TFile("output_files/QCDPythiaTrimmedJetMass.root");
  TFile *file_B_2 = new TFile("output_files/QCDPythiaTrimmedJets.root");

  // signal histograms
  TH2D *h2_nPV_JetMass_S_1 = (TH2D*)file_S_1->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_2 = (TH2D*)file_S_2->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S_1 = (TH2D*)file_S_1->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_2 = (TH2D*)file_S_2->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B_1 = (TH2D*)file_B_1->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_2 = (TH2D*)file_B_2->Get(("rutgersJetAnalyzer/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B_1 = (TH2D*)file_B_1->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_2 = (TH2D*)file_B_2->Get(("rutgersJetAnalyzer/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S_1 = h2_nPV_JetMass_S_1->Integral(fPVLow,fPVHigh,0,201);
  double denom_S_2 = h2_nPV_JetMass_S_2->Integral(fPVLow,fPVHigh,0,201);

  // background denominator counts
  double denom_B_1 = h2_nPV_JetMass_B_1->Integral(fPVLow,fPVHigh,0,201);
  double denom_B_2 = h2_nPV_JetMass_B_2->Integral(fPVLow,fPVHigh,0,201);

  // Default N-subjettiness
  TGraph *g_eff_1 = new TGraph(21);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_1->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_1->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_1->SetPoint(i,(num_B/denom_B_1),(num_S/denom_S_1));
  }

  // Trimmed N-subjettiness
  TGraph *g_eff_2 = new TGraph(21);
  g_eff_2->SetName("g_eff_2");
  g_eff_2->SetLineColor(kRed);
  g_eff_2->SetLineWidth(2);
  g_eff_2->SetLineStyle(2);
  g_eff_2->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S_2->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B_2->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff_2->SetPoint(i,(num_B/denom_B_2),(num_S/denom_S_2));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_1->Draw("L");
  g_eff_2->Draw("L");

  TLegend *legend = new TLegend(.55,.25,.85,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, "Default N-subj.","l");
  legend->AddEntry(g_eff_2, "Trimmed N-subj.","l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.05);
  l1.DrawLatex(0.14,0.97, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());
}


void makePlots()
{
  efficiency_curves_grooming("Pt500to700", 11, 21, "500<p_{T}<700 GeV, 10#leqnPV#leq20", 0, 0.3, "W_eff_grooming_nPV10to20_Pt500to700.eps");
  efficiency_curves_grooming("Pt500to700", 21, 31, "500<p_{T}<700 GeV, 20#leqnPV#leq30", 0, 0.3, "W_eff_grooming_nPV20to30_Pt500to700.eps");
  efficiency_curves_grooming("Pt700to900", 11, 21, "700<p_{T}<900 GeV, 10#leqnPV#leq20", 0, 0.4, "W_eff_grooming_nPV10to20_Pt700to900.eps");
  efficiency_curves_grooming("Pt700to900", 21, 31, "700<p_{T}<900 GeV, 20#leqnPV#leq30", 0, 0.4, "W_eff_grooming_nPV20to30_Pt700to900.eps");

  efficiency_curves_jetpt(11, 21, "Trimmed jet mass, 10#leqnPV#leq20", 0, 0.3, "W_eff_jetpt_nPV10to20_TrimmedJetMass.eps");
  efficiency_curves_jetpt(21, 31, "Trimmed jet mass, 20#leqnPV#leq30", 0, 0.3, "W_eff_jetpt_nPV20to30_TrimmedJetMass.eps");

  efficiency_curves_pileup("Pt500to700", "Trimmed jet mass, 500<p_{T}<700 GeV", 0, 0.3, "W_eff_pileup_Pt500to700_TrimmedJetMass.eps");
  efficiency_curves_pileup("Pt700to900", "Trimmed jet mass, 700<p_{T}<900 GeV", 0, 0.3, "W_eff_pileup_Pt700to900_TrimmedJetMass.eps");

  efficiency_curves_bkgmc("Pt500to700", 11, 21, "500<p_{T}<700 GeV, 10#leqnPV#leq20", 0, 0.3, "W_eff_bkgmc_nPV10to20_Pt500to700_TrimmedJetMass.eps");
  efficiency_curves_bkgmc("Pt500to700", 21, 31, "500<p_{T}<700 GeV, 20#leqnPV#leq30", 0, 0.3, "W_eff_bkgmc_nPV20to30_Pt500to700_TrimmedJetMass.eps");
  efficiency_curves_bkgmc("Pt700to900", 11, 21, "700<p_{T}<900 GeV, 10#leqnPV#leq20", 0, 0.3, "W_eff_bkgmc_nPV10to20_Pt700to900_TrimmedJetMass.eps");
  efficiency_curves_bkgmc("Pt700to900", 21, 31, "700<p_{T}<900 GeV, 20#leqnPV#leq30", 0, 0.3, "W_eff_bkgmc_nPV20to30_Pt700to900_TrimmedJetMass.eps");

  efficiency_curves_nsjgroomed("Pt500to700", 11, 21, "500<p_{T}<700 GeV, 10#leqnPV#leq20", 0, 0.3, "W_eff_nsjgroomed_nPV10to20_Pt500to700_TrimmedJetMass.eps");
  efficiency_curves_nsjgroomed("Pt500to700", 21, 31, "500<p_{T}<700 GeV, 20#leqnPV#leq30", 0, 0.3, "W_eff_nsjgroomed_nPV20to30_Pt500to700_TrimmedJetMass.eps");
  efficiency_curves_nsjgroomed("Pt700to900", 11, 21, "700<p_{T}<900 GeV, 10#leqnPV#leq20", 0, 0.3, "W_eff_nsjgroomed_nPV10to20_Pt700to900_TrimmedJetMass.eps");
  efficiency_curves_nsjgroomed("Pt700to900", 21, 31, "700<p_{T}<900 GeV, 20#leqnPV#leq30", 0, 0.3, "W_eff_nsjgroomed_nPV20to30_Pt700to900_TrimmedJetMass.eps");
}