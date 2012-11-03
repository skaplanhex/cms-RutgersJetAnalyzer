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
#include "TGraphAsymmErrors.h"
#include "exoStyle.C"


using namespace std;


void efficiency_curves_grooming(const string& fInputDir, const string& fFileS, const string& fFileB,
                                const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                                const string& fLeg1, const string& fLeg2, const string& fLeg3, const string& fLeg4,
                                const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S = new TFile((fInputDir + "/"+ fFileS).c_str());

  // background file
  TFile *file_B = new TFile((fInputDir + "/"+ fFileB).c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S          = (TH2D*)file_S->Get(("jetAnalyzerDefaultJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Filtered = (TH2D*)file_S->Get(("jetAnalyzerFilteredJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Pruned   = (TH2D*)file_S->Get(("jetAnalyzerPrunedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_Trimmed  = (TH2D*)file_S->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S          = (TH2D*)file_S->Get(("jetAnalyzerDefaultJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Filtered = (TH2D*)file_S->Get(("jetAnalyzerFilteredJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Pruned   = (TH2D*)file_S->Get(("jetAnalyzerPrunedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_Trimmed  = (TH2D*)file_S->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B          = (TH2D*)file_B->Get(("jetAnalyzerDefaultJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Filtered = (TH2D*)file_B->Get(("jetAnalyzerFilteredJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Pruned   = (TH2D*)file_B->Get(("jetAnalyzerPrunedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_Trimmed  = (TH2D*)file_B->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B          = (TH2D*)file_B->Get(("jetAnalyzerDefaultJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Filtered = (TH2D*)file_B->Get(("jetAnalyzerFilteredJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Pruned   = (TH2D*)file_B->Get(("jetAnalyzerPrunedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_Trimmed  = (TH2D*)file_B->Get(("jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S          = h2_nPV_JetMass_S->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_Filtered = h2_nPV_JetMass_S_Filtered->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_Pruned   = h2_nPV_JetMass_S_Pruned->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_Trimmed  = h2_nPV_JetMass_S_Trimmed->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B = h2_nPV_JetMass_B->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_Filtered = h2_nPV_JetMass_B_Filtered->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_Pruned   = h2_nPV_JetMass_B_Pruned->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_Trimmed  = h2_nPV_JetMass_B_Trimmed->Integral(fPVLow,fPVHigh,0,401);

  // Default jets
  TGraph *g_eff = new TGraph(21);
  g_eff->SetName("g_eff");
  g_eff->SetLineColor(kBlack);
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
  legend->AddEntry(g_eff_Trimmed, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_Filtered, fLeg2.c_str(),"l");
  legend->AddEntry(g_eff_Pruned, fLeg3.c_str(),"l");
  legend->AddEntry(g_eff, fLeg4.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff;
  delete g_eff_Filtered;
  delete g_eff_Pruned;
  delete g_eff_Trimmed;
  delete file_S;
  delete file_B;
}


void efficiency_curves_comp(const string& fInputDir, const string& fFileS, const string& fFileB, const string& fFileDir1, const string& fFileDir2,
                            const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                            const string& fLeg1, const string& fLeg2, const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S  = new TFile((fInputDir + "/"+ fFileS).c_str());

  // background file
  TFile *file_B = new TFile((fInputDir + "/"+ fFileB).c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S_1 = (TH2D*)file_S->Get((fFileDir1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_2 = (TH2D*)file_S->Get((fFileDir2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S_1 = (TH2D*)file_S->Get((fFileDir1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_2 = (TH2D*)file_S->Get((fFileDir2 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B_1 = (TH2D*)file_B->Get((fFileDir1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_2 = (TH2D*)file_B->Get((fFileDir2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B_1 = (TH2D*)file_B->Get((fFileDir1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_2 = (TH2D*)file_B->Get((fFileDir2 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S_1 = h2_nPV_JetMass_S_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_2 = h2_nPV_JetMass_S_2->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B_1 = h2_nPV_JetMass_B_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_2 = h2_nPV_JetMass_B_2->Integral(fPVLow,fPVHigh,0,401);

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
  legend->AddEntry(g_eff_1, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_2;
  delete g_eff_1;
  delete file_S;
  delete file_B;
}


void efficiency_curves_nsj_massdrop(const string& fInputDir, const string& fFileS, const string& fFileB, const string& fFileDir1, const string& fFileDir2,
                                    const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                                    const string& fLeg1, const string& fLeg2, const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S  = new TFile((fInputDir + "/"+ fFileS).c_str());

  // background file
  TFile *file_B = new TFile((fInputDir + "/"+ fFileB).c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S_1 = (TH2D*)file_S->Get((fFileDir1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_2 = (TH2D*)file_S->Get((fFileDir2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S = (TH2D*)file_S->Get((fFileDir1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_MassDrop_S = (TH2D*)file_S->Get((fFileDir2 + "/h2_nPV_MassDrop_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B_1 = (TH2D*)file_B->Get((fFileDir1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_2 = (TH2D*)file_B->Get((fFileDir2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B = (TH2D*)file_B->Get((fFileDir1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_MassDrop_B = (TH2D*)file_B->Get((fFileDir2 + "/h2_nPV_MassDrop_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S_1 = h2_nPV_JetMass_S_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_S_2 = h2_nPV_JetMass_S_2->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B_1 = h2_nPV_JetMass_B_1->Integral(fPVLow,fPVHigh,0,401);
  double denom_B_2 = h2_nPV_JetMass_B_2->Integral(fPVLow,fPVHigh,0,401);

  // Default N-subjettiness
  TGraph *g_eff_1 = new TGraph(21);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num_S = h2_nPV_tau2tau1_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_tau2tau1_B->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

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
    double num_S = h2_nPV_MassDrop_S->Integral(fPVLow,fPVHigh,0,101-(1+i*5));
    double num_B = h2_nPV_MassDrop_B->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

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
  legend->AddEntry(g_eff_1, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_2;
  delete g_eff_1;
  delete file_S;
  delete file_B;
}


void efficiency_vs_cut(const string& fInputFile, const string& fFileDir, const string& fVariable, const string& fPtRange,
                       const int fPVLow, const int fPVHigh, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                       const double fXmin, const double fXmax, const Double_t fYmin, const Double_t fYmax, const string& fOutputFile,
                       const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                       const Double_t fLeftMargin=0.14, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file  = new TFile(fInputFile.c_str());

  // histograms
  TH2D *h2_nPV_JetMass = (TH2D*)file->Get((fFileDir + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_variable = (TH2D*)file->Get((fFileDir + "/h2_nPV_" + fVariable + "_" + fPtRange).c_str());

  // denominator counts
  double denom = h2_nPV_JetMass->Integral(fPVLow,fPVHigh,0,401);

  // efficiency curve
  TGraph *g_eff = new TGraph(21);
  g_eff->SetName("g_eff");
  g_eff->SetLineColor(kGreen+2);
  g_eff->SetLineWidth(2);
  g_eff->SetLineStyle(1);
  g_eff->SetMarkerStyle(20);

  for(int i = 0; i<21; ++i)
  {
    double num = h2_nPV_variable->Integral(fPVLow,fPVHigh,0,101-(1+i*5));

    g_eff->SetPoint(i,(h2_nPV_variable->GetYaxis()->GetBinUpEdge(101-(1+i*5))),(num/denom));
  }

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  g_eff->Draw("L");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(fLeftMargin,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete bkg;
  delete c;
  delete g_eff;
  delete file;
}


void efficiency1D(const string& fInputFile, const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                  const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                  const string& fOutputFile, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                  const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file  = new TFile(fInputFile.c_str());

  // histograms
  TH1D *h1_Pass = (TH1D*)file->Get(fPlotPass.c_str());

  TH1D *h1_Total = (TH1D*)file->Get(fPlotTotal.c_str());

  h1_Pass->Rebin(fRebinX);
  h1_Total->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h1_Pass, h1_Total,"cp");
  g_efficiency->SetLineWidth(2);
  g_efficiency->SetLineColor(kBlue+2);
  g_efficiency->SetMarkerSize(1.);
  g_efficiency->SetMarkerStyle(24);
  g_efficiency->SetMarkerColor(kRed);

  g_efficiency->Draw("LP");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetTextSize(0.05);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin,0.97, fTitle.c_str());

  c->SetLogz();
  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete g_efficiency;
  delete bkg;
  delete c;
  delete file;
}


void efficiency1D_overlay(const string& fInputFile, const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                          const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                          const string& fOutputFile, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                          const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file  = new TFile(fInputFile.c_str());

  // histograms
  TH1D *h1_Pass_CSVL = (TH1D*)file->Get((fPlotPass + "_CSVL").c_str());
  TH1D *h1_Pass_CSVM = (TH1D*)file->Get((fPlotPass + "_CSVM").c_str());
  TH1D *h1_Pass_SubJetCSVL = (TH1D*)file->Get((fPlotPass + "_SubJetCSVL").c_str());
  TH1D *h1_Pass_SubJetCSVM = (TH1D*)file->Get((fPlotPass + "_SubJetCSVM").c_str());
  TH1D *h1_Pass_DoubleB = (TH1D*)file->Get((fPlotPass + "_DoubleB").c_str());

  TH1D *h1_Total = (TH1D*)file->Get(fPlotTotal.c_str());

  h1_Pass_CSVL->Rebin(fRebinX);
  h1_Pass_CSVM->Rebin(fRebinX);
  h1_Pass_SubJetCSVL->Rebin(fRebinX);
  h1_Pass_SubJetCSVM->Rebin(fRebinX);
  h1_Pass_DoubleB->Rebin(fRebinX);
  h1_Total->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  TGraphAsymmErrors *g_eff_CSVL = new TGraphAsymmErrors(h1_Pass_CSVL, h1_Total,"cp");
  g_eff_CSVL->SetLineWidth(2);
  g_eff_CSVL->SetLineColor(kBlue+2);
  g_eff_CSVL->SetMarkerSize(1.);
  g_eff_CSVL->SetMarkerStyle(20);
  g_eff_CSVL->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors *g_eff_CSVM = new TGraphAsymmErrors(h1_Pass_CSVM, h1_Total,"cp");
  g_eff_CSVM->SetLineWidth(2);
  g_eff_CSVM->SetLineColor(kGreen+2);
  g_eff_CSVM->SetMarkerSize(1.);
  g_eff_CSVM->SetMarkerStyle(21);
  g_eff_CSVM->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_eff_SubJetCSVL = new TGraphAsymmErrors(h1_Pass_SubJetCSVL, h1_Total,"cp");
  g_eff_SubJetCSVL->SetLineWidth(2);
  g_eff_SubJetCSVL->SetLineStyle(2);
  g_eff_SubJetCSVL->SetLineColor(kBlue+2);
  g_eff_SubJetCSVL->SetMarkerSize(1.);
  g_eff_SubJetCSVL->SetMarkerStyle(24);
  g_eff_SubJetCSVL->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors *g_eff_SubJetCSVM = new TGraphAsymmErrors(h1_Pass_SubJetCSVM, h1_Total,"cp");
  g_eff_SubJetCSVM->SetLineWidth(2);
  g_eff_SubJetCSVM->SetLineStyle(2);
  g_eff_SubJetCSVM->SetLineColor(kGreen+2);
  g_eff_SubJetCSVM->SetMarkerSize(1.);
  g_eff_SubJetCSVM->SetMarkerStyle(25);
  g_eff_SubJetCSVM->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_eff_DoubleB = new TGraphAsymmErrors(h1_Pass_DoubleB, h1_Total,"cp");
  g_eff_DoubleB->SetLineWidth(2);
  g_eff_DoubleB->SetLineColor(kRed+1);
  g_eff_DoubleB->SetMarkerSize(1.);
  g_eff_DoubleB->SetMarkerStyle(22);
  g_eff_DoubleB->SetMarkerColor(kRed+1);

  g_eff_CSVL->Draw("LP");
  g_eff_CSVM->Draw("LP");
  g_eff_SubJetCSVL->Draw("LP");
  g_eff_SubJetCSVM->Draw("LP");
  g_eff_DoubleB->Draw("LP");

  TLegend *legend = new TLegend(.15,.60,.35,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(g_eff_CSVL, "CSVL","lp");
  legend->AddEntry(g_eff_CSVM, "CSVM","lp");
  legend->AddEntry(g_eff_SubJetCSVL, "Subjet CSVL","lp");
  legend->AddEntry(g_eff_SubJetCSVM, "Subjet CSVM","lp");
  legend->AddEntry(g_eff_DoubleB, "DoubleB","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetTextSize(0.05);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin,0.97, fTitle.c_str());

  c->SetLogz();
  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  //delete legend;
  delete bkg;
  delete c;
  delete file;
}


void makePlots()
{
  //--------------------------------------------------------------------------------------------------------------------
  // W tagging
  efficiency_curves_comp("output_files_v2", "WW500_WTagging.root", "QCDPythia6_WTagging.root",
                         "jetAnalyzerTrimmedJetMass", "jetAnalyzerDefaultJetMass",
                         "Pt500toInf", 0, 52, "W, R=0.6, p_{T}>500 GeV", "Trimmed jet mass", "Default jet mass",
                         0, 0.3, "W_tag_eff_grooming_Pt500toInf_WW500.eps");

  efficiency_curves_nsj_massdrop("output_files_v2", "WW500_WTagging.root", "QCDPythia6_WTagging.root",
                                 "jetAnalyzerTrimmedJetMass", "jetAnalyzerCAPrunedJets",
                                 "Pt500toInf", 0, 52, "W, R=0.6, p_{T}>500 GeV", "N-subjettiness", "Mass drop",
                                 0, 0.3, "W_tag_eff_nsj_massdrop_Pt500toInf_WW500.eps");

  efficiency_vs_cut("output_files_v2/WW500_WTagging.root", "jetAnalyzerCAPrunedJets",
                    "MassDrop", "Pt500toInf", 0, 52, "W, CA R=0.6 pruned, 60<m<90 GeV, p_{T}>500 GeV",
                    "#mu=m_{subjet1}/m_{jet}<", "Tagging efficiency", 0, 1., 0, 0.9,
                    "W_tag_eff_vs_massdrop_cut_Pt500toInf_CAPrunedJets_WW500.eps", 1., 1.05);

  efficiency_vs_cut("output_files_v2/WW500_WTagging.root", "jetAnalyzerTrimmedJetMass",
                    "tau2tau1", "Pt500toInf", 0, 52, "W, anti-k_{T} R=0.6, 65<m<95 GeV (trimmed mass), p_{T}>500 GeV",
                    "#tau_{2}/#tau_{1}<", "Tagging efficiency", 0, 1., 0, 0.9,
                    "W_tag_eff_vs_nsj_cut_Pt500toInf_TrimmedJetMass_WW500.eps", 1., 1.05);

  //--------------------------------------------------------------------------------------------------------------------
  // Higgs tagging
  efficiency_curves_grooming("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "QCDPythia6_HiggsTagging.root",
                             "Pt300toInf", 0, 52, "Higgs, R=0.8, p_{T}>300 GeV", "Trimmed jet mass", "Filtered jet mass",
                             "Pruned jet mass", "Default jet mass", 0, 0.3, "H_tag_eff_grooming_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");


  efficiency_curves_comp("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "QCDPythia6_HiggsTagging.root",
                         "jetAnalyzerTrimmedJetMass", "jetAnalyzerTrimmedJets",
                         "Pt300toInf", 0, 52, "Higgs, R=0.8, p_{T}>300 GeV", "Default N-subj.", "Trimmed N-subj.",
                         0, 0.3, "H_tag_eff_nsjgroomed_Pt300toInf_TrimmedJetMass_BprimeBprimeToBHBHinc_M-800.eps");


  efficiency_curves_nsj_massdrop("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "QCDPythia6_HiggsTagging.root",
                                 "jetAnalyzerTrimmedJetMass", "jetAnalyzerCAPrunedJets",
                                 "Pt300toInf", 0, 52, "Higgs, R=0.8, p_{T}>300 GeV", "N-subjettiness", "Mass drop",
                                 0, 0.3, "H_tag_eff_nsj_massdrop_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");


  efficiency_vs_cut("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerCAPrunedJets",
                    "MassDrop", "Pt300toInf", 0, 52, "Higgs, CA R=0.8 pruned, 75<m<130 GeV, p_{T}>300 GeV",
                    "#mu=m_{subjet1}/m_{jet}<", "Tagging efficiency", 0, 1., 0, 0.9,
                    "H_tag_eff_vs_massdrop_cut_Pt300toInf_CAPrunedJets_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.05);

  efficiency_vs_cut("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass",
                    "tau2tau1", "Pt300toInf", 0, 52, "Higgs, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass), p_{T}>300 GeV",
                    "#tau_{2}/#tau_{1}<", "Tagging efficiency", 0, 1., 0, 0.9,
                    "H_tag_eff_vs_nsj_cut_Pt300toInf_TrimmedJetMass_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.05);

  // Higgs true Pt
  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
               "jetAnalyzerTrimmedJetMass/h1_BosonPt_Matched", "jetAnalyzerTrimmedJetMass/h1_BosonPt_DecaySel",
               "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, #DeltaR(H,jet)<0.5",
               "Higgs true p_{T} [GeV]", "Matching efficiency", 25, 0, 1000, 0, 1, "Matching_eff_dRHjet_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
               "jetAnalyzerTrimmedJetMass/h1_BosonPt_DecayProdMatched", "jetAnalyzerTrimmedJetMass/h1_BosonPt_DecaySel",
               "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, #DeltaR(b,jet)<0.8 & #DeltaR(#bar{b},jet)<0.8",
               "Higgs true p_{T} [GeV]", "Matching efficiency", 25, 0, 1000, 0, 1, "Matching_eff_dRbjet_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  // Jet Pt
  // Tagging efficiency
  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
               "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
               "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass)",
               "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  //efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
  //             "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonDecayProdMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonDecayProdMatched",
  //             "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass)",
  //             "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_HiggsToBBbar_BosonDecayProdMatched_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
                       "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
                       "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass)",
                       "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  //efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
  //                     "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass",
  //                     "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass)",
  //                     "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_JTACone_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
                       "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass_Nsubj", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
                       "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass)",
                       "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_Nsubj_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTaggingWithBTagging.root",
                       "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
                       "H#rightarrowb#bar{b}, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass)",
                       "Jet p_{T} [GeV]", "Higgs tagging efficiency", 40, 0, 1000, 0, 1, "Higgs_tag_eff_total_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);
  
  // Mistag rate
  efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTaggingWithBTagging.root",
                       "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
                       "QCD, anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass)",
                       "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_QCDPythia6.eps", 1., 0.9);
}