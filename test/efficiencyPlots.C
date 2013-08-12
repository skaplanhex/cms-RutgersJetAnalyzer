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


void efficiency_curves_grooming(const string& fFileS, const string& fFileB,
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
  TFile *file_S = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

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


void efficiency_curves_comp(const string& fFileS, const string& fFileB, const string& fFileDirS1, const string& fFileDirS2, const string& fFileDirB1, const string& fFileDirB2,
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
  TFile *file_S  = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S_1 = (TH2D*)file_S->Get((fFileDirS1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_S_2 = (TH2D*)file_S->Get((fFileDirS2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S_1 = (TH2D*)file_S->Get((fFileDirS1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_S_2 = (TH2D*)file_S->Get((fFileDirS2 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B_1 = (TH2D*)file_B->Get((fFileDirB1 + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH2D *h2_nPV_JetMass_B_2 = (TH2D*)file_B->Get((fFileDirB2 + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B_1 = (TH2D*)file_B->Get((fFileDirB1 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());
  TH2D *h2_nPV_tau2tau1_B_2 = (TH2D*)file_B->Get((fFileDirB2 + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

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

  TLegend *legend = new TLegend(.45,.25,.75,.45);
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
  l1.SetTextSize(0.035);
  l1.DrawLatex(0.06,0.96, fTitle.c_str());

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


void efficiency_curve(const string& fFileS, const string& fFileB, const string& fFileDirS, const string& fFileDirB,
                            const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
                            const string& fLeg, const double fXmin, const double fXmax, const string& fOutputFile)
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
  TFile *file_S  = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

  // signal histograms
  TH2D *h2_nPV_JetMass_S = (TH2D*)file_S->Get((fFileDirS + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_S = (TH2D*)file_S->Get((fFileDirS + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // background histograms
  TH2D *h2_nPV_JetMass_B = (TH2D*)file_B->Get((fFileDirB + "/h2_nPV_JetMass_" + fPtRange).c_str());

  TH2D *h2_nPV_tau2tau1_B = (TH2D*)file_B->Get((fFileDirB + "/h2_nPV_tau2tau1_" + fPtRange).c_str());

  // signal denominator counts
  double denom_S = h2_nPV_JetMass_S->Integral(fPVLow,fPVHigh,0,401);

  // background denominator counts
  double denom_B = h2_nPV_JetMass_B->Integral(fPVLow,fPVHigh,0,401);

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

    g_eff_1->SetPoint(i,(num_B/denom_B),(num_S/denom_S));
  }


  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg",";Mistag rate; Tagging efficiency",100,fXmin,fXmax,100,0,0.9);
  bkg->SetTitleOffset(1.05,"Y");
  bkg->Draw();

  g_eff_1->Draw("L");

  TLegend *legend = new TLegend(.45,.25,.75,.45);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, fLeg.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.035);
  l1.DrawLatex(0.06,0.96, fTitle.c_str());

  c->SetGridx();
  c->SetGridy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_1;
  delete file_S;
  delete file_B;
}


void efficiency_curves_nsj_massdrop(const string& fFileS, const string& fFileB, const string& fFileDir1, const string& fFileDir2,
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
  TFile *file_S  = new TFile(fFileS.c_str());

  // background file
  TFile *file_B = new TFile(fFileB.c_str());

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

  TLegend *legend = new TLegend(.45,.25,.75,.45);
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


void efficiency_curves_comp_xrange(const string& fFileS1, const string& fFileS2, const string& fFileB1, const string& fFileB2,
                                   const string& fPlotS1, const string& fPlotS2, const string& fPlotB1, const string& fPlotB2,
                                   const double fXMin, const double fXMax, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                                   const string& fLeg1, const string& fLeg2, const double fXmin, const double fXmax, const double fYmin, const double fYmax,
                                   const string& fOutputFile, const Int_t fLogy=0)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  // signal file
  TFile *file_S1  = new TFile(fFileS1.c_str());
  TFile *file_S2  = new TFile(fFileS2.c_str());

  // background file
  TFile *file_B1 = new TFile(fFileB1.c_str());
  TFile *file_B2 = new TFile(fFileB2.c_str());

  // signal histograms
  TH2D *h2_S_1 = (TH2D*)file_S1->Get(fPlotS1.c_str());
  TH2D *h2_S_2 = (TH2D*)file_S2->Get(fPlotS2.c_str());

  // background histograms
  TH2D *h2_B_1 = (TH2D*)file_B1->Get(fPlotB1.c_str());
  TH2D *h2_B_2 = (TH2D*)file_B2->Get(fPlotB2.c_str());

  // signal denominator counts
  double denom_S_1 = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),0,101);
  double denom_S_2 = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),0,101);

  // background denominator counts
  double denom_B_1 = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),0,101);
  double denom_B_2 = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),0,101);


  TGraph *g_eff_1 = new TGraph(29);
  g_eff_1->SetName("g_eff_1");
  g_eff_1->SetLineColor(kGreen+2);
  g_eff_1->SetLineWidth(2);
  g_eff_1->SetLineStyle(1);
  g_eff_1->SetMarkerStyle(20);

  for(int i = 0; i<5; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),101-i,101);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),101-i,101);

    g_eff_1->SetPoint(i,(num_S/denom_S_1),(num_B/denom_B_1));
  }
  for(int i = 1; i<20; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),101-(i*5),101);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),101-(i*5),101);

    g_eff_1->SetPoint(i+4,(num_S/denom_S_1),(num_B/denom_B_1));
  }
  for(int i = 1; i<6; ++i)
  {
    double num_S = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),6-i,101);
    double num_B = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),6-i,101);

    g_eff_1->SetPoint(i+23,(num_S/denom_S_1),(num_B/denom_B_1));
  }

  TGraph *g_eff_2 = new TGraph(29);
  g_eff_2->SetName("g_eff_2");
  g_eff_2->SetLineColor(kRed);
  g_eff_2->SetLineWidth(2);
  g_eff_2->SetLineStyle(2);
  g_eff_2->SetMarkerStyle(20);

  for(int i = 0; i<5; ++i)
  {
    double num_S = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),101-i,101);
    double num_B = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),101-i,101);

    g_eff_2->SetPoint(i,(num_S/denom_S_2),(num_B/denom_B_2));
  }
  for(int i = 1; i<20; ++i)
  {
    double num_S = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),101-(i*5),101);
    double num_B = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),101-(i*5),101);

    g_eff_2->SetPoint(i+4,(num_S/denom_S_2),(num_B/denom_B_2));
  }
  for(int i = 1; i<6; ++i)
  {
    double num_S = h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),6-i,101);
    double num_B = h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),6-i,101);

    g_eff_2->SetPoint(i+23,(num_S/denom_S_2),(num_B/denom_B_2));
  }

  // CSV loose operating point
  TGraph *g_eff_L = new TGraph(2);
  g_eff_L->SetName("g_eff_L");
  g_eff_L->SetMarkerStyle(31);
  g_eff_L->SetMarkerSize(1.5);

  g_eff_L->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(0.244),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(0.244),101)/denom_B_1));
  g_eff_L->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(0.244),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(0.244),101)/denom_B_2));

  // CSV medium operating point
  TGraph *g_eff_M = new TGraph(2);
  g_eff_M->SetName("g_eff_L");
  g_eff_M->SetMarkerStyle(27);
  g_eff_M->SetMarkerSize(1.5);

  g_eff_M->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(0.679),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(0.679),101)/denom_B_1));
  g_eff_M->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(0.679),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(0.679),101)/denom_B_2));

  // CSV tight operating point
  TGraph *g_eff_T = new TGraph(2);
  g_eff_T->SetName("g_eff_L");
  g_eff_T->SetMarkerStyle(30);
  g_eff_T->SetMarkerSize(1.5);

  g_eff_T->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(0.898),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(0.898),101)/denom_B_1));
  g_eff_T->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(0.898),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(0.898),101)/denom_B_2));


  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(1.1,"Y");
  bkg->Draw();
  c->SetGridx();
  c->SetGridy();

  g_eff_1->Draw("L");
  g_eff_2->Draw("L");

  g_eff_L->Draw("P");
  g_eff_M->Draw("P");
  g_eff_T->Draw("P");

  TLegend *legend = new TLegend(.16,.64,.36,.77);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_eff_1, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLegend *legend2 = new TLegend(.16,.45,.36,.60);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetTextFont(42);
  legend2->SetTextSize(0.03);
  legend2->AddEntry(g_eff_L, "Loose","p");
  legend2->AddEntry(g_eff_M, "Medium","p");
  legend2->AddEntry(g_eff_T, "Tight","p");
  legend2->Draw();

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14+0.03,0.90, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.045);
  l1.SetTextFont(62);
  l1.DrawLatex(0.14,0.96, "CMS Preliminary Simulation, #sqrt{s} = 8 TeV");
  //l1.DrawLatex(0.14,0.97, "CMS Simulation");
  //l1.SetTextFont(42);
  //l1.DrawLatex(0.14+0.40,0.97, "#sqrt{s} = 8 TeV");
  l1.SetTextFont(42);
  l1.SetTextSize(0.04);
  if(fOutputFile.find("JTA")!=string::npos) l1.DrawLatex(0.48,0.18, "JTA = jet-track association");


  if(fLogy) c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete g_eff_2;
  delete g_eff_1;
  delete file_S1;
  delete file_S2;
  delete file_B1;
  delete file_B2;
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

  TCanvas *c = new TCanvas("c", "",1000,800);
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
  gStyle->SetGridColor(kGray);
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
  c->SetGridx();
  c->SetGridy();

  TGraphAsymmErrors *g_efficiency = new TGraphAsymmErrors(h1_Pass, h1_Total,"cp");
  g_efficiency->SetLineWidth(2);
  g_efficiency->SetLineColor(kBlue+2);
  g_efficiency->SetMarkerSize(1.);
  g_efficiency->SetMarkerStyle(24);
  g_efficiency->SetMarkerColor(kBlue+2);

  g_efficiency->Draw("LP");

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetTextSize(0.045);
  l1.SetNDC();
  if( fOutputFile.find("Matching_eff_dR")!=string::npos ) l1.DrawLatex(fLeftMargin+0.03,0.21, fTitle.c_str());
  else                                                    l1.DrawLatex(fLeftMargin+0.03,0.90, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.05);
  l1.SetTextFont(62);
  l1.DrawLatex(fLeftMargin,0.97, "CMS Preliminary Simulation, #sqrt{s} = 8 TeV");
  //l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation");
  //l1.SetTextFont(42);
  //l1.DrawLatex(fLeftMargin+0.35,0.97, "#sqrt{s} = 8 TeV");

  c->SetLogz();
  c->SaveAs(fOutputFile.c_str());

  delete g_efficiency;
  delete bkg;
  delete c;
  delete file;
}


void efficiency1D_overlay(const string& fInputFile, const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                          const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                          const string& fOutputFile, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                          const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8, const string& fSubJetMode="k_{T}")
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
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
  c->SetGridx();
  c->SetGridy();

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

  TLegend *legend = new TLegend(.15+(fLeftMargin-0.12),.67,.35+(fLeftMargin-0.12),.92);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(g_eff_CSVL, "CSVL","lp");
  legend->AddEntry(g_eff_CSVM, "CSVM","lp");
  legend->AddEntry(g_eff_SubJetCSVL, ("Subjet CSVL ("+ fSubJetMode +")").c_str(),"lp");
  legend->AddEntry(g_eff_SubJetCSVM, ("Subjet CSVM ("+ fSubJetMode +")").c_str(),"lp");
  legend->AddEntry(g_eff_DoubleB, "DoubleB","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetTextSize(0.05);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin,0.97, fTitle.c_str());

  c->SetLogz();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete bkg;
  delete c;
  delete file;
}

void efficiency1D_overlayMulti_3(const string& fInputFile1, const string& fInputFile2, const string& fInputFile3,
                               const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                               const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                               const string& fOutputFile, const Int_t fLogy=0, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                               const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file1  = new TFile(fInputFile1.c_str());
  TFile *file2  = new TFile(fInputFile2.c_str());
  TFile *file3  = new TFile(fInputFile3.c_str());

  // histograms
  TH1D *h1_Pass1 = (TH1D*)file1->Get(fPlotPass.c_str());
  TH1D *h1_Pass2 = (TH1D*)file2->Get(fPlotPass.c_str());
  TH1D *h1_Pass3 = (TH1D*)file3->Get(fPlotPass.c_str());

  TH1D *h1_Total1 = (TH1D*)file1->Get(fPlotTotal.c_str());
  TH1D *h1_Total2 = (TH1D*)file2->Get(fPlotTotal.c_str());
  TH1D *h1_Total3 = (TH1D*)file3->Get(fPlotTotal.c_str());

  h1_Pass1->Rebin(fRebinX);
  h1_Pass2->Rebin(fRebinX);
  h1_Pass3->Rebin(fRebinX);

  h1_Total1->Rebin(fRebinX);
  h1_Total2->Rebin(fRebinX);
  h1_Total3->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();
  c->SetGridx();
  c->SetGridy();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  TGraphAsymmErrors *g_efficiency1 = new TGraphAsymmErrors(h1_Pass1, h1_Total1,"cp");
  g_efficiency1->SetLineWidth(2);
  g_efficiency1->SetLineStyle(1);
  g_efficiency1->SetLineColor(kGreen+2);
  g_efficiency1->SetMarkerSize(1.);
  g_efficiency1->SetMarkerStyle(24);
  g_efficiency1->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_efficiency2 = new TGraphAsymmErrors(h1_Pass2, h1_Total2,"cp");
  g_efficiency2->SetLineWidth(2);
  g_efficiency2->SetLineStyle(2);
  g_efficiency2->SetLineColor(kRed);
  g_efficiency2->SetMarkerSize(1.);
  g_efficiency2->SetMarkerStyle(26);
  g_efficiency2->SetMarkerColor(kRed);

  TGraphAsymmErrors *g_efficiency3 = new TGraphAsymmErrors(h1_Pass3, h1_Total3,"cp");
  g_efficiency3->SetLineWidth(2);
  g_efficiency3->SetLineStyle(3);
  g_efficiency3->SetLineColor(kBlue+2);
  g_efficiency3->SetMarkerSize(1.);
  g_efficiency3->SetMarkerStyle(20);
  g_efficiency3->SetMarkerColor(kBlue+2);

  g_efficiency1->Draw("LP");
  g_efficiency2->Draw("LPsame");
  g_efficiency3->Draw("LPsame");

  TLegend *legend = new TLegend(.15+(fLeftMargin-0.12),.67,.35+(fLeftMargin-0.12),.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.045);
  legend->AddEntry(g_efficiency3, "106<m_{jet}<135 GeV/c^{2} (pruned)","lp");
  legend->AddEntry(g_efficiency1, "75<m_{jet}<135 GeV/c^{2} (pruned)","lp");
  legend->AddEntry(g_efficiency2, "75<m_{jet}<106 GeV/c^{2} (pruned)","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetTextSize(0.05);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin+0.03,0.27, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.05);
  l1.SetTextFont(62);
  l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation");
  l1.SetTextFont(42);
  l1.DrawLatex(fLeftMargin+0.35,0.97, "#sqrt{s} = 8 TeV");

  //c->RedrawAxis();
  c->SetLogz();
  if(fLogy) c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete g_efficiency1;
  delete g_efficiency2;
  delete bkg;
  delete c;
  delete file1;
  delete file2;
}


void efficiency1D_overlayMulti_5(const string& fInputFile1, const string& fInputFile2, const string& fInputFile3, const string& fInputFile4, const string& fInputFile5,
                               const string& fPlotPass, const string& fPlotTotal, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                               const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                               const string& fOutputFile, const Int_t fLogy=0, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                               const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  // input files
  TFile *file1  = new TFile(fInputFile1.c_str());
  TFile *file2  = new TFile(fInputFile2.c_str());
  TFile *file3  = new TFile(fInputFile3.c_str());
  TFile *file4  = new TFile(fInputFile4.c_str());
  TFile *file5  = new TFile(fInputFile5.c_str());

  // histograms
  TH1D *h1_Pass1 = (TH1D*)file1->Get(fPlotPass.c_str());
  TH1D *h1_Pass2 = (TH1D*)file2->Get(fPlotPass.c_str());
  TH1D *h1_Pass3 = (TH1D*)file3->Get(fPlotPass.c_str());
  TH1D *h1_Pass4 = (TH1D*)file4->Get(fPlotPass.c_str());
  TH1D *h1_Pass5 = (TH1D*)file5->Get(fPlotPass.c_str());

  TH1D *h1_Total1 = (TH1D*)file1->Get(fPlotTotal.c_str());
  TH1D *h1_Total2 = (TH1D*)file2->Get(fPlotTotal.c_str());
  TH1D *h1_Total3 = (TH1D*)file3->Get(fPlotTotal.c_str());
  TH1D *h1_Total4 = (TH1D*)file4->Get(fPlotTotal.c_str());
  TH1D *h1_Total5 = (TH1D*)file5->Get(fPlotTotal.c_str());

  h1_Pass1->Rebin(fRebinX);
  h1_Pass2->Rebin(fRebinX);
  h1_Pass3->Rebin(fRebinX);
  h1_Pass4->Rebin(fRebinX);
  h1_Pass5->Rebin(fRebinX);

  h1_Total1->Rebin(fRebinX);
  h1_Total2->Rebin(fRebinX);
  h1_Total3->Rebin(fRebinX);
  h1_Total4->Rebin(fRebinX);
  h1_Total5->Rebin(fRebinX);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();
  c->SetGridx();
  c->SetGridy();

  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  TGraphAsymmErrors *g_efficiency1 = new TGraphAsymmErrors(h1_Pass1, h1_Total1,"cp");
  g_efficiency1->SetLineWidth(2);
  g_efficiency1->SetLineStyle(1);
  g_efficiency1->SetLineColor(kGreen+2);
  g_efficiency1->SetMarkerSize(1.);
  g_efficiency1->SetMarkerStyle(24);
  g_efficiency1->SetMarkerColor(kGreen+2);

  TGraphAsymmErrors *g_efficiency2 = new TGraphAsymmErrors(h1_Pass2, h1_Total2,"cp");
  g_efficiency2->SetLineWidth(2);
  g_efficiency2->SetLineStyle(2);
  g_efficiency2->SetLineColor(kRed);
  g_efficiency2->SetMarkerSize(1.);
  g_efficiency2->SetMarkerStyle(26);
  g_efficiency2->SetMarkerColor(kRed);

  TGraphAsymmErrors *g_efficiency3 = new TGraphAsymmErrors(h1_Pass3, h1_Total3,"cp");
  g_efficiency3->SetLineWidth(2);
  g_efficiency3->SetLineStyle(3);
  g_efficiency3->SetLineColor(kBlue+2);
  g_efficiency3->SetMarkerSize(1.);
  g_efficiency3->SetMarkerStyle(20);
  g_efficiency3->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors *g_efficiency4 = new TGraphAsymmErrors(h1_Pass4, h1_Total4,"cp");
  g_efficiency4->SetLineWidth(2);
  g_efficiency4->SetLineStyle(4);
  g_efficiency4->SetLineColor(kOrange+1);
  g_efficiency4->SetMarkerSize(1.);
  g_efficiency4->SetMarkerStyle(22);
  g_efficiency4->SetMarkerColor(kOrange+1);

  TGraphAsymmErrors *g_efficiency5 = new TGraphAsymmErrors(h1_Pass5, h1_Total5,"cp");
  g_efficiency5->SetLineWidth(2);
  g_efficiency5->SetLineStyle(5);
  g_efficiency5->SetLineColor(kMagenta+2);
  g_efficiency5->SetMarkerSize(1.);
  g_efficiency5->SetMarkerStyle(23);
  g_efficiency5->SetMarkerColor(kMagenta+2);

  g_efficiency1->Draw("LP");
  g_efficiency2->Draw("LPsame");
  g_efficiency3->Draw("LPsame");
  g_efficiency4->Draw("LPsame");
  g_efficiency5->Draw("LPsame");

  TLegend *legend = new TLegend(.15+(fLeftMargin-0.12),.55,.35+(fLeftMargin-0.12),.80);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.05);
  legend->AddEntry(g_efficiency1, "H(120)#rightarrowb#bar{b}","lp");
  legend->AddEntry(g_efficiency4, "Z","lp");
  legend->AddEntry(g_efficiency5, "top","lp");
  legend->AddEntry(g_efficiency3, "W","lp");
  legend->AddEntry(g_efficiency2, "QCD","lp");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetTextSize(0.045);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin+0.03,0.26, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.05);
  l1.SetTextFont(62);
  l1.DrawLatex(fLeftMargin,0.97, "CMS Preliminary Simulation, #sqrt{s} = 8 TeV");
  //l1.DrawLatex(fLeftMargin,0.97, "CMS Simulation");
  //l1.SetTextFont(42);
  //l1.DrawLatex(fLeftMargin+0.35,0.97, "#sqrt{s} = 8 TeV");

  //c->RedrawAxis();
  c->SetLogz();
  if(fLogy) c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete g_efficiency1;
  delete g_efficiency2;
  delete bkg;
  delete c;
  delete file1;
  delete file2;
}


void makePlots()
{
  //--------------------------------------------------------------------------------------------------------------------
  //================
  // W tagging
  //================
//   efficiency_curves_grooming("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
//                             "Pt500toInf", 0, 52, "WW vs QCD, AK R=0.6, p_{T}>500 GeV", "Trimmed jet mass", "Filtered jet mass",
//                             "Pruned jet mass", "Default jet mass", 0, 0.3, "W_tag_eff_grooming_Pt500toInf_WW500.eps");
//
//   efficiency_curves_comp("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
//                          "jetAnalyzerTrimmedJetMass", "jetAnalyzerTrimmedJets", "jetAnalyzerTrimmedJetMass", "jetAnalyzerTrimmedJets",
//                          "Pt500toInf", 0, 52, "WW vs QCD, AK R=0.6, p_{T}>500 GeV", "Default N-subjettiness", "Trimmed N-subjettiness",
//                          0, 0.3, "W_tag_eff_nsjgroomed_Pt500toInf_TrimmedJetMass_WW500.eps");
//
//   efficiency_curves_comp("output_files_v2/WW500_WTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_WTagging_JetSubstructure.root",
//                          "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJets", "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJets",
//                          "Pt500toInf", 0, 52, "WW vs QCD, AK R=0.8, p_{T}>500 GeV, 55<m<95 GeV", "Default N-subjettiness", "Pruned N-subjettiness",
//                          0, 0.3, "W_tag_eff_nsjgroomed_AK_Pt500toInf_PrunedJetMass_WW500.eps");
//
//   efficiency_curves_comp("output_files_v2/WW500_WTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_WTagging_JetSubstructure.root",
//                          "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJets", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJets",
//                          "Pt500toInf", 0, 52, "WW vs QCD, CA R=0.8, p_{T}>500 GeV, 55<m<95 GeV", "Default N-subjettiness", "Pruned N-subjettiness",
//                          0, 0.3, "W_tag_eff_nsjgroomed_CA_Pt500toInf_PrunedJetMass_WW500.eps");
//
//   efficiency_curves_comp("output_files_v2/WW500_WTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_WTagging_JetSubstructure.root",
//                          "jetAnalyzerPrunedJetMass", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerPrunedJetMass", "jetAnalyzerCAPrunedJetMass",
//                          "Pt500toInf", 0, 52, "WW vs QCD, R=0.8, p_{T}>500 GeV, 55<m<95 GeV", "N-subjettiness (AK)", "N-subjettiness (CA)",
//                          0, 0.3, "W_tag_eff_AK_vs_CA_Pt500toInf_PrunedJetMass_WW500.eps");
//
//   efficiency_curves_comp("output_files_v2/WW500_WTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_WTagging_JetSubstructure.root",
//                          "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJetMassKtAxes", "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJetMassKtAxes",
//                          "Pt500toInf", 0, 52, "WW vs QCD, AK R=0.8, p_{T}>500 GeV, 55<m<95 GeV", "N-subj. (Onepass k_{T} axes)", "N-subj. (k_{T} axes)",
//                          0, 0.3, "W_tag_eff_AK_onepassktaxes_vs_ktaxes_Pt500toInf_PrunedJetMass_WW500.eps");
//
//   efficiency_curves_nsj_massdrop("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
//                                  "jetAnalyzerTrimmedJetMass", "jetAnalyzerCAPrunedJets",
//                                  "Pt500toInf", 0, 52, "WW vs QCD, R=0.6, p_{T}>500 GeV", "N-subj. (trimmed AK)", "Mass drop (pruned CA)",
//                                  0, 0.3, "W_tag_eff_nsj_massdrop_Pt500toInf_WW500.eps");
//
//   efficiency_curves_nsj_massdrop("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
//                                  "jetAnalyzerPrunedJetMass", "jetAnalyzerCAPrunedJets",
//                                  "Pt500toInf", 0, 52, "WW vs QCD, R=0.6, p_{T}>500 GeV", "N-subj. (pruned AK)", "Mass drop (pruned CA)",
//                                  0, 0.3, "W_tag_eff_nsj_massdrop_pruned_Pt500toInf_WW500.eps");
//
//   efficiency_vs_cut("output_files_v2/WW500_WTagging.root", "jetAnalyzerCAPrunedJets",
//                     "MassDrop", "Pt500toInf", 0, 52, "WW, CA R=0.6 pruned, p_{T}>500 GeV, 60<m<90 GeV",
//                     "#mu=m_{subjet1}/m_{jet}<", "Tagging efficiency", 0, 1., 0, 0.9,
//                     "W_tag_eff_vs_massdrop_cut_Pt500toInf_CAPrunedJets_WW500.eps", 0.9, 1.05);
//
//   efficiency_vs_cut("output_files_v2/WW500_WTagging.root", "jetAnalyzerTrimmedJetMass",
//                     "tau2tau1", "Pt500toInf", 0, 52, "WW, AK R=0.6, p_{T}>500 GeV, 65<m<95 GeV (trimmed)",
//                     "#tau_{2}/#tau_{1}<", "Tagging efficiency", 0, 1., 0, 0.9,
//                     "W_tag_eff_vs_nsj_cut_Pt500toInf_TrimmedJetMass_WW500.eps", 0.9, 1.05);

  //--------------------------------------------------------------------------------------------------------------------
  //================
  // Higgs tagging
  //================
//   efficiency_curves_grooming("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
//                              "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, AK R=0.8, p_{T}>300 GeV", "Trimmed jet mass", "Filtered jet mass",
//                              "Pruned jet mass", "Default jet mass", 0, 0.3, "H_tag_eff_grooming_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");
//
//
//   efficiency_curves_comp("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
//                          "jetAnalyzerTrimmedJetMass", "jetAnalyzerTrimmedJets", "jetAnalyzerTrimmedJetMass", "jetAnalyzerTrimmedJets",
//                          "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, AK R=0.8, p_{T}>300 GeV", "Default N-subjettiness", "Trimmed N-subjettiness",
//                          0, 0.3, "H_tag_eff_nsjgroomed_Pt300toInf_TrimmedJetMass_BprimeBprimeToBHBHinc_M-800.eps");

//   efficiency_curves_comp("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_HiggsTagging_JetSubstructure.root",
//                          "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJets", "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJets",
//                          "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV", "Default N-subjettiness", "Pruned N-subjettiness",
//                          0, 0.3, "H_tag_eff_nsjgroomed_AK_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-800.eps");
//
//   efficiency_curves_comp("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_HiggsTagging_JetSubstructure.root",
//                          "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJets", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJets",
//                          "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV", "Default N-subjettiness", "Pruned N-subjettiness",
//                          0, 0.3, "H_tag_eff_nsjgroomed_CA_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-800.eps");
//
//   efficiency_curves_comp("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_HiggsTagging_JetSubstructure.root",
//                          "jetAnalyzerPrunedJetMass", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerPrunedJetMass", "jetAnalyzerCAPrunedJetMass",
//                          "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, R=0.8, p_{T}>300 GeV, 75<m<135 GeV", "N-subjettiness (AK)", "N-subjettiness (CA)",
//                          0, 0.3, "H_tag_eff_AK_vs_CA_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-800.eps");
//
//   efficiency_curves_comp("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_HiggsTagging_JetSubstructure.root",
//                          "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJetMassKtAxes", "jetAnalyzerPrunedJetMass", "jetAnalyzerPrunedJetMassKtAxes",
//                          "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV", "N-subj. (Onepass k_{T} axes)", "N-subj. (k_{T} axes)",
//                          0, 0.3, "H_tag_eff_AK_onepassktaxes_vs_ktaxes_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-800.eps");

  // N-subjettiness for inclusive QCD compared with udsg and GSP b jets
//   efficiency_curves_comp("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor.root",
//                          "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJetMass_bQuarksGSP",
//                          "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, CA R=0.8, p_{T}>300 GeV, 75<m<135 GeV (pruned)", "Inclusive QCD", "GSP b jets",
//                          0, 0.3, "Nsubjettiness_eff_vs_mistag_InclQCD_bGSP_CA8_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-1000.eps");
//
//   efficiency_curves_comp("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor.root",
//                          "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJetMass_udsQuarks_g",
//                          "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, CA R=0.8, p_{T}>300 GeV, 75<m<135 GeV (pruned)", "Inclusive QCD", "udsg jets",
//                          0, 0.3, "Nsubjettiness_eff_vs_mistag_InclQCD_udsg_CA8_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-1000.eps");

//   efficiency_curve("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root",
//                          "jetAnalyzerCAPrunedJetMass", "jetAnalyzerCAPrunedJetMass", "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, CA R=0.8, p_{T}>300 GeV, 75<m<135 GeV (pruned)", "N-subjettiness",
//                          0, 0.3, "H_tag_eff_vs_mistag_CA_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-1000.eps");

//
//
//   efficiency_curves_nsj_massdrop("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
//                                  "jetAnalyzerTrimmedJetMass", "jetAnalyzerCAPrunedJets",
//                                  "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, R=0.8, p_{T}>300 GeV", "N-subj. (trimmed AK)", "Mass drop (pruned CA)",
//                                  0, 0.3, "H_tag_eff_nsj_massdrop_trimmedAK_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");
//
//   efficiency_curves_nsj_massdrop("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
//                                  "jetAnalyzerPrunedJetMass", "jetAnalyzerCAPrunedJets",
//                                  "Pt300toInf", 0, 52, "H#rightarrowb#bar{b} vs QCD, R=0.8, p_{T}>300 GeV", "N-subj. (pruned AK)", "Mass drop (pruned CA)",
//                                  0, 0.3, "H_tag_eff_nsj_massdrop_prunedAK_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");
//
//   efficiency_vs_cut("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerCAPrunedJets",
//                     "MassDrop", "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, CA R=0.8 pruned, p_{T}>300 GeV, 75<m<135 GeV",
//                     "#mu=m_{subjet1}/m_{jet}<", "Tagging efficiency", 0, 1., 0, 0.9,
//                     "H_tag_eff_vs_massdrop_cut_Pt300toInf_CAPrunedJets_BprimeBprimeToBHBHinc_M-800.eps", 0.9, 1.05);
//
//   efficiency_vs_cut("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass",
//                     "tau2tau1", "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV (trimmed)",
//                     "#tau_{2}/#tau_{1}<", "Tagging efficiency", 0, 1., 0, 0.9,
//                     "H_tag_eff_vs_nsj_cut_Pt300toInf_TrimmedJetMass_BprimeBprimeToBHBHinc_M-800.eps", 0.9, 1.05);

//   efficiency_vs_cut("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_JetSubstructure.root", "jetAnalyzerPrunedJetMass",
//                     "tau2tau1", "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV",
//                     "#tau_{2}/#tau_{1}<", "Tagging efficiency", 0, 1., 0, 0.9,
//                     "H_tag_eff_vs_nsj_cut_AK_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-800.eps", 0.9, 1.05);
//
//   efficiency_vs_cut("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_JetSubstructure.root", "jetAnalyzerPrunedJets",
//                     "tau2tau1", "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV",
//                     "Pruned #tau_{2}/#tau_{1}<", "Tagging efficiency", 0, 1., 0, 0.9,
//                     "H_tag_eff_vs_prunednsj_cut_AK_Pt300toInf_PrunedJetMass_BprimeBprimeToBHBHinc_M-800.eps", 0.9, 1.05);
//
//   efficiency_vs_cut("output_files_v2/QCDPythia6_HiggsTagging_JetSubstructure.root", "jetAnalyzerPrunedJetMass",
//                     "tau2tau1", "Pt300toInf", 0, 52, "QCD, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV",
//                     "#tau_{2}/#tau_{1}<", "Mistag rate", 0, 1., 0, 0.16,
//                     "H_mistag_rate_vs_nsj_cut_AK_Pt300toInf_PrunedJetMass_QCDPythia6.eps", 0.9, 1.05);
//
//   efficiency_vs_cut("output_files_v2/QCDPythia6_HiggsTagging_JetSubstructure.root", "jetAnalyzerPrunedJets",
//                     "tau2tau1", "Pt300toInf", 0, 52, "QCD, AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV",
//                     "Pruned #tau_{2}/#tau_{1}<", "Mistag rate", 0, 1., 0, 0.16,
//                     "H_mistag_rate_vs_prunednsj_cut_AK_Pt300toInf_PrunedJetMass_QCDPythia6.eps", 0.9, 1.05);

  // Higgs true Pt
  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
               "jetAnalyzerCAPrunedJetMass/h1_BosonPt_Matched", "jetAnalyzerCAPrunedJetMass/h1_BosonPt_DecaySel",
               "H(120)#rightarrowb#bar{b}, CA R=0.8, #DeltaR(H,jet)<0.5",
               "Higgs true p_{T} [GeV/c]", "Matching efficiency", 10, 0, 1000, 0, 1, "Matching_eff_dRHjet_BprimeBprimeToBHBHinc_M-1000.eps", 1., 0.9);

  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
               "jetAnalyzerCAPrunedJetMass/h1_BosonPt_DecayProdMatched", "jetAnalyzerCAPrunedJetMass/h1_BosonPt_DecaySel",
               "H(120)#rightarrowb#bar{b}, CA R=0.8, #DeltaR(b,jet)<0.8 & #DeltaR(#bar{b},jet)<0.8",
               "Higgs true p_{T} [GeV/c]", "Matching efficiency", 10, 0, 1000, 0, 1, "Matching_eff_dRbjet_BprimeBprimeToBHBHinc_M-1000.eps", 1., 0.9);

  //--------------------------------------------------------------------------------------------------------------------

//   // Jet mass cut efficiency for trimmed jet mass
//   efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
//                "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);
//
//   efficiency1D("output_files_v2/QCDPythia6_HiggsTagging.root",
//                "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
//                "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_QCDPythia6.eps", 1., 0.9);

//   // Jet mass cut efficiency for pruned jet mass
//   // AK8
//   efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched",
//                "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_PrunedJetMass_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);
//
//   efficiency1D("output_files_v2/QCDPythia6_HiggsTagging.root",
//                "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched",
//                "QCD, AK R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_PrunedJetMass_QCDPythia6.eps", 1., 0.9);

  // CA8
  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
               "#splitline{H(120)#rightarrowb#bar{b}, CA R=0.8, #DeltaR(H,jet)<0.5}{75<m_{jet}<135 GeV/c^{2} (pruned)}",
               "Fat jet p_{T} [GeV/c]", "Jet mass cut efficiency", 10, 0, 1000, 0, 1, "Jet_mass_cut_eff_CApruned_H_matched_BprimeBprimeToBHBHinc_M-1500.eps", 1., 0.9);

  efficiency1D("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
               "#splitline{QCD, CA R=0.8}{75<m_{jet}<135 GeV/c^{2} (pruned)}",
               "Fat jet p_{T} [GeV/c]", "Jet mass cut efficiency", 10, 0, 1000, 0, 1, "Jet_mass_cut_eff_CApruned_QCDPythia6.eps", 1., 0.9);

//   efficiency1D("output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
//                "W, CA R=0.8, #DeltaR(W,jet)<0.5, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 10, 0, 1000, 0, 1, "Jet_mass_cut_eff_CApruned_W_matched_BprimeBprimeToTWTWinc_M-1300.eps", 1., 0.9);
//
//   efficiency1D("output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
//                "Z, CA R=0.8, #DeltaR(Z,jet)<0.5, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 10, 0, 1000, 0, 1, "Jet_mass_cut_eff_CApruned_Z_matched_BprimeBprimeToBZBZinc_M-1200.eps", 1., 0.9);
//
//   efficiency1D("output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
//                "t, CA R=0.8, #DeltaR(t,jet)<0.5, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 10, 0, 1000, 0, 1, "Jet_mass_cut_eff_CApruned_Top_matched_TprimeToTHinc_M-1700.eps", 1., 0.9);


//   // Jet mass cut efficiency for trimmed jet mass with no PF CHS
//   efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_ExplicitJTA_noPFchs.root",
//                "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
//                "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_noPFchs_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);
//
//   efficiency1D("output_files_v2/QCDPythia6_HiggsTagging_ExplicitJTA_noPFchs.root",
//                "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
//                "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_noPFchs_QCDPythia6.eps", 1., 0.9);
//
//   //efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//   //             "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonDecayProdMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonDecayProdMatched",
//   //             "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//   //             "Jet p_{T} [GeV]", "Jet mass cut efficiency", 40, 0, 1000, 0, 1, "Jet_mass_cut_eff_HiggsToBBbar_BosonDecayProdMatched_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);

  //--------------------------------------------------------------------------------------------------------------------

//   // b-tagging efficiency
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.8, "k_{T}");
//   // b-tagging efficiency with filtered subjets
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_FilteredSubJets_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.8, "Filt");
//   // b-tagging efficiency with pruned subjets
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                        "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (pruned)",
//                        "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_PrunedSubJets_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.8, "Pruned");
//   // b-tagging efficiency with enlarged JTA cone
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                      "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass",
//                      "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                      "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_JTACone_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);
//   // b-tagging efficiency with N-subj. cut applied
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                       "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass_Nsubj", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
//                       "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                       "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_Nsubj_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);
//   // b-tagging efficiency (explicit JTA)
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_ExplicitJTA.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_ExplicitJTA_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.8, "k_{T}");
//   // b-tagging efficiency (explicit JTA, no PF CHS)
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_ExplicitJTA_noPFchs.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_ExplicitJTA_noPFchs_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.8, "k_{T}");
//   // b-tagging efficiency with filtered subjets (explicit JTA)
//   //efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_ExplicitJTA.root",
//   //                     "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass",
//   //                     "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//   //                     "Jet p_{T} [GeV]", "b-tagging efficiency", 40, 0, 1000, 0, 1, "b-tag_eff_FilteredSubJets_ExplicitJTA_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.8, "Filt");
//
//   efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PFJTA.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "H#rightarrowb#bar{b}, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "b-tag eff (subjet CSVL, PF JTA)", 10, 0, 1000, 0, 1, "b-tag_eff_SubjetCSVL_PFJTA_CAPrunedJetMass_HiggsToBBbar_BprimeBprimeToBHBHinc_M-1500.eps", 1., 0.9);
//
//   efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF-CSV.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "H#rightarrowb#bar{b}, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "b-tag eff (subjet CSVL, IVF)", 10, 0, 1000, 0, 1, "b-tag_eff_SubjetCSVL_IVF_CAPrunedJetMass_HiggsToBBbar_BprimeBprimeToBHBHinc_M-1500.eps", 1., 0.9);
//
//   //efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_BJM.root",
//   //             "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//   //             "H#rightarrowb#bar{b}, CA R=0.8, 75<m<135 GeV (pruned)",
//   //             "Jet p_{T} [GeV]", "b-tag eff (subjet CSVL)", 10, 0, 1000, 0, 1, "b-tag_eff_SubjetCSVL_CAPrunedJetMass_HiggsToBBbar_BprimeBprimeToBHBHinc_M-1500_BJM.eps", 1., 0.9);


  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
               "#splitline{H(120)#rightarrowb#bar{b}, CA R=0.8, #DeltaR(H,jet)<0.5}{75<m_{jet}<135 GeV/c^{2} (pruned), Subjet CSVL}",
               "Fat jet p_{T} [GeV/c]", "Double-b-tagging efficiency", 10, 0, 1000, 0, 1, "b-tag_eff_SubjetCSVL_CAPrunedJetMass_HiggsToBBbar_BprimeBprimeToBHBHinc_M-1500.eps", 1., 0.9);

  efficiency1D("output_files_v2/RadionToHHTo4B_M-1500_HiggsTagging_dRsubjetBhadron_JetMass75to140_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
               "#splitline{H(125)#rightarrowb#bar{b}, CA R=0.8, #DeltaR(H,jet)<0.5}{75<m_{jet}<140 GeV/c^{2} (pruned), Subjet CSVL}",
               "Fat jet p_{T} [GeV/c]", "Double-b-tagging efficiency", 10, 0, 1000, 0, 1, "b-tag_eff_SubjetCSVL_CAPrunedJetMass_HiggsToBBbar_RadionToHHTo4B_M-1500.eps", 1., 0.9);

  efficiency1D("output_files_v2/RadionToHHTo4B_M-2000_HiggsTagging_dRsubjetBhadron_JetMass75to140_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
               "#splitline{H(125)#rightarrowb#bar{b}, CA R=0.8, #DeltaR(H,jet)<0.5}{75<m_{jet}<140 GeV/c^{2} (pruned), Subjet CSVL}",
               "Fat jet p_{T} [GeV/c]", "Double-b-tagging efficiency", 10, 0, 1000, 0, 1, "b-tag_eff_SubjetCSVL_CAPrunedJetMass_HiggsToBBbar_RadionToHHTo4B_M-2000.eps", 1., 0.9);

  // overlay multiple jet mass selections
  //efficiency1D_overlayMulti_3("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
  //                          "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_JetMass75to106_CA8only.root",
  //                          "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_JetMass106to135_CA8only.root",
  //                          "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
  //                          "#splitline{H#rightarrowb#bar{b}, CA R=0.8, #DeltaR(H,jet)<0.5}{Subjet CSVL}",
  //                          "Fat jet p_{T} [GeV/c]", "b-tagging efficiency", 10, 0, 1000, 0.001, 1,
  //                          "b-tag_eff_SubjetCSVL_CAPrunedJetMass_JetMassCutComparison.eps", 0, 1., 1.);

  // overlay multiple backgrounds
  efficiency1D_overlayMulti_5("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
                            "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
                            "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
                            "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
                            "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
                            "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
                            "#splitline{CA R=0.8}{75<m_{jet}<135 GeV/c^{2} (pruned), Subjet CSVL}",
                            "Fat jet p_{T} [GeV/c]", "Double-b-tagging efficiency", 10, 0, 1000, 0.001, 1,
                            "b-tag_eff_SubjetCSVL_CAPrunedJetMass.eps", 1, 1., 1.);

  //--------------------------------------------------------------------------------------------------------------------

//   // b-tagging mistag rate
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_QCDPythia6.eps", 1., 1., 0.12, 0.07, 0.8, "k_{T}");
//   // mistag rate when b from gluon splitting
//   efficiency1D_overlay("QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMassbquarks/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassbquarks/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 1, "b-tag_mistag_rate_QCDPythia6_bquarks.eps", 1., 1., 0.12, 0.07, 0.8, "k_{T}");
//   // mistag rate when jet is a c quark jet
//   efficiency1D_overlay("QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMasscquarks/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMasscquarks/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.4, "b-tag_mistag_rate_QCDPythia6_cquarks.eps", 1., 1., 0.12, 0.07, 0.8, "k_{T}");
//   // mistag rate when jet is from u,d,s quark or gluon
//   efficiency1D_overlay("QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMassudsquarks/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassudsquarks/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_QCDPythia6_udsquarks.eps", 1., 1., 0.12, 0.07, 0.8, "k_{T}");
//   // b-tagging mistag rate with filtered subjets
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_FilteredSubJets_QCDPythia6.eps", 1., 1., 0.12, 0.07, 0.8, "Filt");
//   // b-tagging mistag rate with pruned subjets
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (pruned)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_PrunedSubJets_QCDPythia6.eps", 1., 1., 0.12, 0.07, 0.8, "Pruned");
//   // b-tagging mistag rate with enlarged JTA cone
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging.root",
//                       "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass",
//                       "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                       "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_JTACone_QCDPythia6.eps", 1., 1.);
//   // b-tagging mistag rate (explicit JTA)
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging_ExplicitJTA.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_ExplicitJTA_QCDPythia6.eps", 1., 1., 0.12, 0.07, 0.8, "k_{T}");
//   // b-tagging mistag rate (explicit JTA, no PF CHS)
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging_ExplicitJTA_noPFchs.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_ExplicitJTA_noPFchs_QCDPythia6.eps", 1., 1., 0.12, 0.07, 0.8, "k_{T}");
//   // b-tagging mistag rate with filtered subjets (explicit JTA)
//   //efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging_ExplicitJTA.root",
//   //                     "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassFilteredSub/h1_JetPt_BosonMatched_JetMass",
//   //                     "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//   //                     "Jet p_{T} [GeV]", "Mistag rate", 40, 0, 1000, 0, 0.2, "b-tag_mistag_rate_FilteredSubJets_ExplicitJTA_QCDPythia6.eps", 1., 1., 0.12, 0.07, 0.8, "Filt");

  // QCD background
  efficiency1D("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
               "#splitline{QCD, CA R=0.8}{75<m_{jet}<135 GeV/c^{2} (pruned), Subjet CSVL}",
               "Fat jet p_{T} [GeV/c]", "Double-b-misid. probability", 10, 0, 1000, 0, 0.2, "b-tag_mistag_rate_SubjetCSVL_CAPrunedJetMass_QCDPythia6.eps", 1., 1.1, 0.13);

//   // W background
//   efficiency1D("output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, W)", 10, 0, 1000, 0, 0.5, "b-tag_mistag_rate_SubjetMinCSVL_CAPrunedJetMass_WBkg.eps", 1., 1.1, 0.13);
//
//   efficiency1D("output_files_v2/RSGravitonToWW_kMpl01_M-1000_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, W)", 10, 0, 1000, 0, 0.5, "b-tag_mistag_rate_SubjetMinCSVL_CAPrunedJetMass_WBkg_RSGravitonToWW_kMpl01_M-1000.eps", 1., 1.1, 0.13);
//
//   efficiency1D("output_files_v2/RSGravitonToWW_kMpl01_M-1500_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, W)", 10, 0, 1000, 0, 0.5, "b-tag_mistag_rate_SubjetMinCSVL_CAPrunedJetMass_WBkg_RSGravitonToWW_kMpl01_M-1500.eps", 1., 1.1, 0.13);
//
//   efficiency1D("output_files_v2/WWtoAnything_ptmin500_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, W)", 10, 0, 1000, 0, 0.5, "b-tag_mistag_rate_SubjetMinCSVL_CAPrunedJetMass_WBkg_WWtoAnything_ptmin500.eps", 1., 1.1, 0.13);

//     // Z background
//     efficiency1D("output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, Z)", 10, 0, 1000, 0, 0.5, "b-tag_mistag_rate_SubjetMinCSVL_CAPrunedJetMass_ZBkg.eps", 1., 1.1, 0.13);

//     // top background
//     efficiency1D("output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, top)", 10, 0, 1000, 0, 0.5, "b-tag_mistag_rate_SubjetMinCSVL_CAPrunedJetMass_TopBkg.eps", 1., 1.1, 0.13);

//     efficiency1D("output_files_v2/QCDPythia6_HiggsTagging_PFJTA.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, PF JTA)", 10, 0, 1000, 0, 0.2, "b-tag_mistag_rate_SubjetCSVL_PFJTA_CAPrunedJetMass_QCDPythia6.eps", 1., 1.1, 0.13);
//
//     efficiency1D("output_files_v2/QCDPythia6_HiggsTagging_IVF-CSV.root",
//                "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass",
//                "QCD, CA R=0.8, 75<m<135 GeV (pruned)",
//                "Jet p_{T} [GeV]", "Mistag rate (subjet CSVL, IVF)", 10, 0, 1000, 0, 0.2, "b-tag_mistag_rate_SubjetCSVL_IVF_CAPrunedJetMass_QCDPythia6.eps", 1., 1.1, 0.13);

  //--------------------------------------------------------------------------------------------------------------------

//   // Total Higgs tagging efficiency (jet mass cut + b tagging)
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Higgs tagging efficiency", 40, 0, 1000, 0, 1, "Higgs_tag_eff_total_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);
//   // Total Higgs tagging efficiency (jet mass cut + b tagging) with pruned subjets
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                        "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (pruned)",
//                        "Jet p_{T} [GeV]", "Higgs tagging efficiency", 40, 0, 1000, 0, 1, "Higgs_tag_eff_total_PrunedSubJets_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.8, "Pruned");
//   // Total Higgs tagging efficiency (jet mass cut + b tagging) with enlarged JTA cone
//   efficiency1D_overlay("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched",
//                        "H#rightarrowb#bar{b}, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Higgs tagging efficiency", 40, 0, 1000, 0, 1, "Higgs_tag_eff_total_JTACone_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9);


  efficiency1D("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
               "#splitline{H(120)#rightarrowb#bar{b}, CA R=0.8, #DeltaR(H,jet)<0.5}{75<m_{jet}<135 GeV/c^{2} (pruned), Subjet CSVL}",
               "Fat jet p_{T} [GeV/c]", "Higgs-tagging efficiency", 10, 0, 1000, 0, 1, "Higgs_tag_eff_total_CAPrunedJetMass_HiggsToBBbar_BprimeBprimeToBHBHinc_M-1500.eps", 1., 0.9);

  // overlay multiple backgrounds
  efficiency1D_overlayMulti_5("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
                            "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
                            "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
                            "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
                            "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
                            "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
                            "#splitline{CA R=0.8}{75<m_{jet}<135 GeV/c^{2} (pruned), Subjet CSVL}",
                            "Fat jet p_{T} [GeV/c]", "Higgs-tagging efficiency", 10, 0, 1000, 0.001, 1,
                            "Higgs_tag_eff_total_CAPrunedJetMass.eps", 1, 1., 1.);

  //--------------------------------------------------------------------------------------------------------------------

//   // Total Higgs mistag rate (jet mass cut + b tagging)
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMass/h1_JetPt_BosonMatched",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Higgs mistag rate", 40, 0, 1000, 0, 0.05, "Higgs_mistag_rate_total_QCDPythia6.eps", 1., 1.15, 0.14, 0.07, 0.8);
//   // Total Higgs mistag rate (jet mass cut + b tagging) with pruned subjets
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerPrunedJetMass/h1_JetPt_BosonMatched",
//                        "QCD, AK R=0.8, 75<m<135 GeV (pruned)",
//                        "Jet p_{T} [GeV]", "Higgs mistag rate", 40, 0, 1000, 0, 0.05, "Higgs_mistag_rate_total_PrunedSubJets_QCDPythia6.eps", 1., 1.15, 0.14, 0.07, 0.8, "Pruned");
//   // Total Higgs mistag rate (jet mass cut + b tagging) with enlarged JTA cone
//   efficiency1D_overlay("output_files_v2/QCDPythia6_HiggsTagging.root",
//                        "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched_JetMass", "jetAnalyzerTrimmedJetMassJTACone/h1_JetPt_BosonMatched",
//                        "QCD, AK R=0.8, 75<m<135 GeV (trimmed)",
//                        "Jet p_{T} [GeV]", "Higgs mistag rate", 40, 0, 1000, 0, 0.05, "Higgs_mistag_rate_total_JTACone_QCDPythia6.eps", 1., 1.15, 0.14, 0.07, 0.8);

  efficiency1D("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
               "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL", "jetAnalyzerCAPrunedJetMass/h1_JetPt_BosonMatched",
               "#splitline{QCD, CA R=0.8}{75<m_{jet}<135 GeV/c^{2} (pruned), Subjet CSVL}",
               "Fat jet p_{T} [GeV/c]", "Higgs-misid. probability", 10, 0, 1000, 0, 0.02, "Higgs_mistag_rate_total_CAPrunedJetMass_QCDPythia6.eps", 1., 1.2, 0.14);

  //--------------------------------------------------------------------------------------------------------------------

  // QCD background
  // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
  efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
                         "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
                         "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV", "Subjet CSV",
                         0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass.eps", 1);

  efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
                         "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
                         "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV", "Subjet CSV",
                         0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass.eps", 1);

//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_Pt700toInf_JetMass.eps", 1);
//
//   // QCD background, uds jets
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, uds)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_udsJets_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, uds)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_udsJets_Pt700toInf_JetMass.eps", 1);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, uds)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_udsJets_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, uds)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_udsJets_Pt700toInf_JetMass.eps", 1);
//
//   // QCD background, gluon jets
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, g)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_gJets_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, g)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_gJets_Pt700toInf_JetMass.eps", 1);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, g)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_gJets_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, g)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_gJets_Pt700toInf_JetMass.eps", 1);
//
//   // QCD background, charm jets
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, c)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_cJets_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, c)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_cJets_Pt700toInf_JetMass.eps", 1);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, c)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_cJets_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_cQuarks/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, c)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-4, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_cJets_Pt700toInf_JetMass.eps", 1);
//
//   // QCD background, ME b jets
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, ME b)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_MEbJets_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, ME b)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_MEbJets_Pt700toInf_JetMass.eps", 0);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, ME b)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_MEbJets_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, ME b)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_MEbJets_Pt700toInf_JetMass.eps", 0);
//
//   // QCD background, GSP b jets
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, GSP b)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_GSPbJets_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, GSP b)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_QCDBkg_GSPbJets_Pt700toInf_JetMass.eps", 0);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, GSP b)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_GSPbJets_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD, GSP b)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_QCDBkg_GSPbJets_Pt700toInf_JetMass.eps", 0);
//
//   // W background
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (W)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_WBkg_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (W)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_WBkg_Pt700toInf_JetMass.eps", 1);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (W)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_WBkg_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToTWTWinc_M-1300_HiggsTagging_WBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (W)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_WBkg_Pt700toInf_JetMass.eps", 1);
//
//   // Z background
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (Z)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_ZBkg_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (Z)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_ZBkg_Pt700toInf_JetMass.eps", 0);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (Z)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_ZBkg_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBZBZinc_M-1200_HiggsTagging_ZBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (Z)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_ZBkg_Pt700toInf_JetMass.eps", 0);
//
//   // top background
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (top)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_TopBkg_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (top)","Fat jet CSV", "Subjet minCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMinCSVL_TopBkg_Pt700toInf_JetMass.eps", 0);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, single b-tagged subjet)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (top)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_TopBkg_Pt300to500_JetMass.eps", 0);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron_CA8only.root",
//                          "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root", "output_files_v2/TprimeToTHinc_M-1700_HiggsTagging_TopBkg_dRsubjetBhadron_CA8only.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (top)","Fat jet CSV", "Subjet maxCSV",
//                          0, 1, 0, 1, "b-tag_eff_vs_mistag_CA8_SubJetMaxCSVL_TopBkg_Pt700toInf_JetMass.eps", 0);

  //--------------------------------------------------------------------------------------------------------------------

  // b-tagging efficiency vs mistag rate (CA8 pruned jets with enlarged JTA cone)
  efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
                         "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root",
                         "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (JTA #DeltaR<0.8)", "Subjet CSV",
                         0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_JTA_QCDBkg_Pt300to500_JetMass.eps", 1);

  efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
                         "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root",
                         "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMassJTACone/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
                         700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (JTA #DeltaR<0.8)", "Subjet CSV",
                         0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_JTA_QCDBkg_Pt700toInf_JetMass.eps", 1);

//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, PF JTA for fat jets)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PFJTA.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_PFJTA.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Fat jet CSV", "Fat jet CSV (PF JTA)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_fat_jet_PFJTA_BoostedH_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PFJTA.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_PFJTA.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Fat jet CSV", "Fat jet CSV (PF JTA)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_fat_jet_PFJTA_BoostedH_Pt700toInf_JetMass.eps", 1);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, PF JTA for subjets)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PFJTA.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_PFJTA.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Subjet CSV", "Subjet CSV (PF JTA)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_subjet_PFJTA_BoostedH_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_PFJTA.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_PFJTA.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Subjet CSV", "Subjet CSV (PF JTA)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_subjet_PFJTA_BoostedH_Pt700toInf_JetMass.eps", 1);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, IVF CSV for fat jets)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF-CSV.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_IVF-CSV.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Fat jet CSV", "Fat jet CSV (IVF)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_fat_jet_IVF_BoostedH_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF-CSV.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_IVF-CSV.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_JetCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Fat jet CSV", "Fat jet CSV (IVF)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_fat_jet_IVF_BoostedH_Pt700toInf_JetMass.eps", 1);
//
//   // b-tagging efficiency vs mistag rate (CA8 pruned jets, IVF CSV for subjets)
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF-CSV.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_IVF-CSV.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          300, 500, "CA R=0.8, 300<p_{T}<500 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Subjet CSV", "Subjet CSV (IVF)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_subjet_IVF_BoostedH_Pt300to500_JetMass.eps", 1);
//
//   efficiency_curves_comp_xrange("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF-CSV.root",
//                          "output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/QCDPythia6_HiggsTagging_IVF-CSV.root",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//                          700, 1100, "CA R=0.8, p_{T}>700 GeV, 75<m<135 GeV (pruned)", "Tagging efficiency (H#rightarrowb#bar{b})", "Mistag rate (QCD)","Subjet CSV", "Subjet CSV (IVF)",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_subjet_IVF_BoostedH_Pt700toInf_JetMass.eps", 1);

  //--------------------------------------------------------------------------------------------------------------------
}
