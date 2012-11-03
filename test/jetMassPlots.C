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
#include "exoStyle.C"


using namespace std;


void jet_mass(const string& fInputDir, const string& fFile, const string& fFileDir,
              const string& fPtRange, const int fPVLow, const int fPVHigh, const string& fTitle,
              const double fXmin, const double fXmax, const string& fOutputFile)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat("nemruoi");
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetStatX(0.94);
  gStyle->SetStatY(0.93);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  //gStyle->UseCurrentStyle();
  gROOT->ForceStyle();

  TFile *file = new TFile((fInputDir + "/"+ fFile).c_str());

  TH2D *h2_nPV_JetMass = (TH2D*)file->Get((fFileDir + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH1D *h1_JetMass = h2_nPV_JetMass->ProjectionY("_py",fPVLow,fPVHigh);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h1_JetMass->Rebin(2);
  h1_JetMass->GetYaxis()->SetTitle("Entries");
  h1_JetMass->GetXaxis()->SetRangeUser(fXmin,fXmax);
  h1_JetMass->SetTitleOffset(1.0,"X");
  h1_JetMass->SetTitleOffset(1.2,"Y");
  h1_JetMass->SetLineColor(kBlue+2);

  h1_JetMass->Draw("hists");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void makePlots()
{
  // nPV<=10
  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerDefaultJetMass",
           "Pt300toInf", 0, 11, "H#rightarrowb#bar{b}, Default, p_{T}>300 GeV, #DeltaR(H,jet)<0.5, nPV#leq10",
           0, 299.5, "Jet_mass_default_H_matched_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerFilteredJetMass",
           "Pt300toInf", 0, 11, "H#rightarrowb#bar{b}, Filtered, p_{T}>300 GeV, #DeltaR(H,jet)<0.5, nPV#leq10",
           0, 299.5, "Jet_mass_filtered_H_matched_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerPrunedJetMass",
           "Pt300toInf", 0, 11, "H#rightarrowb#bar{b}, Pruned, p_{T}>300 GeV, #DeltaR(H,jet)<0.5, nPV#leq10",
           0, 299.5, "Jet_mass_pruned_H_matched_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass",
           "Pt300toInf", 0, 11, "H#rightarrowb#bar{b}, Trimmed, p_{T}>300 GeV, #DeltaR(H,jet)<0.5, nPV#leq10",
           0, 299.5, "Jet_mass_trimmed_H_matched_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerCAPrunedJets",
           "Pt300toInf", 0, 11, "H#rightarrowb#bar{b}, CA pruned, p_{T}>300 GeV, #DeltaR(H,jet)<0.5, nPV#leq10",
           0, 299.5, "Jet_mass_CApruned_H_matched_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  // nPV inclusive
  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerDefaultJetMass",
           "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, Default, p_{T}>300 GeV, #DeltaR(H,jet)<0.5",
           0, 299.5, "Jet_mass_default_H_matched_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerFilteredJetMass",
           "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, Filtered, p_{T}>300 GeV, #DeltaR(H,jet)<0.5",
           0, 299.5, "Jet_mass_filtered_H_matched_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerPrunedJetMass",
           "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, Pruned, p_{T}>300 GeV, #DeltaR(H,jet)<0.5",
           0, 299.5, "Jet_mass_pruned_H_matched_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass",
           "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, Trimmed, p_{T}>300 GeV, #DeltaR(H,jet)<0.5",
           0, 299.5, "Jet_mass_trimmed_H_matched_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerCAPrunedJets",
           "Pt300toInf", 0, 52, "H#rightarrowb#bar{b}, CA pruned, p_{T}>300 GeV, #DeltaR(H,jet)<0.5",
           0, 299.5, "Jet_mass_CApruned_H_matched_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");
}