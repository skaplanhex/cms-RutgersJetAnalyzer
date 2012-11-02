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
  gROOT->ForceStyle();

  TFile *file = new TFile((fInputDir + "/"+ fFile).c_str());

  TH2D *h2_nPV_JetMass = (TH2D*)file->Get((fFileDir + "/h2_nPV_JetMass_" + fPtRange).c_str());
  TH1D *h1_JetMass = h2_nPV_JetMass->ProjectionY("_py",fPVLow,fPVHigh);
  h1_JetMass->GetXaxis()->SetRangeUser(fXmin,fXmax);
  h1_JetMass->SetTitleOffset(1.0,"X");

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();

  h1_JetMass->Draw("hists");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14,0.96, fTitle.c_str());

  c->SaveAs(fOutputFile.c_str());
}


void makePlots()
{
  // nPV<=10
  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerDefaultJetMass",
	   "Pt300toInf", 0, 11, "Default jet mass, p_{T}>300 GeV, nPV#leq10",
	   0, 299.5, "Higgs_matched_default_jet_mass_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerFilteredJetMass",
	   "Pt300toInf", 0, 11, "Filtered jet mass, p_{T}>300 GeV, nPV#leq10",
	   0, 299.5, "Higgs_matched_filtered_jet_mass_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerPrunedJetMass",
	   "Pt300toInf", 0, 11, "Pruned jet mass, p_{T}>300 GeV, nPV#leq10",
	   0, 299.5, "Higgs_matched_pruned_jet_mass_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass",
	   "Pt300toInf", 0, 11, "Trimmed jet mass, p_{T}>300 GeV, nPV#leq10",
	   0, 299.5, "Higgs_matched_trimmed_jet_mass_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerCAPrunedJets",
           "Pt300toInf", 0, 11, "CA pruned jet mass, p_{T}>300 GeV, nPV#leq10",
           0, 299.5, "Higgs_matched_CApruned_jet_mass_Pt300toInf_nPV0to10_BprimeBprimeToBHBHinc_M-800.eps");

  // nPV inclusive
  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerDefaultJetMass",
	   "Pt300toInf", 0, 52, "Default jet mass, p_{T}>300 GeV",
	   0, 299.5, "Higgs_matched_default_jet_mass_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerFilteredJetMass",
	   "Pt300toInf", 0, 52, "Filtered jet mass, p_{T}>300 GeV",
	   0, 299.5, "Higgs_matched_filtered_jet_mass_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerPrunedJetMass",
	   "Pt300toInf", 0, 52, "Pruned jet mass, p_{T}>300 GeV",
	   0, 299.5, "Higgs_matched_pruned_jet_mass_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass",
	   "Pt300toInf", 0, 52, "Trimmed jet mass, p_{T}>300 GeV",
	   0, 299.5, "Higgs_matched_trimmed_jet_mass_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");

  jet_mass("output_files_v2", "BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerCAPrunedJets",
           "Pt300toInf", 0, 52, "CA pruned jet mass, p_{T}>300 GeV",
           0, 299.5, "Higgs_matched_CApruned_jet_mass_Pt300toInf_BprimeBprimeToBHBHinc_M-800.eps");
}