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


void plot1D(const string& fInputFile, const string& fPlot, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
            const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const string& fOutputFile,
            const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0, const Double_t fLeftMargin=0.16, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat("nemruoi");
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gStyle->SetStatX(fLeftMargin+fPlotWidth);
  gStyle->SetStatY(1.-fTopMargin);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH1D *h1_plot = (TH1D*)file->Get(fPlot.c_str());

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h1_plot->Rebin(fRebinX);
  h1_plot->GetXaxis()->SetRangeUser(fXmin,fXmax);
  h1_plot->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_plot->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_plot->SetTitleOffset(fTitleOffsetX,"X");
  h1_plot->SetTitleOffset(fTitleOffsetY,"Y");
  h1_plot->SetLineColor(kBlue+2);

  h1_plot->Draw("hists");

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetTextSize(0.05);
  l1.SetNDC();
  l1.DrawLatex(fLeftMargin,0.97, fTitle.c_str());

  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void plot2D(const string& fInputFile, const string& fPlot, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
            const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Int_t fRebinY, const Double_t fYmin, const Double_t fYmax,
            const string& fOutputFile, const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
            const Double_t fLeftMargin=0.16, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetPalette(1);
  gStyle->SetOptStat("nemruoi");
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gStyle->SetStatX(fLeftMargin+fPlotWidth);
  gStyle->SetStatY(1.-fTopMargin);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h2_plot = (TH2D*)file->Get(fPlot.c_str());

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h2_plot->Rebin2D(fRebinX,fRebinY);
  h2_plot->GetXaxis()->SetRangeUser(fXmin,fXmax);
  h2_plot->GetYaxis()->SetRangeUser(fYmin,fYmax);
  h2_plot->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h2_plot->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h2_plot->SetTitleOffset(fTitleOffsetX,"X");
  h2_plot->SetTitleOffset(fTitleOffsetY,"Y");
  h2_plot->SetLineColor(kBlue+2);

  h2_plot->Draw("colz");

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

  delete c;
  delete file;
}

void overlay1D(const string& fInputFileS, const string& fInputFileB, const string& fPlot,
               const Int_t fPVLow, const Int_t fPVHigh, const Int_t fRebinX, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
               const string& fLeg1, const string& fLeg2, const string& fOutputFile,
               const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
               const Double_t fLeftMargin=0.12, const Double_t fTopMargin=0.07, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  //gStyle->UseCurrentStyle();
  gROOT->ForceStyle();

  TFile *file_S = new TFile(fInputFileS.c_str());
  TFile *file_B = new TFile(fInputFileB.c_str());

  TH2D *h2_S = (TH2D*)file_S->Get(fPlot.c_str());
  TH1D *h1_S = h2_S->ProjectionY("_py1",fPVLow,fPVHigh);
  TH2D *h2_B = (TH2D*)file_B->Get(fPlot.c_str());
  TH1D *h1_B = h2_B->ProjectionY("_py2",fPVLow,fPVHigh);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h1_S->Rebin(fRebinX);
  h1_S->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  h1_S->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  h1_S->SetTitleOffset(fTitleOffsetX,"X");
  h1_S->SetTitleOffset(fTitleOffsetY,"Y");
  h1_S->SetLineColor(kBlue+2);

  h1_B->Rebin(fRebinX);
  h1_B->SetLineColor(kRed);

  h1_S->DrawNormalized("hist");
  h1_B->DrawNormalized("histsame");

  TLegend *legend = new TLegend(.65,.65,.95,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_S, fLeg1.c_str(),"l");
  legend->AddEntry(h1_B, fLeg2.c_str(),"l");
  legend->Draw();
  
  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(fLeftMargin,0.96, fTitle.c_str());

  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete c;
  delete file_S;
  delete file_B;
}


void makePlots()
{
  //--------------------------------------------------------------------------------------------------------------------
  // W tagging
  overlay1D("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
            "jetAnalyzerCAPrunedJets/h2_nPV_MassDrop_Pt500toInf", 0, 52, 2,
            "CA R=0.6 pruned, 60<m<90 GeV, p_{T}>500 GeV", "#mu=m_{subjet1}/m_{jet}", "Entries", "Boosted W", "QCD",
            "Mass_drop_boostedW_QCD_Pt500toInf.eps", 0.9);

  overlay1D("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
            "jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_Pt500toInf", 0, 52, 2,
            "Anti-k_{T} R=0.6, 65<m<95 GeV (trimmed mass), p_{T}>500 GeV", "#tau_{2}/#tau_{1}", "Entries", "Boosted W", "QCD",
            "tau2tau1_boostedW_QCD_Pt500toInf.eps", 0.9);

  //--------------------------------------------------------------------------------------------------------------------
  // Higgs tagging
  plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_nPV", "",
         "Primary Vertex Multiplicity", "Events", 1, 0, 50, "nPV_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);

  plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonPt", "H#rightarrowanything",
         "Higgs true p_{T} [GeV]", "Entries", 10, 0, 1000, "Pt_Higgs_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);

  plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonEta", "H#rightarrowanything",
         "Higgs true #eta", "Entries", 1, -4, 4, "eta_Higgs_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);

  plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonPt_DecaySel", "H#rightarrowb#bar{b}",
         "Higgs true p_{T} [GeV]", "Entries", 10, 0, 1000, "Pt_HiggsToBBbar_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);


  plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h2_BosonPt_dRdecay", "H#rightarrowb#bar{b}",
         "Higgs true p_{T} [GeV]", "#DeltaR(b,#bar{b})", 10, 0, 1000, 1, 0, 5, "Pt_HiggsToBBbar_dRdecay_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.11, 0.07, 0.77);


  overlay1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
            "jetAnalyzerCAPrunedJets/h2_nPV_MassDrop_Pt300toInf", 0, 52, 2,
            "CA R=0.8 pruned, 75<m<130 GeV, p_{T}>300 GeV", "#mu=m_{subjet1}/m_{jet}", "Entries", "Boosted H", "QCD",
            "Mass_drop_boostedH_QCD_Pt300toInf.eps", 0.9);

  overlay1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
            "jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_Pt300toInf", 0, 52, 2,
            "Anti-k_{T} R=0.8, 80<m<130 GeV (trimmed mass), p_{T}>300 GeV", "#tau_{2}/#tau_{1}", "Entries", "Boosted H", "QCD",
            "tau2tau1_boostedH_QCD_Pt300toInf.eps", 0.9);
}