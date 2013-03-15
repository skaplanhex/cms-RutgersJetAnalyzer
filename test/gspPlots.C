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


void plot_ProjX(const string& fInputFile, const string& fPlot, const double fPRmin, const double fPRmax,
                const string& fTitle, const string& fTitleX, const int fRebinX, const double fXmin, const double fXmax, const double fYmin, const double fYmax,
                const string& fLeg1, const string& fLeg2, const string& fLeg3, const string& fLeg4, const string& fOutputFile,
                const int fCommonNorm=0, const int fLogy=0, const Double_t fLegendOffsetX=0., const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                const Double_t fLeftMargin=0.1, const Double_t fTopMargin=0.06, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  //gStyle->SetStatX(fLeftMargin+fPlotWidth);
  //gStyle->SetStatY(1.-fTopMargin);
  //gStyle->SetStatH(0.25);
  //gStyle->SetStatW(0.25);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h2_bQuarksGSP = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_bQuarksGSP/" + fPlot).c_str());
  TH2D *h2_bQuarksME = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_bQuarksME/" + fPlot).c_str());
  TH2D *h2_udsQuarks = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_udsQuarks/" + fPlot).c_str());
  TH2D *h2_gluons = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_gluons/" + fPlot).c_str());

  TH1D *h1_bQuarksGSP = h2_bQuarksGSP->ProjectionX("_py1",h2_bQuarksGSP->GetYaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetYaxis()->FindBin(fPRmax));
  TH1D *h1_bQuarksME  = h2_bQuarksME->ProjectionX("_py2",h2_bQuarksGSP->GetYaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetYaxis()->FindBin(fPRmax));
  TH1D *h1_udsQuarks  = h2_udsQuarks->ProjectionX("_py3",h2_bQuarksGSP->GetYaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetYaxis()->FindBin(fPRmax));
  TH1D *h1_gluons     = h2_gluons->ProjectionX("_py4",h2_bQuarksGSP->GetYaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetYaxis()->FindBin(fPRmax));

  h1_bQuarksGSP->Rebin(fRebinX);
  h1_bQuarksME->Rebin(fRebinX);
  h1_udsQuarks->Rebin(fRebinX);
  h1_gluons->Rebin(fRebinX);

  h1_bQuarksGSP->Scale(1./( fCommonNorm ? h1_gluons->Integral() : h1_bQuarksGSP->Integral() ));
  h1_bQuarksME->Scale(1./( fCommonNorm ? h1_gluons->Integral() : h1_bQuarksME->Integral() ));
  h1_udsQuarks->Scale(1./( fCommonNorm ? h1_gluons->Integral() : h1_udsQuarks->Integral() ));
  h1_gluons->Scale(1./h1_gluons->Integral());

  h1_bQuarksGSP->SetLineColor(kRed);
  h1_bQuarksGSP->SetLineWidth(2);
  h1_bQuarksGSP->SetFillColor(kRed);
  h1_bQuarksGSP->SetFillStyle(3004);

  h1_bQuarksME->SetLineColor(kBlue);
  h1_bQuarksME->SetLineWidth(2);
  h1_bQuarksME->SetFillColor(kBlue);
  h1_bQuarksME->SetFillStyle(3005);

  h1_udsQuarks->SetLineColor(kGreen+2);
  h1_udsQuarks->SetLineWidth(2);
  h1_udsQuarks->SetFillColor(kGreen+2);
  h1_udsQuarks->SetFillStyle(3003);

  h1_gluons->SetLineColor(kBlack);
  h1_gluons->SetLineWidth(2);
  h1_gluons->SetFillColor(kBlack);
  h1_gluons->SetFillStyle(3006);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h1_bQuarksGSP->GetXaxis()->SetTitle(fTitleX.c_str());
  h1_bQuarksGSP->SetTitleOffset(fTitleOffsetX,"X");
  h1_bQuarksGSP->SetTitleOffset(fTitleOffsetY,"Y");
  h1_bQuarksGSP->GetXaxis()->SetRangeUser(fXmin,fXmax);
  h1_bQuarksGSP->GetYaxis()->SetRangeUser(fYmin,fYmax);

  h1_bQuarksGSP->Draw("hist");
  h1_bQuarksME->Draw("histsame");
  h1_udsQuarks->Draw("histsame");
  h1_gluons->Draw("histsame");

  TLegend *legend = new TLegend(fLegendOffsetX+.15,.70,fLegendOffsetX+.35,.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(h1_bQuarksGSP, fLeg1.c_str(),"lf");
  legend->AddEntry(h1_bQuarksME, fLeg2.c_str(),"lf");
  legend->AddEntry(h1_udsQuarks, fLeg3.c_str(),"lf");
  legend->AddEntry(h1_gluons, fLeg4.c_str(),"lf");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.1,0.96, fTitle.c_str());

  if(fLogy) c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void plot_ProjY(const string& fInputFile, const string& fPlot, const double fPRmin, const double fPRmax,
                const string& fTitle, const string& fTitleX, const int fRebinX, const double fXmin, const double fXmax, const double fYmin, const double fYmax,
                const string& fLeg1, const string& fLeg2, const string& fLeg3, const string& fLeg4, const string& fOutputFile,
                const int fCommonNorm=0, const int fLogy=0, const Double_t fLegendOffsetX=0., const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0,
                const Double_t fLeftMargin=0.1, const Double_t fTopMargin=0.06, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  //gStyle->SetStatX(fLeftMargin+fPlotWidth);
  //gStyle->SetStatY(1.-fTopMargin);
  //gStyle->SetStatH(0.25);
  //gStyle->SetStatW(0.25);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h2_bQuarksGSP = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_bQuarksGSP/" + fPlot).c_str());
  TH2D *h2_bQuarksME = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_bQuarksME/" + fPlot).c_str());
  TH2D *h2_udsQuarks = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_udsQuarks/" + fPlot).c_str());
  TH2D *h2_gluons = (TH2D*)file->Get(("jetAnalyzerCAPrunedJetMass_gluons/" + fPlot).c_str());

  TH1D *h1_bQuarksGSP = h2_bQuarksGSP->ProjectionY("_py1",h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax));
  TH1D *h1_bQuarksME  = h2_bQuarksME->ProjectionY("_py2",h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax));
  TH1D *h1_udsQuarks  = h2_udsQuarks->ProjectionY("_py3",h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax));
  TH1D *h1_gluons     = h2_gluons->ProjectionY("_py4",h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax));

  h1_bQuarksGSP->Rebin(fRebinX);
  h1_bQuarksME->Rebin(fRebinX);
  h1_udsQuarks->Rebin(fRebinX);
  h1_gluons->Rebin(fRebinX);

  h1_bQuarksGSP->Scale(1./( fCommonNorm ? h1_gluons->Integral() : h1_bQuarksGSP->Integral() ));
  h1_bQuarksME->Scale(1./( fCommonNorm ? h1_gluons->Integral() : h1_bQuarksME->Integral() ));
  h1_udsQuarks->Scale(1./( fCommonNorm ? h1_gluons->Integral() : h1_udsQuarks->Integral() ));
  h1_gluons->Scale(1./h1_gluons->Integral());

  h1_bQuarksGSP->SetLineColor(kRed);
  h1_bQuarksGSP->SetLineWidth(2);
  h1_bQuarksGSP->SetFillColor(kRed);
  h1_bQuarksGSP->SetFillStyle(3004);

  h1_bQuarksME->SetLineColor(kBlue);
  h1_bQuarksME->SetLineWidth(2);
  h1_bQuarksME->SetFillColor(kBlue);
  h1_bQuarksME->SetFillStyle(3005);

  h1_udsQuarks->SetLineColor(kGreen+2);
  h1_udsQuarks->SetLineWidth(2);
  h1_udsQuarks->SetFillColor(kGreen+2);
  h1_udsQuarks->SetFillStyle(3003);

  h1_gluons->SetLineColor(kBlack);
  h1_gluons->SetLineWidth(2);
  h1_gluons->SetFillColor(kBlack);
  h1_gluons->SetFillStyle(3006);

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  h1_bQuarksGSP->GetXaxis()->SetTitle(fTitleX.c_str());
  h1_bQuarksGSP->SetTitleOffset(fTitleOffsetX,"X");
  h1_bQuarksGSP->SetTitleOffset(fTitleOffsetY,"Y");
  h1_bQuarksGSP->GetXaxis()->SetRangeUser(fXmin,fXmax);
  h1_bQuarksGSP->GetYaxis()->SetRangeUser(fYmin,fYmax);

  h1_bQuarksGSP->Draw("hist");
  h1_bQuarksME->Draw("histsame");
  h1_udsQuarks->Draw("histsame");
  h1_gluons->Draw("histsame");

  TLegend *legend = new TLegend(fLegendOffsetX+.15,.70,fLegendOffsetX+.35,.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(h1_bQuarksGSP, fLeg1.c_str(),"lf");
  legend->AddEntry(h1_bQuarksME, fLeg2.c_str(),"lf");
  legend->AddEntry(h1_udsQuarks, fLeg3.c_str(),"lf");
  legend->AddEntry(h1_gluons, fLeg4.c_str(),"lf");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.1,0.96, fTitle.c_str());

  if(fLogy) c->SetLogy();
  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void efficiency_vs_cut(const string& fInputFile, const string& fTitle, const double fPRmin, const double fPRmax,
               const string& fLeg1, const string& fLeg2, const string& fLeg3, const string& fLeg4, const string& fOutputFile,
               const Double_t fLegendOffsetX=0., const Double_t fTitleOffsetX=1.0, const Double_t fTitleOffsetY=1.0, const Double_t fLeftMargin=0.1, const Double_t fTopMargin=0.06, const Double_t fPlotWidth=0.8)
{
  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(fTopMargin);
  gStyle->SetPadBottomMargin(1.-fTopMargin-fPlotWidth);
  gStyle->SetPadLeftMargin(fLeftMargin);
  gStyle->SetPadRightMargin(1.-fLeftMargin-fPlotWidth);
  gROOT->ForceStyle();

  TFile *file = new TFile(fInputFile.c_str());

  TH2D *h2_bQuarksGSP = (TH2D*)file->Get("jetAnalyzerCAPrunedJetMass_bQuarksGSP/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass");
  TH2D *h2_bQuarksME = (TH2D*)file->Get("jetAnalyzerCAPrunedJetMass_bQuarksME/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass");
  TH2D *h2_udsQuarks = (TH2D*)file->Get("jetAnalyzerCAPrunedJetMass_udsQuarks/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass");
  TH2D *h2_gluons = (TH2D*)file->Get("jetAnalyzerCAPrunedJetMass_gluons/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass");

  // denominator counts
  double denom_bQuarksGSP = h2_bQuarksGSP->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),0,101);
  double denom_bQuarksME = h2_bQuarksME->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),0,101);
  double denom_udsQuarks = h2_udsQuarks->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),0,101);
  double denom_gluons = h2_gluons->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),0,101);

  TGraph *g_eff_bQuarksGSP = new TGraph(20);
  g_eff_bQuarksGSP->SetName("g_eff_bQuarksGSP");
  g_eff_bQuarksGSP->SetLineColor(kRed);
  g_eff_bQuarksGSP->SetLineWidth(2);
  g_eff_bQuarksGSP->SetLineStyle(1);

  TGraph *g_eff_bQuarksME = new TGraph(20);
  g_eff_bQuarksME->SetName("g_eff_bQuarksME");
  g_eff_bQuarksME->SetLineColor(kBlue);
  g_eff_bQuarksME->SetLineWidth(2);
  g_eff_bQuarksME->SetLineStyle(2);

  TGraph *g_eff_udsQuarks = new TGraph(20);
  g_eff_udsQuarks->SetName("g_eff_udsQuarks");
  g_eff_udsQuarks->SetLineColor(kGreen+2);
  g_eff_udsQuarks->SetLineWidth(2);
  g_eff_udsQuarks->SetLineStyle(3);

  TGraph *g_eff_gluons = new TGraph(20);
  g_eff_gluons->SetName("g_eff_gluons");
  g_eff_gluons->SetLineColor(kBlack);
  g_eff_gluons->SetLineWidth(2);
  g_eff_gluons->SetLineStyle(4);

  for(int i = 0; i<21; ++i)
  {
    double num_bQuarksGSP = h2_bQuarksGSP->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),101-(i*5),101);
    double num_bQuarksME  = h2_bQuarksME->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),101-(i*5),101);
    double num_udsQuarks  = h2_udsQuarks->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),101-(i*5),101);
    double num_gluons     = h2_gluons->Integral(h2_bQuarksGSP->GetXaxis()->FindBin(fPRmin),h2_bQuarksGSP->GetXaxis()->FindBin(fPRmax),101-(i*5),101);

    g_eff_bQuarksGSP->SetPoint(i,(h2_bQuarksGSP->GetYaxis()->GetBinLowEdge(101-(i*5))),(num_bQuarksGSP/denom_bQuarksGSP));
    g_eff_bQuarksME->SetPoint(i,(h2_bQuarksGSP->GetYaxis()->GetBinLowEdge(101-(i*5))),(num_bQuarksME/denom_bQuarksME));
    g_eff_udsQuarks->SetPoint(i,(h2_bQuarksGSP->GetYaxis()->GetBinLowEdge(101-(i*5))),(num_udsQuarks/denom_udsQuarks));
    g_eff_gluons->SetPoint(i,(h2_bQuarksGSP->GetYaxis()->GetBinLowEdge(101-(i*5))),(num_gluons/denom_gluons));
  }

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg","",100,0.,1.,100,0.,1.);
  bkg->GetXaxis()->SetTitle("CSV discriminator cut");
  bkg->GetYaxis()->SetTitle("Efficiency");
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  g_eff_bQuarksGSP->Draw("Lsame");
  g_eff_bQuarksME->Draw("Lsame");
  g_eff_udsQuarks->Draw("Lsame");
  g_eff_gluons->Draw("Lsame");

  TLegend *legend = new TLegend(fLegendOffsetX+.15,.70,fLegendOffsetX+.35,.90);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(g_eff_bQuarksGSP, fLeg1.c_str(),"l");
  legend->AddEntry(g_eff_bQuarksME, fLeg2.c_str(),"l");
  legend->AddEntry(g_eff_udsQuarks, fLeg3.c_str(),"l");
  legend->AddEntry(g_eff_gluons, fLeg4.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.12,0.96, fTitle.c_str());

  c->SaveAs(fOutputFile.c_str());

  delete c;
  delete file;
}


void makePlots()
{
  plot_ProjX("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_tau2tau1_Pt250to450", -1, 2,
             "CA R=0.8, 250<p_{T}<450 GeV, |#eta|<1.5", "m_{jet} [GeV]", 1, 0., 200., 0., 0.18,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "JetMass_Pt250to450.eps", 0, 0, 0.4);

//   plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_tau2tau1_Pt50toInf", 50, 500,
//              "CA R=0.8, p_{T}>50 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "#tau_{2}/#tau_{1}", 2, 0., 1., 0., 0.11,
//              "b jets from gluon splitting", "b jets", "uds jets", "g jets", "tau2tau1_Pt50toInf_m50toInf.eps");

  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_tau2tau1_Pt250to450", 50, 500,
             "CA R=0.8, 250<p_{T}<450 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "#tau_{2}/#tau_{1}", 2, 0., 1., 0., 0.11,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "tau2tau1_Pt250to450_m50toInf.eps");


//   plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_nTracks_Pt50toInf", 50, 500,
//              "CA R=0.8, p_{T}>50 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "Number of associated tracks", 1, 0., 80., 0., 0.1,
//              "b jets from gluon splitting", "b jets", "uds jets", "g jets", "nTracks_Pt50toInf_m50toInf.eps", 0, 0, 0.4);

  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_nTracks_Pt250to450", 50, 500,
             "CA R=0.8, 250<p_{T}<450 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "Number of associated tracks", 1, 0., 80., 0., 0.1,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "nTracks_Pt250to450_m50toInf.eps", 0, 0, 0.4);


//   plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_nSelectedTracks_Pt50toInf", 50, 500,
//              "CA R=0.8, p_{T}>50 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "Number of selected associated tracks", 1, 0., 40., 0., 0.2,
//              "b jets from gluon splitting", "b jets", "uds jets", "g jets", "nSelectedTracks_Pt50toInf_m50toInf.eps", 0, 0, 0.4);

  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_nSelectedTracks_Pt250to450", 50, 500,
             "CA R=0.8, 250<p_{T}<450 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "Number of selected associated tracks", 1, 0., 40., 0., 0.2,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "nSelectedTracks_Pt250to450_m50toInf.eps", 0, 0, 0.4);


//   plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_nSelectedTracks_Pt50toInf", 50, 500,
//              "CA R=0.8, p_{T}>50 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "Number of selected associated tracks", 1, 0., 40., 0., 0.2,
//              "b jets from gluon splitting", "b jets", "uds jets", "g jets", "nSelectedTracks_Pt50toInf_m50toInf.eps", 0, 0, 0.4);

  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_nSelectedTracks_Pt250to450", 50, 500,
             "CA R=0.8, 250<p_{T}<450 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "Number of selected associated tracks", 1, 0., 40., 0., 0.2,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "nSelectedTracks_Pt250to450_m50toInf.eps", 0, 0, 0.4);


  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_SubJetMinCSVL_Pt50toInf", 0, 500,
             "CA R=0.8, p_{T}>50 GeV, |#eta|<1.5, pruned m_{jet}>0 GeV", "SubJet min CSV Discr", 2, 0., 40., 1e-6, 1,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "SubJetMinCSVL_Pt50toInf_m0toInf.eps", 1, 1, 0.2);

  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_SubJetMinCSVL_Pt250to450", 0, 500,
             "CA R=0.8, 250<p_{T}<450 GeV, |#eta|<1.5, pruned m_{jet}>0 GeV", "SubJet min CSV Discr", 2, 0., 40., 1e-6, 1,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "SubJetMinCSVL_Pt250to450_m0toInf.eps", 1, 1, 0.2);

  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_SubJetMinCSVL_Pt50toInf", 50, 500,
             "CA R=0.8, p_{T}>50 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "SubJet min CSV Discr", 2, 0., 40., 1e-6, 1.,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "SubJetMinCSVL_Pt50toInf_m50toInf.eps", 1, 1, 0.2);

  plot_ProjY("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin50.root", "h2_JetMass_SubJetMinCSVL_Pt250to450", 50, 500,
             "CA R=0.8, 250<p_{T}<450 GeV, |#eta|<1.5, pruned m_{jet}>50 GeV", "SubJet min CSV Discr", 2, 0., 40., 1e-6, 1.,
             "b jets from gluon splitting", "b jets", "uds jets", "g jets", "SubJetMinCSVL_Pt250to450_m50toInf.eps", 1, 1, 0.2);



//   efficiency_vs_cut("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor_JetMass0toInf_JetPtMin0.root",
//                     "CA R=0.8, p_{T}>50 GeV, |#eta|<1.5", 50., 1100.,
//                     "b jets from gluon splitting", "b jets", "uds jets", "g jets", "efficiency_vs_cut_Pt50toInf_test.eps", 0.4, 1, 0.9, 0.12);

//   efficiency_vs_cut("output_files_v2/QCDPythia6_HiggsTagging_dRsubjetBhadron_jetFlavor.root",
//                     "CA R=0.8, p_{T}>300 GeV, |#eta|<1.5", 300., 500.,
//                     "b jets from gluon splitting", "b jets", "uds jets", "g jets", "efficiency_vs_cut_Pt300to500_test.eps", 0.4, 1, 0.9, 0.12);

}