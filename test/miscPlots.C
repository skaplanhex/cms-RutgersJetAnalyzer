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
  gStyle->SetStatH(0.25);
  gStyle->SetStatW(0.25);
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
  gStyle->SetStatH(0.25);
  gStyle->SetStatW(0.25);
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

  TH2D *bkg = new TH2D("bkg", "bkg", 100, 0, 1, 100, 0., 0.08);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  h1_S->Rebin(fRebinX);
  h1_S->SetLineColor(kBlue+2);

  h1_B->Rebin(fRebinX);
  h1_B->SetLineColor(kRed);

  h1_S->Scale(1./h1_S->Integral());
  h1_B->Scale(1./h1_B->Integral());
  
  h1_S->Draw("histsame");
  h1_B->Draw("histsame");

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


void overlay1D_dR(const string& fInputFile1, const string& fInputFile2, const string& fPlot1, const string& fPlot2,
                  const Int_t fPtMin1, const Int_t fPtMax1, const Int_t fPtMin2, const Int_t fPtMax2,
                  const Int_t fRebinX, const Double_t fXmin, const Double_t fXmax, const Double_t fYmin, const Double_t fYmax,
                  const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
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

  TFile *file1 = new TFile(fInputFile1.c_str());
  TFile *file2 = new TFile(fInputFile2.c_str());

  TH2D *h2_1 = (TH2D*)file1->Get(fPlot1.c_str());
  TH1D *h1_1 = h2_1->ProjectionY("_py1",h2_1->GetXaxis()->FindBin(fPtMin1),h2_1->GetXaxis()->FindBin(fPtMax1));
  TH2D *h2_2 = (TH2D*)file2->Get(fPlot2.c_str());
  TH1D *h1_2 = h2_2->ProjectionY("_py2",h2_2->GetXaxis()->FindBin(fPtMin2),h2_2->GetXaxis()->FindBin(fPtMax2));

  TCanvas *c = new TCanvas("c", "",1000,800);
  c->cd();

  TH2D *bkg = new TH2D("bkg", "bkg", 100, fXmin, fXmax, 100, fYmin, fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(fTitleOffsetX,"X");
  bkg->SetTitleOffset(fTitleOffsetY,"Y");
  bkg->Draw();

  h1_1->Rebin(fRebinX);
  h1_1->SetLineWidth(2);
  h1_1->SetLineColor(kBlue+2);

  h1_2->Rebin(fRebinX);
  h1_2->SetLineWidth(2);
  h1_2->SetLineColor(kRed);
  h1_2->SetLineStyle(2);

  h1_1->Scale(1./h1_1->Integral());
  h1_2->Scale(1./h1_2->Integral());

  h1_1->Draw("histsame");
  h1_2->Draw("histsame");

  TLegend *legend = new TLegend(.4,.65,.9,.85);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(h1_1, fLeg1.c_str(),"l");
  legend->AddEntry(h1_2, fLeg2.c_str(),"l");
  legend->Draw();

  TLatex l1;
  l1.SetTextAlign(12);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(fLeftMargin-0.05,0.96, fTitle.c_str());

  c->SaveAs(fOutputFile.c_str());

  delete legend;
  delete c;
  delete file1;
  delete file2;
}


void makePlots()
{
  //--------------------------------------------------------------------------------------------------------------------
  // W tagging
//   overlay1D("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
//             "jetAnalyzerCAPrunedJets/h2_nPV_MassDrop_Pt500toInf", 0, 52, 2,
//             "CA R=0.6 pruned, p_{T}>500 GeV, 60<m<90 GeVV", "#mu=m_{subjet1}/m_{jet}", "Entries", "Boosted W", "QCD",
//             "Mass_drop_BoostedW_QCD_Pt500toInf.eps", 0.9);
// 
//   overlay1D("output_files_v2/WW500_WTagging.root", "output_files_v2/QCDPythia6_WTagging.root",
//             "jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_Pt500toInf", 0, 52, 2,
//             "AK R=0.6, p_{T}>500 GeV, 65<m<95 GeV (trimmed)", "#tau_{2}/#tau_{1}", "Entries", "Boosted W", "QCD",
//             "tau2tau1_BoostedW_QCD_Pt500toInf.eps", 0.9);
// 
//   overlay1D("output_files_v2/WW500_WTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_WTagging_JetSubstructure.root",
//             "jetAnalyzerPrunedJetMass/h2_nPV_tau2tau1_Pt500toInf", 0, 52, 2,
//             "AK R=0.8, p_{T}>500 GeV, 55<m<95 GeV (pruned), onepass_kt", "#tau_{2}/#tau_{1}", "Entries", "Boosted W", "QCD",
//             "tau2tau1_onepassktaxes_BoostedW_QCD_Pt500toInf.eps", 0.9);
// 
//   overlay1D("output_files_v2/WW500_WTagging_JetSubstructure.root", "output_files_v2/QCDPythia6_WTagging_JetSubstructure.root",
//             "jetAnalyzerPrunedJetMassKtAxes/h2_nPV_tau2tau1_Pt500toInf", 0, 52, 2,
//             "AK R=0.8, p_{T}>500 GeV, 55<m<95 GeV (pruned), kt", "#tau_{2}/#tau_{1}", "Entries", "Boosted W", "QCD",
//             "tau2tau1_ktaxes_BoostedW_QCD_Pt500toInf.eps", 0.9);

  //--------------------------------------------------------------------------------------------------------------------
  // Higgs tagging
//   plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_nPV", "",
//          "Primary Vertex Multiplicity", "Events", 1, 0, 50, "nPV_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);
// 
//   plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonPt", "H#rightarrowanything",
//          "Higgs true p_{T} [GeV]", "Entries", 10, 0, 1000, "Pt_Higgs_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);
// 
//   plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonEta", "H#rightarrowanything",
//          "Higgs true #eta", "Entries", 1, -4, 4, "eta_Higgs_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);
// 
//   plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonPt_Isolated", "H#rightarrowanything, #DeltaR(H,other b' decay products)>0.8",
//          "Higgs true p_{T} [GeV]", "Entries", 10, 0, 1000, "Pt_Higgs_Isolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);
// 
//   plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonEta_Isolated", "H#rightarrowanything, #DeltaR(H,other b' decay products)>0.8",
//          "Higgs true #eta", "Entries", 1, -4, 4, "eta_Higgs_Isolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);
// 
//   plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonPt_DecaySel", "H#rightarrowb#bar{b}, #DeltaR(H,other b' decay products)>0.8",
//          "Higgs true p_{T} [GeV]", "Entries", 10, 0, 1000, "Pt_HiggsToBBbar_Isolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);
// 
//   plot1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h1_BosonEta_DecaySel", "H#rightarrowb#bar{b}, #DeltaR(H,other b' decay products)>0.8",
//          "Higgs true #eta", "Entries", 1, -4, 4, "eta_HiggsToBBbar_Isolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 1.3);
// 
// 
//   plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h2_BosonPt_dRdecay", "H#rightarrowb#bar{b}, #DeltaR(H,other b' decay products)>0.8",
//          "Higgs true p_{T} [GeV]", "#DeltaR(b,#bar{b})", 10, 0, 1000, 1, 0, 5, "Pt_HiggsToBBbar_Isolated_dRdecay_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.11, 0.07, 0.77);
// 
//   plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_NonIsolatedBosons_TrimmedJetMass.root", "jetAnalyzerTrimmedJetMass/h2_JetPt_JetPtOverBosonPt", "H#rightarrowb#bar{b} (non-isolated), AK R=0.8, #DeltaR(H,jet)<0.5",
//          "p_{T}^{jet} [GeV]", "p_{T}^{jet}/p_{T}^{boson}", 10, 300, 1000, 2, 0, 2, "JetPt_JetPtOverBosonPt_HiggsToBBbar_NonIsolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.76);
// 
//   plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h2_JetPt_JetPtOverBosonPt", "H#rightarrowb#bar{b} (isolated), AK R=0.8, #DeltaR(H,jet)<0.5",
//          "p_{T}^{jet} [GeV]", "p_{T}^{jet}/p_{T}^{boson}", 5, 300, 1000, 1, 0, 2, "JetPt_JetPtOverBosonPt_HiggsToBBbar_Isolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.76);
// 
//   plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_NonIsolatedBosons_TrimmedJetMass.root", "jetAnalyzerTrimmedJetMass/h2_JetPt_JetPtOverGenJetPt_BosonMatched", "H#rightarrowb#bar{b} (non-isolated), AK R=0.8, #DeltaR(H,jet)<0.5",
//          "p_{T}^{jet} [GeV]", "p_{T}^{jet}/p_{T}^{genjet}", 10, 300, 1000, 2, 0, 2, "JetPt_JetPtOverGenJetPt_HiggsToBBbar_NonIsolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.76);
// 
//   plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h2_JetPt_JetPtOverGenJetPt_BosonMatched", "H#rightarrowb#bar{b} (isolated), AK R=0.8, #DeltaR(H,jet)<0.5",
//          "p_{T}^{jet} [GeV]", "p_{T}^{jet}/p_{T}^{genjet}", 5, 300, 1000, 1, 0, 2, "JetPt_JetPtOverGenJetPt_HiggsToBBbar_Isolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.76);
// 
//   plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging_NonIsolatedBosons_TrimmedJetMass.root", "jetAnalyzerTrimmedJetMass/h2_JetPt_JetMass_BosonMatched", "H#rightarrowb#bar{b} (non-isolated), AK R=0.8, #DeltaR(H,jet)<0.5",
//          "p_{T}^{jet} [GeV]", "m_{jet} (trimmed) [GeV]", 10, 300, 1000, 4, 0, 400, "JetPt_JetMass_HiggsToBBbar_NonIsolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.76);
// 
//   plot2D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "jetAnalyzerTrimmedJetMass/h2_JetPt_JetMass_BosonMatched", "H#rightarrowb#bar{b} (isolated), AK R=0.8, #DeltaR(H,jet)<0.5",
//          "p_{T}^{jet} [GeV]", "m_{jet} (trimmed) [GeV]", 5, 300, 1000, 2, 0, 400, "JetPt_JetMass_HiggsToBBbar_Isolated_BprimeBprimeToBHBHinc_M-800.eps", 1., 0.9, 0.12, 0.07, 0.76);
// 
// 
//   overlay1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
//             "jetAnalyzerCAPrunedJets/h2_nPV_MassDrop_Pt300toInf", 0, 52, 2,
//             "CA R=0.8 pruned, p_{T}>300 GeV, 75<m<135 GeV", "#mu=m_{subjet1}/m_{jet}", "Entries", "Boosted H", "QCD",
//             "Mass_drop_BoostedH_QCD_Pt300toInf.eps", 0.9);
// 
//   overlay1D("output_files_v2/BprimeBprimeToBHBHinc_M-800_HiggsTagging.root", "output_files_v2/QCDPythia6_HiggsTagging.root",
//             "jetAnalyzerTrimmedJetMass/h2_nPV_tau2tau1_Pt300toInf", 0, 52, 2,
//             "AK R=0.8, p_{T}>300 GeV, 75<m<135 GeV (trimmed)", "#tau_{2}/#tau_{1}", "Entries", "Boosted H", "QCD",
//             "tau2tau1_BoostedH_QCD_Pt300toInf.eps", 0.9);

  // dR(subjet,Bhadron)
  // Pruned subjets
  // 300<pT<500 GeV
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron", "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass",
//                300, 500, 300, 500, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, pruned subjets", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "No m cut", "75<m<135 GeV (pruned)",
//                "mindRSubjet1Bhadron_AK8pruned_BoostedH_Pt300to500_JetMassCut.eps", 0.95);
// 
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron", "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass",
//                300, 500, 300, 500, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, pruned subjets", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "No m cut", "75<m<135 GeV (pruned)",
//                "mindRSubjet2Bhadron_AK8pruned_BoostedH_Pt300to500_JetMassCut.eps", 0.95);
// 
//   // 500<pT<700 GeV
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron", "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass",
//                500, 700, 500, 700, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, 500<p_{T}<700 GeV, #DeltaR(H,jet)<0.5, pruned subjets", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "No m cut", "75<m<135 GeV (pruned)",
//                "mindRSubjet1Bhadron_AK8pruned_BoostedH_Pt500to700_JetMassCut.eps", 0.95);
// 
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron", "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass",
//                500, 700, 500, 700, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, 500<p_{T}<700 GeV, #DeltaR(H,jet)<0.5, pruned subjets", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "No m cut", "75<m<135 GeV (pruned)",
//                "mindRSubjet2Bhadron_AK8pruned_BoostedH_Pt500to700_JetMassCut.eps", 0.95);
// 
//   // pT>700 GeV
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron", "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass",
//                700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, pruned subjets", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "No m cut", "75<m<135 GeV (pruned)",
//                "mindRSubjet1Bhadron_AK8pruned_BoostedH_Pt700toInf_JetMassCut.eps", 0.95);
// 
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron", "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass",
//                700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, pruned subjets", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "No m cut", "75<m<135 GeV (pruned)",
//                "mindRSubjet2Bhadron_AK8pruned_BoostedH_Pt700toInf_JetMassCut.eps", 0.95);

  // Pruned vs kT subjets
  // 300<pT<500 GeV
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass", "jetAnalyzerPrunedJetMassKtSub/h2_JetPt_mindRSubjet1Bhadron_JetMass",
//                300, 500, 300, 500, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "Pruned subjets", "k_{T} subjets",
//                "mindRSubjet1Bhadron_AK8_BoostedH_Pt300to500_JetMass_pruned_vs_kT.eps", 0.95);
// 
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass", "jetAnalyzerPrunedJetMassKtSub/h2_JetPt_mindRSubjet2Bhadron_JetMass",
//                300, 500, 300, 500, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "Pruned subjets", "k_{T} subjets",
//                "mindRSubjet2Bhadron_AK8_BoostedH_Pt300to500_JetMass_pruned_vs_kT.eps", 0.95);
// 
//   // pT>700 GeV
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass", "jetAnalyzerPrunedJetMassKtSub/h2_JetPt_mindRSubjet1Bhadron_JetMass",
//                700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "Pruned subjets", "k_{T} subjets",
//                "mindRSubjet1Bhadron_AK8_BoostedH_Pt700toInf_JetMass_pruned_vs_kT.eps", 0.95);
// 
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass", "jetAnalyzerPrunedJetMassKtSub/h2_JetPt_mindRSubjet2Bhadron_JetMass",
//                700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "Pruned subjets", "k_{T} subjets",
//                "mindRSubjet2Bhadron_AK8_BoostedH_Pt700toInf_JetMass_pruned_vs_kT.eps", 0.95);

  // Pruned vs filtered subjets
  // 300<pT<500 GeV
  overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
               "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass", "jetAnalyzerPrunedJetMassFilteredSub/h2_JetPt_mindRSubjet1Bhadron_JetMass",
               300, 500, 300, 500, 1, 0, 1, 0, 0.8,
               "AK R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "Pruned subjets", "Filtered subjets",
               "mindRSubjet1Bhadron_AK8_BoostedH_Pt300to500_JetMass_pruned_vs_filtered.eps", 0.95);

  overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
               "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass", "jetAnalyzerPrunedJetMassFilteredSub/h2_JetPt_mindRSubjet2Bhadron_JetMass",
               300, 500, 300, 500, 1, 0, 1, 0, 0.8,
               "AK R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "Pruned subjets", "Filtered subjets",
               "mindRSubjet2Bhadron_AK8_BoostedH_Pt300to500_JetMass_pruned_vs_filtered.eps", 0.95);

  // pT>700 GeV
  overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
               "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass", "jetAnalyzerPrunedJetMassFilteredSub/h2_JetPt_mindRSubjet1Bhadron_JetMass",
               700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
               "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "Pruned subjets", "Filtered subjets",
               "mindRSubjet1Bhadron_AK8_BoostedH_Pt700toInf_JetMass_pruned_vs_filtered.eps", 0.95);

  overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
               "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass", "jetAnalyzerPrunedJetMassFilteredSub/h2_JetPt_mindRSubjet2Bhadron_JetMass",
               700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
               "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "Pruned subjets", "Filtered subjets",
               "mindRSubjet2Bhadron_AK8_BoostedH_Pt700toInf_JetMass_pruned_vs_filtered.eps", 0.95);

  // Pruned subjets AK vs CA
  // 300<pT<500 GeV
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass",
//                300, 500, 300, 500, 1, 0, 1, 0, 0.8,
//                "R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "AK pruned subjets", "CA pruned subjets",
//                "mindRSubjet1Bhadron_BoostedH_Pt300to500_JetMass_pruned_AK8_vs_CA8.eps", 0.95);
// 
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1000_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass",
//                300, 500, 300, 500, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, 300<p_{T}<500 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "AK pruned subjets", "CA pruned subjets",
//                "mindRSubjet2Bhadron_BoostedH_Pt300to500_JetMass_pruned_AK8_vs_CA8.eps", 0.95);
// 
//   // pT>700 GeV
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_mindRSubjet1Bhadron_JetMass",
//                700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{1},B hadron)", "Relative fraction", "AK pruned subjets", "CA pruned subjets",
//                "mindRSubjet1Bhadron_BoostedH_Pt700toInf_JetMass_pruned_AK8_vs_CA8.eps", 0.95);
// 
//   overlay1D_dR("output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root", "output_files_v2/BprimeBprimeToBHBHinc_M-1500_HiggsTagging_dRsubjetBhadron.root",
//                "jetAnalyzerPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass", "jetAnalyzerCAPrunedJetMass/h2_JetPt_mindRSubjet2Bhadron_JetMass",
//                700, 1100, 700, 1100, 1, 0, 1, 0, 0.8,
//                "AK R=0.8, p_{T}>700 GeV, #DeltaR(H,jet)<0.5, 75<m<135 GeV (pruned)", "min #DeltaR(subjet_{2},B hadron)", "Relative fraction", "AK pruned subjets", "CA pruned subjets",
//                "mindRSubjet2Bhadron_BoostedH_Pt700toInf_JetMass_pruned_AK8_vs_CA8.eps", 0.95);
  
}