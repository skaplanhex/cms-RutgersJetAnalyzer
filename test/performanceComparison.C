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
#include <map>
#include "Rtypes.h"

//CINT sucks
#ifdef __MAKECINT__
#pragma link C++ class std::map<std::string,TGraph*>+;
#pragma link C++ operators std::map<std::string,TGraph*>::const_iterator;
#endif

using namespace std;

TGraph* getEfficiencyCurve(const string& fFileS1, const string& fFileB1,const string& fPlot,const double fXMin, const double fXMax){
  //get files and histograms
  TFile *file_S1  = new TFile(fFileS1.c_str());
  TFile *file_B1 = new TFile(fFileB1.c_str());

  TH2D *h2_S_1 = (TH2D*)file_S1->Get(fPlot.c_str());
  TH2D *h2_B_1 = (TH2D*)file_B1->Get(fPlot.c_str());

  //total jet count for denominator of efficiency calculation
  double denom_S_1 = h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),0,101);
  double denom_B_1 = h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),0,101);

  TGraph *g_eff_1 = new TGraph(29);

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
  return g_eff_1;
}

void formatGraph(TGraph* graph, int graphNum){
  short colors[7]={ kGreen+2, kRed, kBlue, kBlack, kMagenta, kOrange+2, kCyan };
  short graphColor = colors[graphNum % 7];
  int lineStyle = (graphNum % 11) + 1;
  graph->SetLineColor(graphColor);
  graph->SetLineStyle(lineStyle);
  graph->SetLineWidth(2);
}
std::string getHistName(std::string grooming = "Pruned",std::string algo = "CSV", bool subjetPlot = false){
  if (!subjetPlot)
    return "jetAnalyzerCA8FatJets_"+grooming+"Subjets/h2_JetPt_Jet"+algo+"_BosonMatched_JetMass";
  else
    return "jetAnalyzerCA8FatJets_"+grooming+"Subjets/h2_JetPt_SubJetMin"+algo+"_BosonMatched_JetMass"; 
}
void plotEfficiencyCurves(std::map< std::string,TGraph* > &graphs, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle, const string& fExtraInfo, const string& fOutputFile, const double fXmin, const double fXmax, const double fYmin, const double fYmax,const Int_t fLogy=0){

  gROOT->SetBatch(kTRUE);
  setEXOStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gROOT->ForceStyle();

  TCanvas *c = new TCanvas("c", "",800,800);
  c->cd();
  TH2D *bkg = new TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax);
  bkg->GetXaxis()->SetTitle(fXAxisTitle.c_str());
  bkg->GetYaxis()->SetTitle(fYAxisTitle.c_str());
  bkg->SetTitleOffset(1.1,"Y");
  bkg->Draw();
  c->SetGridx();
  c->SetGridy();

  TLegend *legend = new TLegend(.16,.64,.36,.77);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.021);

  int graphCounter = 0;
  for (std::map<std::string,TGraph*>::const_iterator it = graphs.begin(); it != graphs.end(); ++it){
    std::string label = it->first;
    TGraph* graph = it->second;
    legend->AddEntry(graph, label.c_str(),"l");
    formatGraph(graph,graphCounter);
    graph->Draw("L");
    graphCounter++;
  }

  if (fLogy) c->SetLogy();
  legend->Draw();
  TLatex l1;
  l1.SetTextAlign(13);
  l1.SetTextFont(42);
  l1.SetNDC();
  l1.SetTextSize(0.04);
  l1.DrawLatex(0.14+0.03,0.90, fTitle.c_str());

  l1.SetTextAlign(12);
  l1.SetTextSize(0.045);
  l1.SetTextFont(62);
  l1.DrawLatex(0.14,0.96, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");

  l1.SetTextFont(42);
  l1.SetTextSize(0.025);
  l1.DrawLatex(0.48,0.18, fExtraInfo.c_str());
  c->SaveAs(fOutputFile.c_str());
  graphs.clear();
  delete c;
  delete legend;
  delete bkg;
}

void makePlots()
{
  //for multiple plots on the same canvas

  //maps to hold legend entries and TGraph*s
  std::map< std::string,TGraph* > graphsPt300To500;
  std::map< std::string,TGraph* > graphsPt700ToInf;

  // no IVF, no Explicit JTA, no SV Clustering

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Pruned","JP",true),300,500);

  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Pruned","JP",true),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf.png",0, 1, 1E-3, 1,1);

  // IVF, Explicit JTA

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JBP",true),300,500);

  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JBP",true),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF}{Explicit JTA}","btagperfcomp_Pt300to500_IVF_ExplicitJTA.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF}{Explicit JTA}","btagperfcomp_Pt700toInf_IVF_ExplicitJTA.png",0, 1, 1E-3, 1,1);

  // IVF, Explicit JTA, SV Clustering

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JBP",true),300,500);

  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JBP",true),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SV Clustering}{Explicit JTA}","btagperfcomp_Pt300to500_IVF_ExplicitJTA_SVClustering.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SV Clustering}{Explicit JTA}","btagperfcomp_Pt700toInf_IVF_ExplicitJTA_SVClustering.png",0, 1, 1E-3, 1,1);

  // IVF, Explicit JTA, SV Clustering (with SV direction taken along its momentum)

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Pruned","JBP",true),300,500);

  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Pruned","JP",true),700,1100);
  graphsPt700ToInf["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering_SVMomentum.root",getHistName("Pruned","JBP",true),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SV Clustering By SV Momentum}{Explicit JTA}","btagperfcomp_Pt300to500_IVF_ExplicitJTA_SVClustering_SVMomentum.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SV Clustering By SV Momentum}{Explicit JTA}","btagperfcomp_Pt700toInf_IVF_ExplicitJTA_SVClustering_SVMomentum.png",0, 1, 1E-3, 1,1);

  // IVF

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Pruned","JP",true),300,500);

  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Pruned","JP",true),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF}{}","btagperfcomp_Pt300to500_IVF.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF}{}","btagperfcomp_Pt700toInf_IVF.png",0, 1, 1E-3, 1,1);

  // IVF, SV Clustering

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Pruned","JP",true),300,500);

  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging_IVF_SVClustering.root","QCDPythia6_HiggsTagging_IVF_SVClustering.root",getHistName("Pruned","JP",true),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SVClustering}{}","btagperfcomp_Pt300to500_IVF_SVClustering.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SVClustering}{}","btagperfcomp_Pt700toInf_IVF_SVClustering.png",0, 1, 1E-3, 1,1);

  // IVF, Explicit JTA, PFchs

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JBP",true),300,500);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, PFchs}{Explicit JTA}","btagperfcomp_Pt300to500_IVF_ExplicitJTA_PATTuple_v3.png",0, 1, 1E-3, 1,1);

  // IVF, Explicit JTA, SV Clustering, PFchs

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JBP",true),300,500);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SV Clustering, PFchs}{Explicit JTA}","btagperfcomp_Pt300to500_IVF_ExplicitJTA_SVClustering_PATTuple_v3.png",0, 1, 1E-3, 1,1);

  // IVF, PFchs

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_PATTuple_v3.root","QCDPythia6_HiggsTagging_IVF.root",getHistName("Pruned","JP",true),300,500);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, PFchs}{}","btagperfcomp_Pt300to500_IVF_PATTuple_v3.png",0, 1, 1E-3, 1,1);

  // PFchs

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","QCDPythia6_MDPlots.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","QCDPythia6_MDPlots.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","QCDPythia6_MDPlots.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","QCDPythia6_MDPlots.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","QCDPythia6_MDPlots.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3.root","QCDPythia6_MDPlots.root",getHistName("Pruned","JP",true),300,500);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","PFchs","btagperfcomp_Pt300to500_PATTuple_v3.png",0, 1, 1E-3, 1,1);

  // IVF, Explicit JTA, Extended PFchs

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JBP",true),300,500);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, Extended PFchs}{Explicit JTA}","btagperfcomp_Pt300to500_IVF_ExplicitJTA_PATTuple_v3_ExtPFchs.png",0, 1, 1E-3, 1,1);

  // IVF, Explicit JTA, SV Clustering, Extended PFchs

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JP",true),300,500);
  graphsPt300To500["Subjet JBP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.root","QCDPythia6_HiggsTagging_IVF_ExplicitJTA_SVClustering.root",getHistName("Pruned","JBP",true),300,500);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","#splitline{IVF, SV Clustering, Extended PFchs}{Explicit JTA}","btagperfcomp_Pt300to500_IVF_ExplicitJTA_SVClustering_PATTuple_v3_ExtPFchs.png",0, 1, 1E-3, 1,1);

  // Extended PFchs

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3_ExtPFchs.root","QCDPythia6_MDPlots.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3_ExtPFchs.root","QCDPythia6_MDPlots.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3_ExtPFchs.root","QCDPythia6_MDPlots.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3_ExtPFchs.root","QCDPythia6_MDPlots.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3_ExtPFchs.root","QCDPythia6_MDPlots.root",getHistName("KtBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet JP (Pruned)"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging_PATTuple_v3_ExtPFchs.root","QCDPythia6_MDPlots.root",getHistName("Pruned","JP",true),300,500);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Extended PFchs","btagperfcomp_Pt300to500_PATTuple_v3_ExtPFchs.png",0, 1, 1E-3, 1,1);

  // Fat Jet Plots
  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Fat Jet CSV"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root","QCDPythia6_HiggsTagging.root",getHistName("Pruned","CSV",false),300,500);
  graphsPt300To500["Fat Jet JP"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root","QCDPythia6_HiggsTagging.root",getHistName("Pruned","JP",false),300,500);
  graphsPt300To500["Fat Jet JBP"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root","QCDPythia6_HiggsTagging.root",getHistName("Pruned","JBP",false),300,500);
  graphsPt300To500["Fat Jet CSV (Explicit JTA)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JP",false),300,500);
  graphsPt300To500["Fat Jet JP (Explicit JTA)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JP",false),300,500);
  graphsPt300To500["Fat Jet JBP (Explicit JTA)"]=getEfficiencyCurve("BprimeM1000_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JBP",false),300,500);

  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Fat Jet CSV"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root","QCDPythia6_HiggsTagging.root",getHistName("Pruned","CSV",false),700,1100);
  graphsPt700ToInf["Fat Jet JP"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root","QCDPythia6_HiggsTagging.root",getHistName("Pruned","JP",false),700,1100);
  graphsPt700ToInf["Fat Jet JBP"]=getEfficiencyCurve("BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root","QCDPythia6_HiggsTagging.root",getHistName("Pruned","JBP",false),700,1100);
  graphsPt700ToInf["Fat Jet CSV (Explicit JTA)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","CSV",false),700,1100);
  graphsPt700ToInf["Fat Jet JP (Explicit JTA)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JP",false),700,1100);
  graphsPt700ToInf["Fat Jet JBP (Explicit JTA)"]=getEfficiencyCurve("BprimeM1500_IVF_ExplicitJTA_Plots.root","QCDPythia6_IVF_ExplicitJTA_Plots.root",getHistName("Pruned","JBP",false),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt300to500_FatJetComparison.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","","btagperfcomp_Pt700toInf_FatJetComparison.png",0, 1, 1E-3, 1,1);

}
