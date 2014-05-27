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
void plotEfficiencyCurves(std::map< std::string,TGraph* > graphs, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle, const string& fOutputFile, const double fXmin, const double fXmax, const double fYmin, const double fYmax,const Int_t fLogy=0, const double fOPL1=0.244, const double fOPM1=0.679, const double fOPT1=0.898){

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
  legend->SetTextSize(0.025);

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
  c->SaveAs(fOutputFile.c_str());
  delete c;
  delete legend;
  delete bkg;
}

void efficiency_curves_comp_xrange(const string& fFileS1, const string& fFileS2, const string& fFileB1, const string& fFileB2,
                                   const string& fPlotS1, const string& fPlotS2, const string& fPlotB1, const string& fPlotB2,
                                   const double fXMin, const double fXMax, const string& fTitle, const string& fXAxisTitle, const string& fYAxisTitle,
                                   const string& fLeg1, const string& fLeg2, const double fXmin, const double fXmax, const double fYmin, const double fYmax,
                                   const string& fOutputFile, const Int_t fLogy=0, const double fOPL1=0.244, const double fOPM1=0.679, const double fOPT1=0.898,
                                   const double fOPL2=0.244, const double fOPM2=0.679, const double fOPT2=0.898)
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

  // loose operating point
  TGraph *g_eff_L = new TGraph(2);
  g_eff_L->SetName("g_eff_L");
  g_eff_L->SetMarkerStyle(31);
  g_eff_L->SetMarkerSize(1.5);

  g_eff_L->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(fOPL1),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(fOPL1),101)/denom_B_1));
  g_eff_L->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(fOPL2),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(fOPL2),101)/denom_B_2));

  // medium operating point
  TGraph *g_eff_M = new TGraph(2);
  g_eff_M->SetName("g_eff_L");
  g_eff_M->SetMarkerStyle(27);
  g_eff_M->SetMarkerSize(1.5);

  g_eff_M->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(fOPM1),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(fOPM1),101)/denom_B_1));
  g_eff_M->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(fOPM2),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(fOPM2),101)/denom_B_2));

  // tight operating point
  TGraph *g_eff_T = new TGraph(2);
  g_eff_T->SetName("g_eff_L");
  g_eff_T->SetMarkerStyle(30);
  g_eff_T->SetMarkerSize(1.5);

  g_eff_T->SetPoint(0,(h2_S_1->Integral(h2_S_1->GetXaxis()->FindBin(fXMin),h2_S_1->GetXaxis()->FindBin(fXMax),h2_S_1->GetYaxis()->FindBin(fOPT1),101)/denom_S_1),(h2_B_1->Integral(h2_B_1->GetXaxis()->FindBin(fXMin),h2_B_1->GetXaxis()->FindBin(fXMax),h2_B_1->GetYaxis()->FindBin(fOPT1),101)/denom_B_1));
  g_eff_T->SetPoint(1,(h2_S_2->Integral(h2_S_2->GetXaxis()->FindBin(fXMin),h2_S_2->GetXaxis()->FindBin(fXMax),h2_S_2->GetYaxis()->FindBin(fOPT2),101)/denom_S_2),(h2_B_2->Integral(h2_B_2->GetXaxis()->FindBin(fXMin),h2_B_2->GetXaxis()->FindBin(fXMax),h2_B_2->GetYaxis()->FindBin(fOPT2),101)/denom_B_2));


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
  l1.DrawLatex(0.14,0.96, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");
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

void makePlots()
{
  //for multiple plots on the same canvas

  std::map< std::string,TGraph* > graphsPt300To500;
  std::map< std::string,TGraph* > graphsPt700ToInf;

  graphsPt300To500["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1000plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),300,500);
  graphsPt300To500["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Pruned","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Filtered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Kt","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("MDBDRSFiltered","CSV",true),300,500);
  graphsPt300To500["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeM1000_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("KtBDRSFiltered","CSV",true),300,500);


  graphsPt700ToInf["Fat Jet CSV (BTV-13-001)"]=getEfficiencyCurve("BPrimeM1500plots.root","QCDPythia6Plots.root",getHistName("Pruned","CSVL",false),700,1100);
  graphsPt700ToInf["Subjet CSV (Pruned)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Pruned","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (Filtered)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Filtered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kT)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("Kt","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (MDBDRS)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("MDBDRSFiltered","CSV",true),700,1100);
  graphsPt700ToInf["Subjet CSV (kTBDRS)"]=getEfficiencyCurve("BprimeM1500_MDPlots.root","QCDPythia6_MDPlots.root",getHistName("KtBDRSFiltered","CSV",true),700,1100);

  plotEfficiencyCurves(graphsPt300To500,"#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","test.png",0, 1, 1E-3, 1,1);
  plotEfficiencyCurves(graphsPt700ToInf,"#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","test2.png",0, 1, 1E-3, 1,1);


// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000plots_newcuts.root",
//          "QCDPythia6Plots.root", "QCDPythia6Plots_NewCuts.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_Pruned.png", 1);
// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500plots_newcuts.root",
//          "QCDPythia6Plots.root", "QCDPythia6Plots_NewCuts.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_Pruned.png", 1);


// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000plots_newcuts.root",
//          "QCDPythia6Plots.root", "QCDPythia6Plots_NewCuts.root",
//          "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_Filtered.png", 1);
// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500plots_newcuts.root",
//          "QCDPythia6Plots.root", "QCDPythia6Plots_NewCuts.root",
//          "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_Filtered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000plots_newcuts.root",
//          "QCDPythia6Plots.root", "QCDPythia6Plots_NewCuts.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_Kt.png", 1);
// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500plots_newcuts.root",
//          "QCDPythia6Plots.root", "QCDPythia6Plots_NewCuts.root",
//          "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_Kt.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BprimeM1000_MDPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_MDPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_MDBDRSFiltered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BprimeM1500_MDPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_MDPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_MDBDRSFiltered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BprimeM1000_MDPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_MDPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_KtBDRSFiltered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BprimeM1500_MDPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_MDPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_KtBDRSFiltered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_IVF_KtBDRSFiltered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_IVF_KtBDRSFiltered.png", 1);
// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_IVF_MDBDRSFiltered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_IVF_MDBDRSFiltered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_IVF_Kt.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_IVF_Kt.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_IVF_Filtered.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_FilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_IVF_Filtered.png", 1);
// efficiency_curves_comp_xrange("BPrimeM1000plots.root", "BPrimeM1000_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_IVF_Pruned.png", 1);

// efficiency_curves_comp_xrange("BPrimeM1500plots.root", "BPrimeM1500_IVFPlots.root",
//          "QCDPythia6Plots.root", "QCDPythia6_IVFPlots.root",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetCSVL_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV (BTV-13-001)", "Subjet CSV",
//          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_IVF_Pruned.png", 1);



//   //CSV
//   //b-tagging efficiency vs mistag rate (CA8 jets, both subjets b-tagged)
//   efficiency_curves_comp_xrange("BprimeM1000_MDPlots.root", "BprimeM1000_MDPlots.root",
//                          "QCDPythia6_MDPlots.root", "QCDPythia6_MDPlots.root",
//                          "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_JetCSV_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//                          "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_JetCSV_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//                          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV", "Subjet CSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_KtBDRSFiltered.png", 1);

//   efficiency_curves_comp_xrange("BprimeM1500_MDPlots.root", "BprimeM1500_MDPlots.root",
//                          "QCDPythia6_MDPlots.root", "QCDPythia6_MDPlots.root",
//                          "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_JetCSV_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//                          "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_JetCSV_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",
//                          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet CSV", "Subjet CSV",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_KtBDRSFiltered.png", 1);

//   //JP
//   efficiency_curves_comp_xrange("BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root", "BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root",
//                          "QCDPythia6_HiggsTagging.root", "QCDPythia6_HiggsTagging.root",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJP_BosonMatched_JetMass;",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJP_BosonMatched_JetMass",
//                          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet JP", "Subjet JP",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_JPPruned.png", 1,0.275,0.545,0.790,0.275,0.545,0.790);
//   efficiency_curves_comp_xrange("BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root", "BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root",
//                          "QCDPythia6_HiggsTagging.root", "QCDPythia6_HiggsTagging.root",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJP_BosonMatched_JetMass;",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJP_BosonMatched_JetMass",
//                          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet JP", "Subjet JP",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_JPPruned.png", 1,0.275,0.545,0.790,0.275,0.545,0.790);

// //JPB
//   efficiency_curves_comp_xrange("BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root", "BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root",
//                          "QCDPythia6_HiggsTagging.root", "QCDPythia6_HiggsTagging.root",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJBP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJBP_BosonMatched_JetMass;",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJBP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJBP_BosonMatched_JetMass",
//                          300, 500, "#splitline{CA R=0.8, 300<p_{T}<500 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet JBP", "Subjet JBP",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt300to500_JetMass_JBPPruned.png", 1,1.33, 2.55, 3.74, 1.33, 2.55, 3.74);
//   efficiency_curves_comp_xrange("BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root", "BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root",
//                          "QCDPythia6_HiggsTagging.root", "QCDPythia6_HiggsTagging.root",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJBP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJBP_BosonMatched_JetMass;",
//                          "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_JetJBP_BosonMatched_JetMass", "jetAnalyzerCA8FatJets_PrunedSubjets/h2_JetPt_SubJetMinJBP_BosonMatched_JetMass",
//                          700, 1100, "#splitline{CA R=0.8, p_{T}>700 GeV/c}{75<m_{jet}<135 GeV/c^{2} (pruned)}", "b-tagging efficiency (H(120)#rightarrowb#bar{b})", "Misidentification probability (QCD)","Fat jet JBP", "Subjet JBP",
//                          0, 1, 1E-3, 1, "b-tag_eff_vs_mistag_CA8_QCDBkg_Pt700toInf_JetMass_JBPPruned.png", 1,1.33, 2.55, 3.74, 1.33, 2.55, 3.74);
}
