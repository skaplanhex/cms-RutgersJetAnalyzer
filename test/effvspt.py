from ROOT import *
from numpy import array

gROOT.SetBatch(kTRUE)

CSVL = 0.244
CSVM = 0.679
JPL = 0.275
JPM = 0.545
JBPL = 1.33
JBPM = 2.55

def makePlotsJP(analyzer,grooming,isJPB):
	algo = "JP"
	l=JPL
	m=JPM
	if isJPB:
		algo="JBP"
		l=JBPL
		m=JBPM
	########################### SIGNAL #########################################

	f1000 = TFile("BprimeBprimeToBHBHinc_M-1000_HiggsTagging.root")
	f1500 = TFile("BprimeBprimeToBHBHinc_M-1500_HiggsTagging.root ")
	fqcd = TFile("QCDPythia6_HiggsTagging.root")

	h1000j = f1000.Get(analyzer+"/h2_JetPt_Jet"+algo+"_BosonMatched_JetMass")
	h1000sj = f1000.Get(analyzer+"/h2_JetPt_SubJetMin"+algo+"_BosonMatched_JetMass")
	h1500j = f1500.Get(analyzer+"/h2_JetPt_Jet"+algo+"_BosonMatched_JetMass")
	h1500sj = f1500.Get(analyzer+"/h2_JetPt_SubJetMin"+algo+"_BosonMatched_JetMass")
	hqcdj = fqcd.Get(analyzer+"/h2_JetPt_Jet"+algo+"_BosonMatched_JetMass")
	hqcdsj = fqcd.Get(analyzer+"/h2_JetPt_SubJetMin"+algo+"_BosonMatched_JetMass")

	c = TCanvas()
	graphL = TGraph()
	graphM = TGraph()
	graphSJL = TGraph()
	graphSJM = TGraph()
	#pt binned 0-1000 with each bin = 4 GeV --> 250 bins.  300 GeV starts at bin 75
	#CSV binned 0-1 with 100 bins
	for i in range(7):
		pt=350+100*i
		proj = h1000j.ProjectionY("proj",76+25*i,100+25*i)
		projsj = h1000sj.ProjectionY("projsj",76+25*i,100+25*i)
		if (pt > 700):
			proj = h1500j.ProjectionY("proj",76+25*i,100+25*i)
			projsj = h1500sj.ProjectionY("projsj",76+25*i,100+25*i)
		effL = float( proj.Integral( proj.GetXaxis().FindBin(l),100))/float(proj.Integral())
		effM = float( proj.Integral(proj.GetXaxis().FindBin(m),100))/float(proj.Integral())
		effSJL = float( projsj.Integral(projsj.GetXaxis().FindBin(l),100))/float(proj.Integral())
		effSJM = float( projsj.Integral(proj.GetXaxis().FindBin(m),100))/float(proj.Integral())
		graphL.SetPoint(i,pt,effL)
		graphM.SetPoint(i,pt,effM)
		graphSJL.SetPoint(i,pt,effSJL)
		graphSJM.SetPoint(i,pt,effSJM)
	graphL.GetYaxis().SetRangeUser(0.,1.)
	graphM.GetYaxis().SetRangeUser(0.,1.)
	graphSJL.GetYaxis().SetRangeUser(0.,1.)
	graphSJM.GetYaxis().SetRangeUser(0.,1.)
	graphL.GetXaxis().SetTitle("Fat Jet pT (GeV/c)")
	graphL.GetXaxis().SetTitleSize(0.05)
	graphL.GetXaxis().SetTitleOffset(0.85)
	graphL.GetYaxis().SetTitle("b-tagging efficiency (h#rightarrowb#bar{b})")
	graphL.GetYaxis().SetTitleSize(0.05)
	graphL.GetYaxis().SetTitleOffset(0.85)
	graphL.SetLineColor(kRed)
	graphM.SetLineColor(kGreen)
	graphSJL.SetLineColor(kBlue)
	graphSJM.SetLineColor(kMagenta)
	for iGraph in (graphL,graphM,graphSJL,graphSJM):
		iGraph.SetLineWidth(2)
	graphL.Draw("AC")
	graphM.Draw("C")
	graphSJL.Draw("C")
	graphSJM.Draw("C")
	legend = TLegend(.16,.12,.36,.25);
	legend.SetBorderSize(0);
	legend.SetFillColor(0);
	legend.SetFillStyle(0);
	legend.SetTextFont(42);
	legend.SetTextSize(0.04);
	legend.AddEntry(graphL, "Fat Jet " +algo+"L","l");
	legend.AddEntry(graphM, "Fat Jet " +algo+"M","l");
	legend.AddEntry(graphSJL, "Subjet "+algo+"L","l");
	legend.AddEntry(graphSJM, "Subjet "+algo+"M","l");
	legend.Draw();

	l1 = TLatex();
	l1.SetNDC();
	l1.SetTextAlign(12);
	l1.SetTextSize(0.045);
	l1.SetTextFont(62);
	l1.DrawLatex(0.14,0.925, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");

	l1.SetTextFont(42);
	l1.SetTextSize(0.04);
	l1.DrawLatex(0.63,0.20, "75 < m_{fat jet} < 135");
	l1.DrawLatex(0.63,0.15, grooming+" Subjets");
	c.SaveAs("EffVsPt_Signal_"+grooming+algo+".png")
	c.Clear()

	########################### BACKGROUND #########################################

	graphL = TGraph()
	graphM = TGraph()
	graphSJL = TGraph()
	graphSJM = TGraph()
	#pt binned 0-1000 with each bin = 4 GeV --> 250 bins.  300 GeV starts at bin 75
	#CSV binned 0-1 with 100 bins
	for i in range(7):
		pt=350+100*i
		proj = hqcdj.ProjectionY("proj",76+25*i,100+25*i)
		projsj = hqcdsj.ProjectionY("projsj",76+25*i,100+25*i)
		effL = float( proj.Integral(proj.GetXaxis().FindBin(l),100))/float(proj.Integral())
		effM = float( proj.Integral(proj.GetXaxis().FindBin(m),100))/float(proj.Integral())
		effSJL = float( projsj.Integral(projsj.GetXaxis().FindBin(l),100))/float(proj.Integral())
		effSJM = float( projsj.Integral(projsj.GetXaxis().FindBin(m),100))/float(proj.Integral())
		graphL.SetPoint(i,pt,effL)
		graphM.SetPoint(i,pt,effM)
		graphSJL.SetPoint(i,pt,effSJL)
		graphSJM.SetPoint(i,pt,effSJM)
	graphL.GetYaxis().SetRangeUser(0.,0.75)
	graphM.GetYaxis().SetRangeUser(0.,0.75)
	graphSJL.GetYaxis().SetRangeUser(0.,0.75)
	graphSJM.GetYaxis().SetRangeUser(0.,0.75)
	graphL.GetXaxis().SetTitle("Fat Jet pT (GeV/c)")
	graphL.GetXaxis().SetTitleSize(0.05)
	graphL.GetXaxis().SetTitleOffset(0.85)
	graphL.GetYaxis().SetTitle("Mistag Rate (QCD)")
	graphL.GetYaxis().SetTitleSize(0.05)
	graphL.GetYaxis().SetTitleOffset(0.85)
	graphL.SetLineColor(kRed)
	graphM.SetLineColor(kGreen)
	graphSJL.SetLineColor(kBlue)
	graphSJM.SetLineColor(kMagenta)
	for iGraph in (graphL,graphM,graphSJL,graphSJM):
		iGraph.SetLineWidth(2)
	graphL.Draw("AC")
	graphM.Draw("C")
	graphSJL.Draw("C")
	graphSJM.Draw("C")
	legend = TLegend(.16,.67,.36,.80);
	legend.SetBorderSize(0);
	legend.SetFillColor(0);
	legend.SetFillStyle(0);
	legend.SetTextFont(42);
	legend.SetTextSize(0.04);
	legend.AddEntry(graphL, "Fat Jet "+algo+"L","l");
	legend.AddEntry(graphM, "Fat Jet "+algo+"M","l");
	legend.AddEntry(graphSJL, "Subjet "+algo+"L","l");
	legend.AddEntry(graphSJM, "Subjet "+algo+"M","l");
	legend.Draw();

	l1 = TLatex();
	l1.SetNDC();
	l1.SetTextAlign(12);
	l1.SetTextSize(0.045);
	l1.SetTextFont(62);
	l1.DrawLatex(0.14,0.925, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");

	l1.SetTextFont(42);
	l1.SetTextSize(0.04);
	l1.DrawLatex(0.63,0.68, grooming+" Subjets");
	l1.DrawLatex(0.63,0.73, "75 < m_{fat jet} < 135");
	c.SaveAs("EffVsPt_Background_"+grooming+algo+".png")
	c.Clear()

def makePlots(analyzer,grooming,massDrop):
	########################### SIGNAL #########################################

	f1000 = TFile("oldstuff/BPrimeM1000plots_newcuts.root")
	f1500 = TFile("oldstuff/BPrimeM1500plots_newcuts.root")
	fqcd = TFile("oldstuff/QCDPythia6Plots_NewCuts.root")

	h1000j = f1000.Get(analyzer+"/h2_JetPt_JetCSVL_BosonMatched_JetMass")
	h1000sj = f1000.Get(analyzer+"/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass")
	h1500j = f1500.Get(analyzer+"/h2_JetPt_JetCSVL_BosonMatched_JetMass")
	h1500sj = f1500.Get(analyzer+"/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass")
	hqcdj = fqcd.Get(analyzer+"/h2_JetPt_JetCSVL_BosonMatched_JetMass")
	hqcdsj = fqcd.Get(analyzer+"/h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass")

	if massDrop:
		f1000 = TFile("BprimeM1000_MDPlots.root")
		f1500 = TFile("BprimeM1500_MDPlots.root")
		fqcd = TFile("QCDPythia6_MDPlots.root")

		h1000j = f1000.Get(analyzer+"/h2_JetPt_JetCSV_BosonMatched_JetMass")
		h1000sj = f1000.Get(analyzer+"/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass")
		h1500j = f1500.Get(analyzer+"/h2_JetPt_JetCSV_BosonMatched_JetMass")
		h1500sj = f1500.Get(analyzer+"/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass")
		hqcdj = fqcd.Get(analyzer+"/h2_JetPt_JetCSV_BosonMatched_JetMass")
		hqcdsj = fqcd.Get(analyzer+"/h2_JetPt_SubJetMinCSV_BosonMatched_JetMass")

	c = TCanvas()
	graphL = TGraph()
	graphM = TGraph()
	graphSJL = TGraph()
	graphSJM = TGraph()
	#pt binned 0-1000 with each bin = 4 GeV --> 250 bins.  300 GeV starts at bin 75
	#CSV binned 0-1 with 100 bins
	for i in range(7):
		pt=350+100*i
		proj = h1000j.ProjectionY("proj",76+25*i,100+25*i)
		projsj = h1000sj.ProjectionY("projsj",76+25*i,100+25*i)
		if (pt > 700):
			proj = h1500j.ProjectionY("proj",76+25*i,100+25*i)
			projsj = h1500sj.ProjectionY("projsj",76+25*i,100+25*i)
		effL = float( proj.Integral(24,100))/float(proj.Integral())
		effM = float( proj.Integral(68,100))/float(proj.Integral())
		effSJL = float( projsj.Integral(24,100))/float(proj.Integral())
		effSJM = float( projsj.Integral(68,100))/float(proj.Integral())
		graphL.SetPoint(i,pt,effL)
		graphM.SetPoint(i,pt,effM)
		graphSJL.SetPoint(i,pt,effSJL)
		graphSJM.SetPoint(i,pt,effSJM)
	graphL.GetYaxis().SetRangeUser(0.,1.)
	graphM.GetYaxis().SetRangeUser(0.,1.)
	graphSJL.GetYaxis().SetRangeUser(0.,1.)
	graphSJM.GetYaxis().SetRangeUser(0.,1.)
	graphL.GetXaxis().SetTitle("Fat Jet pT (GeV/c)")
	graphL.GetXaxis().SetTitleSize(0.05)
	graphL.GetXaxis().SetTitleOffset(0.85)
	graphL.GetYaxis().SetTitle("b-tagging efficiency (h#rightarrowb#bar{b})")
	graphL.GetYaxis().SetTitleSize(0.05)
	graphL.GetYaxis().SetTitleOffset(0.85)
	graphL.SetLineColor(kRed)
	graphM.SetLineColor(kGreen)
	graphSJL.SetLineColor(kBlue)
	graphSJM.SetLineColor(kMagenta)
	for iGraph in (graphL,graphM,graphSJL,graphSJM):
		iGraph.SetLineWidth(2)
	graphL.Draw("AC")
	graphM.Draw("C")
	graphSJL.Draw("C")
	graphSJM.Draw("C")
	legend = TLegend(.16,.12,.36,.25);
	legend.SetBorderSize(0);
	legend.SetFillColor(0);
	legend.SetFillStyle(0);
	legend.SetTextFont(42);
	legend.SetTextSize(0.04);
	legend.AddEntry(graphL, "Fat Jet CSVL","l");
	legend.AddEntry(graphM, "Fat Jet CSVM","l");
	legend.AddEntry(graphSJL, "Subjet CSVL","l");
	legend.AddEntry(graphSJM, "Subjet CSVM","l");
	legend.Draw();

	l1 = TLatex();
	l1.SetNDC();
	l1.SetTextAlign(12);
	l1.SetTextSize(0.045);
	l1.SetTextFont(62);
	l1.DrawLatex(0.14,0.925, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");

	l1.SetTextFont(42);
	l1.SetTextSize(0.04);
	l1.DrawLatex(0.63,0.20, "75 < m_{fat jet} < 135");
	l1.DrawLatex(0.63,0.15, grooming+" Subjets");
	c.SaveAs("EffVsPt_Signal_"+grooming+".png")
	c.Clear()

	########################### BACKGROUND #########################################

	graphL = TGraph()
	graphM = TGraph()
	graphSJL = TGraph()
	graphSJM = TGraph()
	#pt binned 0-1000 with each bin = 4 GeV --> 250 bins.  300 GeV starts at bin 75
	#CSV binned 0-1 with 100 bins
	for i in range(7):
		pt=350+100*i
		proj = hqcdj.ProjectionY("proj",76+25*i,100+25*i)
		projsj = hqcdsj.ProjectionY("projsj",76+25*i,100+25*i)
		effL = float( proj.Integral(24,100))/float(proj.Integral())
		effM = float( proj.Integral(68,100))/float(proj.Integral())
		effSJL = float( projsj.Integral(24,100))/float(proj.Integral())
		effSJM = float( projsj.Integral(68,100))/float(proj.Integral())
		graphL.SetPoint(i,pt,effL)
		graphM.SetPoint(i,pt,effM)
		graphSJL.SetPoint(i,pt,effSJL)
		graphSJM.SetPoint(i,pt,effSJM)
	graphL.GetYaxis().SetRangeUser(0.,0.2)
	graphM.GetYaxis().SetRangeUser(0.,0.2)
	graphSJL.GetYaxis().SetRangeUser(0.,0.2)
	graphSJM.GetYaxis().SetRangeUser(0.,0.2)
	graphL.GetXaxis().SetTitle("Fat Jet pT (GeV/c)")
	graphL.GetXaxis().SetTitleSize(0.05)
	graphL.GetXaxis().SetTitleOffset(0.85)
	graphL.GetYaxis().SetTitle("Mistag Rate (QCD)")
	graphL.GetYaxis().SetTitleSize(0.05)
	graphL.GetYaxis().SetTitleOffset(0.85)
	graphL.SetLineColor(kRed)
	graphM.SetLineColor(kGreen)
	graphSJL.SetLineColor(kBlue)
	graphSJM.SetLineColor(kMagenta)
	for iGraph in (graphL,graphM,graphSJL,graphSJM):
		iGraph.SetLineWidth(2)
	graphL.Draw("AC")
	graphM.Draw("C")
	graphSJL.Draw("C")
	graphSJM.Draw("C")
	legend = TLegend(.16,.67,.36,.80);
	legend.SetBorderSize(0);
	legend.SetFillColor(0);
	legend.SetFillStyle(0);
	legend.SetTextFont(42);
	legend.SetTextSize(0.04);
	legend.AddEntry(graphL, "Fat Jet CSVL","l");
	legend.AddEntry(graphM, "Fat Jet CSVM","l");
	legend.AddEntry(graphSJL, "Subjet CSVL","l");
	legend.AddEntry(graphSJM, "Subjet CSVM","l");
	legend.Draw();

	l1 = TLatex();
	l1.SetNDC();
	l1.SetTextAlign(12);
	l1.SetTextSize(0.045);
	l1.SetTextFont(62);
	l1.DrawLatex(0.14,0.925, "CMS Simulation Preliminary, #sqrt{s} = 8 TeV");

	l1.SetTextFont(42);
	l1.SetTextSize(0.04);
	l1.DrawLatex(0.63,0.68, grooming+" Subjets");
	l1.DrawLatex(0.63,0.73, "75 < m_{fat jet} < 135");
	c.SaveAs("EffVsPt_Background_"+grooming+".png")
	c.Clear()

makePlotsJP("jetAnalyzerCA8FatJets_PrunedSubjets","Pruned",False)
makePlotsJP("jetAnalyzerCA8FatJets_PrunedSubjets","Pruned",True)

# makePlots("jetAnalyzerCA8FatJets_PrunedSubjets","Pruned",False)
# makePlots("jetAnalyzerCA8FatJets_FilteredSubjets","Filtered",False)
# makePlots("jetAnalyzerCA8FatJets_BDRSFilteredSubjets","BDRSFiltered",False)
# makePlots("jetAnalyzerCA8FatJets_KtSubjets","Kt",False)
# makePlots("jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets","MassDrop_BDRSFiltered",True)
# makePlots("jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets","Kt_BDRSFiltered",True)
