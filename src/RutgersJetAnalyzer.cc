// -*- C++ -*-
//
// Package:    RutgersJetAnalyzer
// Class:      RutgersJetAnalyzer
//
/**\class RutgersJetAnalyzer RutgersJetAnalyzer.cc RutgersSandbox/RutgersJetAnalyzer/src/RutgersJetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Dinko Ferencek
//         Created:  Fri Jul 20 12:32:38 CDT 2012
// $Id: RutgersJetAnalyzer.cc,v 1.7.2.15 2012/12/10 03:29:02 mzientek Exp $
//
//


// system include files
#include <memory>

// FastJet include files
#include "fastjet/PseudoJet.hh"
// N-subjettiness include files
#include "RutgersSandbox/RutgersJetAnalyzer/plugins/Njettiness.hh"
// user include files
#include <boost/shared_ptr.hpp>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TH2D.h"

//
// class declaration
//

struct ordering_dR {
    const pat::Jet* mJet;
    ordering_dR(const pat::Jet* fJet) : mJet(fJet) {}
    bool operator ()(const pat::Jet* const& a, const pat::Jet* const& b) {
      return reco::deltaR( mJet->p4(), a->p4() ) < reco::deltaR( mJet->p4(), b->p4() );
    }
};

struct ordering_Pt {
    const std::string mCorrLevel;
    ordering_Pt(const std::string fCorrLevel) : mCorrLevel(fCorrLevel) {}
    bool operator ()(const pat::Jet* const& a, const pat::Jet* const& b) {
      if( mCorrLevel=="Uncorrected" )
        return a->correctedJet("Uncorrected").pt() > b->correctedJet("Uncorrected").pt();
      else
        return a->pt() > b->pt();
    }
};

class RutgersJetAnalyzer : public edm::EDAnalyzer {
public:
    explicit RutgersJetAnalyzer(const edm::ParameterSet&);
    ~RutgersJetAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    // typedefs
    typedef std::vector<pat::Jet> PatJetCollection;


private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    // ----------member data ---------------------------
    const bool          useEventWeight;
    const edm::InputTag genParticleTag;
    const edm::InputTag jetsTag;
    const bool          useGroomedJets;
    const edm::InputTag groomedJetsTag;
    const bool          useSubJets;
    const edm::InputTag subJetsTag;
    const std::string   subJetMode;
    const edm::InputTag pvTag;
    const double        jetRadius;          // radius for jet clustering
    const bool          doBosonMatching;    // parameter for deciding if matching is on or off
    const double        bosonMatchingRadius;
    const int           bosonPdgId;
    const bool          applyBosonIsolation;
    const bool          doBosonDecayProdSelection;
    const std::vector<int> bosonDecayProdPdgIds;
    const bool          useMassDrop;
    const double        jetPtMin;
    const unsigned      jetPtBins;
    const double        jetPtBinWidth;
    const double        jetAbsEtaMax;
    const double        jetMassMin;
    const double        jetMassMax;
    const double        nsubjCut;
    bool                useGroomedJetSubstr;
    bool                useUncorrMassForMassDrop;
    std::string		bdiscriminator;
    bool		findGluonSplitting;
    bool		doJetFlavor;
    const std::vector<int> jetFlavorPdgId;
    bool		findMatrixElement;

    Njettiness nsubjettinessCalculator;

    edm::Service<TFileService> fs;

    TH1D *h1_nPV;

    TH1D *h1_BosonPt;
    TH1D *h1_BosonEta;
    TH1D *h1_BosonPt_Isolated;
    TH1D *h1_BosonEta_Isolated;
    TH1D *h1_BosonPt_DecaySel;
    TH1D *h1_BosonEta_DecaySel;
    TH1D *h1_BosonPt_Matched;
    TH1D *h1_BosonPt_DecayProdMatched;

    TH2D *h2_BosonPt_dRdecay;

    TH1D *h1_JetPt;
    TH1D *h1_JetPt_BosonMatched;
    TH1D *h1_JetPt_BosonMatched_JetMass;
    TH1D *h1_JetPt_BosonMatched_JetMass_CSVL;
    TH1D *h1_JetPt_BosonMatched_JetMass_CSVM;
    TH1D *h1_JetPt_BosonMatched_JetMass_SubJetCSVL;
    TH1D *h1_JetPt_BosonMatched_JetMass_SubJetCSVM;
    TH1D *h1_JetPt_BosonMatched_JetMass_DoubleB;
    TH1D *h1_JetPt_BosonMatched_JetMass_Nsubj_CSVL;
    TH1D *h1_JetPt_BosonMatched_JetMass_Nsubj_CSVM;
    TH1D *h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVL;
    TH1D *h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVM;
    TH1D *h1_JetPt_BosonMatched_JetMass_Nsubj_DoubleB;
    TH1D *h1_JetPt_BosonDecayProdMatched;
    TH1D *h1_JetPt_BosonDecayProdMatched_JetMass;
    TH1D *h1_JetEta;
    TH1D *h1_JetEta_BosonMatched;
    TH1D *h1_JetEta_BosonMatched_JetMass;

    TH2D *h2_JetPt_JetPtOverBosonPt;
    TH2D *h2_JetPt_JetPtOverGenJetPt;
    TH2D *h2_JetPt_JetPtOverGenJetPt_BosonMatched;
    TH2D *h2_JetPt_JetMass;
    TH2D *h2_JetPt_JetMass_BosonMatched;

    TH1D *h1_JetCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_SubJetCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_JetDoubleBDiscr_BosonMatched_JetMass;

    std::map<std::string, TH2D*> h2_nPV_JetMass_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_MassDrop_Pt;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RutgersJetAnalyzer::RutgersJetAnalyzer(const edm::ParameterSet& iConfig) :

  useEventWeight(iConfig.getParameter<bool>("UseEventWeight")),
  genParticleTag(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
  jetsTag(iConfig.getParameter<edm::InputTag>("JetsTag")),
  useGroomedJets(iConfig.getParameter<bool>("UseGroomedJets")),
  groomedJetsTag(iConfig.getParameter<edm::InputTag>("GroomedJetsTag")),
  useSubJets(iConfig.getParameter<bool>("UseSubJets")),
  subJetsTag(iConfig.getParameter<edm::InputTag>("SubJetsTag")),
  subJetMode(iConfig.getParameter<std::string>("SubJetMode")),
  pvTag(iConfig.getParameter<edm::InputTag>("PvTag")),
  jetRadius(iConfig.getParameter<double>("JetRadius")),
  doBosonMatching(iConfig.getParameter<bool>("DoBosonMatching")),
  bosonMatchingRadius(iConfig.getParameter<double>("BosonMatchingRadius")),
  bosonPdgId(iConfig.getParameter<int>("BosonPdgId")),
  applyBosonIsolation(iConfig.getParameter<bool>("ApplyBosonIsolation")),
  doBosonDecayProdSelection(iConfig.getParameter<bool>("DoBosonDecayProdSelection")),
  bosonDecayProdPdgIds(iConfig.getParameter<std::vector<int> >("BosonDecayProdPdgIds")),
  useMassDrop(iConfig.getParameter<bool>("UseMassDrop")),
  jetPtMin(iConfig.getParameter<double>("JetPtMin")),
  jetPtBins(iConfig.getParameter<unsigned>("JetPtBins")),
  jetPtBinWidth(iConfig.getParameter<double>("JetPtBinWidth")),
  jetAbsEtaMax(iConfig.getParameter<double>("JetAbsEtaMax")),
  jetMassMin(iConfig.getParameter<double>("JetMassMin")),
  jetMassMax(iConfig.getParameter<double>("JetMassMax")),
  nsubjCut(iConfig.getParameter<double>("NsubjCut")),
  useGroomedJetSubstr(false),
  useUncorrMassForMassDrop(false),
  bdiscriminator(iConfig.getParameter<std::string>("Bdiscriminator")),
  findGluonSplitting(iConfig.getParameter<bool>("FindGluonSplitting")),
  doJetFlavor(iConfig.getParameter<bool>("DoJetFlavor")),
  jetFlavorPdgId(iConfig.getParameter<std::vector<int> >("JetFlavorPdgId")),
  nsubjettinessCalculator(Njettiness::onepass_kt_axes, NsubParameters(1.0, jetRadius, jetRadius)),
  findMatrixElement(iConfig.getParameter<bool>("FindMatrixElement"))

{
    //now do what ever initialization is needed
    if ( iConfig.exists("UseGroomedJetSubstructure") )
      useGroomedJetSubstr = iConfig.getParameter<bool>("UseGroomedJetSubstructure");
    if ( iConfig.exists("UseUncorrMassForMassDrop") )
      useUncorrMassForMassDrop = iConfig.getParameter<bool>("UseUncorrMassForMassDrop");

    int pvBins=51;
    double pvMin=-0.5, pvMax=50.5;
    int ptBins=1000;
    double ptMin=0., ptMax=1000.;
    int dRBins=100;
    double dRMin=0., dRMax=5.;
    int etaBins=160;
    double etaMin=-4., etaMax=4.;
    int massBins=400;
    int massMin=0., massMax=400.;
    int tauBins=100;
    double tauMin=0., tauMax=1.;
    int massDropBins=100;
    double massDropMin=0., massDropMax=1.;

    h1_nPV = fs->make<TH1D>("h1_nPV","PV Multiplicity;nPV;",pvBins,pvMin,pvMax);
    h1_nPV->Sumw2();
    h1_nPV->SetDefaultSumw2(kTRUE);

    h1_BosonPt           = fs->make<TH1D>("h1_BosonPt",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_BosonEta          = fs->make<TH1D>("h1_BosonEta",";#eta;",etaBins,etaMin,etaMax);
    h1_BosonPt_Isolated  = fs->make<TH1D>("h1_BosonPt_Isolated",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_BosonEta_Isolated = fs->make<TH1D>("h1_BosonEta_Isolated",";#eta;",etaBins,etaMin,etaMax);
    h1_BosonPt_DecaySel  = fs->make<TH1D>("h1_BosonPt_DecaySel",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_BosonEta_DecaySel = fs->make<TH1D>("h1_BosonEta_DecaySel",";#eta;",etaBins,etaMin,etaMax);
    h1_BosonPt_Matched  = fs->make<TH1D>("h1_BosonPt_Matched",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_BosonPt_DecayProdMatched  = fs->make<TH1D>("h1_BosonPt_DecayProdMatched",";p_{T} [GeV];",ptBins,ptMin,ptMax);

    h2_BosonPt_dRdecay = fs->make<TH2D>("h2_BosonPt_dRdecay",";p_{T} [GeV];#DeltaR",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);

    h1_JetPt = fs->make<TH1D>("h1_JetPt",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched = fs->make<TH1D>("h1_JetPt_BosonMatched",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_CSVL = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_CSVL",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_CSVM = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_CSVM",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_SubJetCSVL = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_SubJetCSVL",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_SubJetCSVM = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_SubJetCSVM",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_DoubleB = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_DoubleB",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_Nsubj_CSVL = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_Nsubj_CSVL",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_Nsubj_CSVM = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_Nsubj_CSVM",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVL = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVL",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVM = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVM",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_Nsubj_DoubleB = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_Nsubj_DoubleB",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonDecayProdMatched = fs->make<TH1D>("h1_JetPt_BosonDecayProdMatched",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonDecayProdMatched_JetMass  = fs->make<TH1D>("h1_JetPt_BosonDecayProdMatched_JetMass",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetEta = fs->make<TH1D>("h1_JetEta",";#eta;",etaBins,etaMin,etaMax);
    h1_JetEta_BosonMatched = fs->make<TH1D>("h1_JetEta_BosonMatched",";#eta;",etaBins,etaMin,etaMax);
    h1_JetEta_BosonMatched_JetMass = fs->make<TH1D>("h1_JetEta_BosonMatched_JetMass",";#eta;",etaBins,etaMin,etaMax);

    h2_JetPt_JetPtOverBosonPt = fs->make<TH2D>("h2_JetPt_JetPtOverBosonPt",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{boson}",ptBins/2,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetPtOverGenJetPt = fs->make<TH2D>("h2_JetPt_JetPtOverGenJetPt",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{genjet}",ptBins/2,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetPtOverGenJetPt_BosonMatched = fs->make<TH2D>("h2_JetPt_JetPtOverGenJetPt_BosonMatched",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{genjet}",ptBins/2,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetMass = fs->make<TH2D>("h2_JetPt_JetMass",";p_{T} [GeV];m_{jet} [GeV]",ptBins/2,ptMin,ptMax,massBins/2,massMin,massMax);
    h2_JetPt_JetMass_BosonMatched = fs->make<TH2D>("h2_JetPt_JetMass_BosonMatched",";p_{T} [GeV];m_{jet} [GeV]",ptBins/2,ptMin,ptMax,massBins/2,massMin,massMax);

    h1_JetCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetCSVDiscr_BosonMatched_JetMass",";Jet CSV Discr;",100,0.,1.);
    h1_SubJetCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_SubJetCSVDiscr_BosonMatched_JetMass",";SubJet CSV Discr;",100,0.,1.);
    h1_JetDoubleBDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetDoubleBDiscr_BosonMatched_JetMass",";Jet DoubleB Discr;",100,0.,10.);

    for(unsigned i=0; i<=(jetPtBins+1); ++i)
    {
      std::string suffix, title;

      if(i==0)
      {
        suffix = Form("%.0ftoInf",(jetPtMin + jetPtBinWidth*i));
        title = Form("p_{T}>%.0f GeV",(jetPtMin + jetPtBinWidth*i));

        h2_nPV_JetMass_Pt[suffix]  = fs->make<TH2D>(("h2_nPV_JetMass_Pt" + suffix).c_str(),(title + ";nPV;m_{jet} [GeV]").c_str(),pvBins,pvMin,pvMax,massBins,massMin,massMax);
        h2_nPV_tau1_Pt[suffix]     = fs->make<TH2D>(("h2_nPV_tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2_Pt[suffix]     = fs->make<TH2D>(("h2_nPV_tau2_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nPV_tau2tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}/#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_MassDrop_Pt[suffix] = fs->make<TH2D>(("h2_nPV_MassDrop_Pt" + suffix).c_str(),(title + ";nPV;m_{subjet1}/m_{jet}").c_str(),pvBins,pvMin,pvMax,massDropBins,massDropMin,massDropMax);
      }
      else if(i==(jetPtBins+1))
      {
        suffix = Form("%.0ftoInf",(jetPtMin + jetPtBinWidth*(i-1)));
        title = Form("p_{T}>%.0f GeV",(jetPtMin + jetPtBinWidth*(i-1)));

        h2_nPV_JetMass_Pt[suffix]  = fs->make<TH2D>(("h2_nPV_JetMass_Pt" + suffix).c_str(),(title + ";nPV;m_{jet} [GeV]").c_str(),pvBins,pvMin,pvMax,massBins,massMin,massMax);
        h2_nPV_tau1_Pt[suffix]     = fs->make<TH2D>(("h2_nPV_tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2_Pt[suffix]     = fs->make<TH2D>(("h2_nPV_tau2_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nPV_tau2tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}/#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_MassDrop_Pt[suffix] = fs->make<TH2D>(("h2_nPV_MassDrop_Pt" + suffix).c_str(),(title + ";nPV;m_{subjet1}/m_{jet}").c_str(),pvBins,pvMin,pvMax,massDropBins,massDropMin,massDropMax);
      }
      else
      {
        suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*(i-1)),(jetPtMin + jetPtBinWidth*i));
        title = Form("%.0f<p_{T}<%.0f GeV",(jetPtMin + jetPtBinWidth*(i-1)),(jetPtMin + jetPtBinWidth*i));

        h2_nPV_JetMass_Pt[suffix]  = fs->make<TH2D>(("h2_nPV_JetMass_Pt" + suffix).c_str(),(title + ";nPV;m_{jet} [GeV]").c_str(),pvBins,pvMin,pvMax,massBins,massMin,massMax);
        h2_nPV_tau1_Pt[suffix]     = fs->make<TH2D>(("h2_nPV_tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2_Pt[suffix]     = fs->make<TH2D>(("h2_nPV_tau2_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nPV_tau2tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}/#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_MassDrop_Pt[suffix] = fs->make<TH2D>(("h2_nPV_MassDrop_Pt" + suffix).c_str(),(title + ";nPV;m_{subjet1}/m_{jet}").c_str(),pvBins,pvMin,pvMax,massDropBins,massDropMin,massDropMax);
      }
    }
}


RutgersJetAnalyzer::~RutgersJetAnalyzer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RutgersJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(genParticleTag,genParticles);

    edm::Handle<PatJetCollection> jets;
    iEvent.getByLabel(jetsTag,jets);

    edm::Handle<PatJetCollection> groomedJets;
    if( useGroomedJets ) iEvent.getByLabel(groomedJetsTag,groomedJets);

    edm::Handle<PatJetCollection> subJets;
    if( useSubJets ) iEvent.getByLabel(subJetsTag,subJets);

    edm::Handle<reco::VertexCollection> PVs;
    iEvent.getByLabel(pvTag,PVs);

    double eventWeight = 1.;
    if( !iEvent.isRealData() && useEventWeight )
    {
      edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
      iEvent.getByLabel("generator", genEvtInfoProduct);

      eventWeight = genEvtInfoProduct->weight();
    }

    int nPV = PVs->size();
    // fill histogram of the number of reconstructed PVs
    h1_nPV->Fill(nPV, eventWeight);

    // vector of pointers to status=3 b' decay products
    std::vector<const reco::GenParticle*> bPrimeDecayProducts;
    // vector of pointers to bosons
    std::vector<const reco::GenParticle*> bosons;
    // map to vectors of pointers to boson decay products
    std::map<const reco::GenParticle*,std::vector<const reco::Candidate*> > decayProducts;

    bool isGluonSplitting = false;
    bool isMatrixElement = false;
    bool isAllJets = true;
    if( findGluonSplitting || findMatrixElement)
    {
      bool bFoundS3Quark = false;
      bool bFoundS2Quark = false;
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
      {
	for(std::vector<int>::const_iterator pdgIdIt = jetFlavorPdgId.begin(); pdgIdIt != jetFlavorPdgId.end(); ++pdgIdIt)
	{
	if( abs(it->pdgId()) == abs(*pdgIdIt) && it->status()==3 ) bFoundS3Quark = true;
	if( abs(it->pdgId()) == abs(*pdgIdIt) && it->status()==2 ) bFoundS2Quark = true;
	}
      }
      //if no status 3 quarks but status 2
      if( (!bFoundS3Quark) && bFoundS2Quark) isGluonSplitting = true;
      //if status 2 quark found
      if( bFoundS2Quark ) isMatrixElement = true;
      if( ((!bFoundS3Quark) && bFoundS2Quark) || bFoundS2Quark ) isAllJets = false;
    }
    else
    {
      //isGluonSplitting = true;
      //isMatrixElement = true;
    }

    if( doBosonMatching )
    {
      int bPrimeCount = 0;
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
      {
        if( abs(it->pdgId()) == 7 && it->status() == 3 ) bPrimeCount++;

        if( bPrimeCount>1 && it->status() == 3 ) bPrimeDecayProducts.push_back(&(*it));
      }

      // loop over GenParticles and select bosons isolated from other b' decay products
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
      {
        if( abs(it->pdgId()) == abs(bosonPdgId) && it->status() == 3 )
        {
          h1_BosonPt->Fill( it->pt(), eventWeight );
          h1_BosonEta->Fill( it->eta(), eventWeight );

          bool isIsolated = true;
          if( applyBosonIsolation )
          {
            for(std::vector<const reco::GenParticle*>::const_iterator pIt = bPrimeDecayProducts.begin(); pIt != bPrimeDecayProducts.end(); ++pIt)
            {
              if( &(*it)==(*pIt) ) continue; // skip the boson itself
              bool isBosonDecayProduct = false;
              for(unsigned i=0; i<it->numberOfDaughters(); ++i)
              {
                if( it->daughter(i) == (*pIt) )
                {
                  isBosonDecayProduct = true;
                  break;
                }
              }
              if( isBosonDecayProduct ) continue; // skip the boson decay products

              if( reco::deltaR( it->p4(), (*pIt)->p4() ) < jetRadius ) isIsolated = false;
            }
          }

          if( !isIsolated ) continue;

          h1_BosonPt_Isolated->Fill( it->pt(), eventWeight );
          h1_BosonEta_Isolated->Fill( it->eta(), eventWeight );

          if( doBosonDecayProdSelection )
          {
            bool decayProductsFound = false;

            for(unsigned i=0; i<it->numberOfDaughters(); ++i)
            {
              //std::cout << "Daughter " << i << " PDG ID: " << it->daughter(i)->pdgId() << std::endl;
              for(std::vector<int>::const_iterator pdgIdIt = bosonDecayProdPdgIds.begin(); pdgIdIt != bosonDecayProdPdgIds.end(); ++pdgIdIt)
              {
                if( abs(it->daughter(i)->pdgId()) == abs(*pdgIdIt) )
                {
                  decayProductsFound = true;
                  decayProducts[&(*it)].push_back(it->daughter(i));
                }
              }
            }

            if( decayProductsFound )
            {
              if( decayProducts[&(*it)].size()>2 ) edm::LogError("TooManyDecayProducts") << "More than two boson decay products found.";
              if( decayProducts[&(*it)].size()<2 ) edm::LogError("TooFewDecayProducts") << "Less than two boson decay products found.";

              bosons.push_back(&(*it));
              h1_BosonPt_DecaySel->Fill( it->pt(), eventWeight );
              h1_BosonEta_DecaySel->Fill( it->eta(), eventWeight );
              if( decayProducts[&(*it)].size()>1 )
                h2_BosonPt_dRdecay->Fill( it->pt(), reco::deltaR( decayProducts[&(*it)].at(0)->p4(), decayProducts[&(*it)].at(1)->p4() ), eventWeight );
            }
          }
          else
          {
            bosons.push_back(&(*it));
            h1_BosonPt_DecaySel->Fill( it->pt(), eventWeight );
            h1_BosonEta_DecaySel->Fill( it->eta(), eventWeight );
          }
        }
      }
    }

    // for studying matching efficiency
    for(std::vector<const reco::GenParticle*>::const_iterator bosonIt = bosons.begin(); bosonIt != bosons.end(); ++bosonIt)
    {
      for(PatJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it)
      {
        if( reco::deltaR( (*bosonIt)->p4(), it->p4() ) < bosonMatchingRadius )
        {
          h1_BosonPt_Matched->Fill( (*bosonIt)->pt(), eventWeight );
          break;
        }
      }

      if( decayProducts[*bosonIt].size()>1 )
      {
        for(PatJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it)
        {
          if( reco::deltaR( decayProducts[*bosonIt].at(0)->p4(), it->p4() ) < jetRadius &&
              reco::deltaR( decayProducts[*bosonIt].at(1)->p4(), it->p4() ) < jetRadius )
          {
            h1_BosonPt_DecayProdMatched->Fill( (*bosonIt)->pt(), eventWeight );
            break;
          }
        }
      }
    }

    // loop over jets
    for(PatJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it)
    {
      double jetPt = it->pt();
      // skip the jet if it does not pass pT and eta cuts
      if( !(jetPt > jetPtMin && fabs(it->eta()) < jetAbsEtaMax) ) continue;

      h1_JetPt->Fill(jetPt, eventWeight);
      h1_JetEta->Fill(it->eta(), eventWeight);

      double jetMass = it->mass();

      bool isRightFlavor_gluon = false;
      bool isRightFlavor_matrix = false;
      bool isRightFlavor_all = false;
      if( doJetFlavor )
      {
	if( isGluonSplitting )
	{
	  double jetflavor = it->partonFlavour();
	  for(std::vector<int>::const_iterator pdgIdIt = jetFlavorPdgId.begin(); pdgIdIt != jetFlavorPdgId.end(); ++pdgIdIt)
	  {
	    if( abs(jetflavor) == abs(*pdgIdIt)) isRightFlavor_gluon = true;
	  }	 
	}
	if( isMatrixElement )
	{
	  double jetflavor = it->partonFlavour();
	  for(std::vector<int>::const_iterator pdgIdIt = jetFlavorPdgId.begin(); pdgIdIt != jetFlavorPdgId.end(); ++pdgIdIt)
	  {
	    if( abs(jetflavor) == abs(*pdgIdIt)) isRightFlavor_matrix = true;
	  }	 
	}
	if( isAllJets )
	{
	  double jetflavor = it->partonFlavour();
	  for(std::vector<int>::const_iterator pdgIdIt = jetFlavorPdgId.begin(); pdgIdIt != jetFlavorPdgId.end(); ++pdgIdIt)
	  {
	    if( abs(jetflavor) == abs(*pdgIdIt)) isRightFlavor_all = true;
	  }	 
	}
      }
      else 
      {
	isRightFlavor_all = true;
      }



      PatJetCollection::const_iterator groomedJetMatch;
      bool groomedJetMatchFound = false;
      if( useGroomedJets )
      {
        double dR = jetRadius;
        for(PatJetCollection::const_iterator gjIt = groomedJets->begin(); gjIt != groomedJets->end(); ++gjIt)
        {
          if( reco::deltaR( it->p4(), gjIt->p4() ) < dR )
          {
            groomedJetMatchFound = true;
            dR = reco::deltaR( it->p4(), gjIt->p4() );
            jetMass = gjIt->mass();
            groomedJetMatch = gjIt;
          }
        }
        if( !groomedJetMatchFound ) edm::LogError("NoMatchingGroomedJet") << "Matching groomed jet not found. Using the original jet mass.";
      }

      h2_JetPt_JetPtOverGenJetPt->Fill(jetPt, (it->genJet()!=0 ? jetPt/(it->genJet()->pt()) : -10.), eventWeight);
      h2_JetPt_JetMass->Fill(jetPt, jetMass, eventWeight);

      bool isBosonMatched = false;
      // perform boson matching
      if( doBosonMatching )
      {
        for(std::vector<const reco::GenParticle*>::const_iterator bosonIt = bosons.begin(); bosonIt != bosons.end(); ++bosonIt)
        {
          if( reco::deltaR( (*bosonIt)->p4(), it->p4() ) < bosonMatchingRadius )
          {
            isBosonMatched = true;
            h2_JetPt_JetPtOverBosonPt->Fill( jetPt, jetPt/((*bosonIt)->pt()), eventWeight );
            break;
          }
        }
        for(std::vector<const reco::GenParticle*>::const_iterator bosonIt = bosons.begin(); bosonIt != bosons.end(); ++bosonIt)
        {
          if( decayProducts[*bosonIt].size()>1 )
          {
            if( reco::deltaR( decayProducts[*bosonIt].at(0)->p4(), it->p4() ) < jetRadius &&
                reco::deltaR( decayProducts[*bosonIt].at(1)->p4(), it->p4() ) < jetRadius )
            {
              h1_JetPt_BosonDecayProdMatched->Fill( jetPt, eventWeight );
              if( jetMass > jetMassMin && jetMass < jetMassMax )
                h1_JetPt_BosonDecayProdMatched_JetMass->Fill( jetPt, eventWeight );
              break;
            }
          }
        }
      }
      else
        isBosonMatched = true;

      // skip the jet if it is not matched to a boson or not right flavor
      if( !isBosonMatched || (findGluonSplitting && !isRightFlavor_gluon) || (findMatrixElement && !isRightFlavor_matrix) || ((!findGluonSplitting && !findMatrixElement) && !isRightFlavor_all) ) continue;

      // vector of pointers to subjets
      std::vector<const pat::Jet*> subjets;
      if( useSubJets )
      {
        for(PatJetCollection::const_iterator sjIt = subJets->begin(); sjIt != subJets->end(); ++sjIt)
        {
          if( reco::deltaR( it->p4(), sjIt->p4() ) < jetRadius )
          {
            subjets.push_back(&(*sjIt));
          }
        }

        if( subjets.size()<2 )
          edm::LogError("TooFewSubjets") << "Less than two subjets found.";
        else
        {
          if( subJetMode=="Kt" )
          {
            if( subjets.size()>2 )
            {
              edm::LogError("TooManySubjets") << "More than two subjets found. Will take the two subjets closest to the jet axis.";
              std::sort(subjets.begin(), subjets.end(), ordering_dR(&(*it)));
              subjets.erase(subjets.begin()+2,subjets.end());
              //for(unsigned i=0; i<subjets.size(); ++i)
                //std::cout << "dR(jet,subjet) for subjet" << i << ": " << reco::deltaR( it->p4(), subjets.at(i)->p4() ) << std::endl;
            }
          }
          else if( subJetMode=="Filtered" )
          {
            std::sort(subjets.begin(), subjets.end(), ordering_Pt("Uncorrected"));
            //for(unsigned i=0; i<subjets.size(); ++i)
              //std::cout << "Uncorrected Pt for subjet" << i << ": " << subjets.at(i)->correctedJet("Uncorrected").pt() << std::endl;
          }
          else
            edm::LogError("IllegalSubJetMode") << "Allowed subjet modes are Kt and Filtered.";
        }
      }

      h1_JetPt_BosonMatched->Fill(jetPt, eventWeight);
      h1_JetEta_BosonMatched->Fill(it->eta(), eventWeight);
      h2_JetPt_JetPtOverGenJetPt_BosonMatched->Fill(jetPt, (it->genJet()!=0 ? jetPt/(it->genJet()->pt()) : -10.), eventWeight);
      h2_JetPt_JetMass_BosonMatched->Fill(jetPt, jetMass, eventWeight);

      // fill nPV_JetMass histograms
      std::string suffix = Form("%.0ftoInf",jetPtMin);
      h2_nPV_JetMass_Pt[suffix]->Fill(nPV, jetMass, eventWeight);
      for(unsigned i=0; i<jetPtBins; ++i)
      {
        if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
        {
          suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
          h2_nPV_JetMass_Pt[suffix]->Fill(nPV, jetMass, eventWeight);
        }
      }
      if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
      {
        suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
        h2_nPV_JetMass_Pt[suffix]->Fill(nPV, jetMass, eventWeight);
      }

      // skip the jet if it does not pass the invariant mass cut
      if( !(jetMass > jetMassMin && jetMass < jetMassMax) ) continue;

      h1_JetPt_BosonMatched_JetMass->Fill(jetPt, eventWeight);
      h1_JetEta_BosonMatched_JetMass->Fill(it->eta(), eventWeight);

      // get b-tag discriminators
      //double jet_CSV_discr = it->bDiscriminator("combinedSecondaryVertexBJetTags");
      double jet_CSV_discr = it->bDiscriminator((bdiscriminator).c_str());
      double subJet1_CSV_discr = -999., subJet2_CSV_discr = -999.;
      if( subjets.size()>1 )
      {
//      subJet1_CSV_discr = subjets.at(0)->bDiscriminator("combinedSecondaryVertexBJetTags");
//      subJet2_CSV_discr = subjets.at(1)->bDiscriminator("combinedSecondaryVertexBJetTags");
        subJet1_CSV_discr = subjets.at(0)->bDiscriminator((bdiscriminator).c_str());
        subJet2_CSV_discr = subjets.at(1)->bDiscriminator((bdiscriminator).c_str());
      }
      double jet_DoubleB_discr = it->bDiscriminator("doubleSecondaryVertexHighEffBJetTags");

      h1_JetCSVDiscr_BosonMatched_JetMass->Fill( jet_CSV_discr, eventWeight);
      h1_SubJetCSVDiscr_BosonMatched_JetMass->Fill( subJet1_CSV_discr, eventWeight);
      h1_SubJetCSVDiscr_BosonMatched_JetMass->Fill( subJet2_CSV_discr, eventWeight);
      h1_JetDoubleBDiscr_BosonMatched_JetMass->Fill( jet_DoubleB_discr, eventWeight);

      if( jet_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_CSVL->Fill(jetPt, eventWeight);
      if( jet_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_CSVM->Fill(jetPt, eventWeight);
      if( subJet1_CSV_discr>0.244 && subJet2_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_SubJetCSVL->Fill(jetPt, eventWeight);
      if( subJet1_CSV_discr>0.679 && subJet2_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_SubJetCSVM->Fill(jetPt, eventWeight);
      if( jet_DoubleB_discr>0. ) h1_JetPt_BosonMatched_JetMass_DoubleB->Fill(jetPt, eventWeight);

      PatJetCollection::const_iterator substructJet = it;
      if( useGroomedJetSubstr && groomedJetMatchFound ) substructJet = groomedJetMatch;
      if( !useMassDrop )
      {
        std::vector<fastjet::PseudoJet> fjConstituents;
        std::vector<edm::Ptr<reco::PFCandidate> > constituents = substructJet->getPFConstituents();
        std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator m;
        for ( m = constituents.begin(); m != constituents.end(); ++m )
        {
          reco::PFCandidatePtr constit = *m;
          if (constit->pt() == 0)
          {
            edm::LogError("NullTransverseMomentum") << "dropping input candidate with pt=0";
            continue;
          }
          fjConstituents.push_back(fastjet::PseudoJet(constit->px(),constit->py(),constit->pz(),constit->energy()));
          fjConstituents.back().set_user_index(m - constituents.begin());
        }

        double tau1 = nsubjettinessCalculator.getTau(1,fjConstituents);
        double tau2 = nsubjettinessCalculator.getTau(2,fjConstituents);
        double tau2overtau1 = (tau1>0 ? tau2/tau1 : -10.);

        // fill nPV_tau histograms
        suffix = Form("%.0ftoInf",jetPtMin);
        h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
        h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
        h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);
        for(unsigned i=0; i<jetPtBins; ++i)
        {
          if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
          {
            suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
            h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
            h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
            h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);
          }
        }
        if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
        {
          suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
          h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
          h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
          h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);
        }

        if( tau2/tau1<nsubjCut )
        {
          if( jet_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_Nsubj_CSVL->Fill(jetPt, eventWeight);
          if( jet_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_Nsubj_CSVM->Fill(jetPt, eventWeight);
          if( subJet1_CSV_discr>0.244 && subJet2_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVL->Fill(jetPt, eventWeight);
          if( subJet1_CSV_discr>0.679 && subJet2_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVM->Fill(jetPt, eventWeight);
          if( jet_DoubleB_discr>0. ) h1_JetPt_BosonMatched_JetMass_Nsubj_DoubleB->Fill(jetPt, eventWeight);
        }
      }
      else
      {
        double fatJetMass = jetMass;
        if( useUncorrMassForMassDrop ) fatJetMass = substructJet->correctedJet("Uncorrected").mass();
        double massDrop = ( (jetMass>0. && substructJet->numberOfDaughters()>1) ? std::max( substructJet->daughter(0)->mass(), substructJet->daughter(1)->mass() ) / fatJetMass : -10.);
        // fill nPV_MassDrop histograms
        suffix = Form("%.0ftoInf",jetPtMin);
        h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
        for(unsigned i=0; i<jetPtBins; ++i)
        {
          if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
          {
            suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
            h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
          }
        }
        if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
        {
          suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
          h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
        }
      }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
RutgersJetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
RutgersJetAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
RutgersJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
RutgersJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
RutgersJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
RutgersJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RutgersJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RutgersJetAnalyzer);
