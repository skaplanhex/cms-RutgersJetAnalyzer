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
// $Id: RutgersJetAnalyzer.cc,v 1.16 2013/03/15 04:19:41 ferencek Exp $
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
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
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
    const bool             useEventWeight;
    const edm::InputTag    genParticleTag;
    const edm::InputTag    jetsTag;
    const bool             useGroomedJets;
    const edm::InputTag    groomedJetsTag;
    const bool             useSubJets;
    const edm::InputTag    groomedBasicJetsTag;
    const std::string      subJetMode;
    const edm::InputTag    pvTag;
    const double           jetRadius;          // radius for jet clustering
    const bool             doBosonMatching;    // parameter for deciding if matching is on or off
    const double           bosonMatchingRadius;
    const int              bosonPdgId;
    const bool             applyBosonIsolation;
    const bool             doBosonDecayProdSelection;
    const std::vector<int> bosonDecayProdPdgIds;
    const bool             useMassDrop;
    const double           jetPtMin;
    const unsigned         jetPtBins;
    const double           jetPtBinWidth;
    const double           jetAbsEtaMax;
    const double           jetMassMin;
    const double           jetMassMax;
    const double           nsubjCut;
    const bool             useOnePassKtAxes;
    const bool             useGroomedJetSubstr;
    const bool             useUncorrMassForMassDrop;
    const std::string      bdiscriminator;
    const bool             doJetFlavor;
    const std::vector<int> jetFlavorPdgIds;
    const bool             useAltGSPbDef;
    const bool             findGluonSplitting;
    const bool             findMatrixElement;

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
    TH2D *h2_JetPt_dRsubjets_BosonMatched;
    TH2D *h2_JetPt_dRsubjets_BosonMatched_JetMass;

    TH2D *h2_JetPt_mindRSubjet1Bhadron_BosonMatched;
    TH2D *h2_JetPt_mindRSubjet2Bhadron_BosonMatched;
    TH2D *h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass;

    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched;
    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass;

    TH1D *h1_JetCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_SubJetMinCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_JetDoubleBDiscr_BosonMatched_JetMass;

    TH2D *h2_JetPt_JetCSVL_BosonMatched_JetMass;
    TH2D *h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass;

    std::map<std::string, TH2D*> h2_nPV_JetMass_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_MassDrop_Pt;

    std::map<std::string, TH2D*> h2_JetMass_nTracks_Pt;
    std::map<std::string, TH2D*> h2_JetMass_nSelectedTracks_Pt;
    std::map<std::string, TH2D*> h2_JetMass_tau2tau1_Pt;
    std::map<std::string, TH2D*> h2_JetMass_SubJetMinCSVL_Pt;

    std::map<std::string, TH2D*> h2_nTracks_tau2tau1_Pt;
    std::map<std::string, TH2D*> h2_nSelectedTracks_tau2tau1_Pt;

    std::map<std::string, TH2D*> h2_nTracks_SubJetMinCSVL_Pt;
    std::map<std::string, TH2D*> h2_nSelectedTracks_SubJetMinCSVL_Pt;
    std::map<std::string, TH2D*> h2_tau2tau1_SubJetMinCSVL_Pt;
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
  groomedBasicJetsTag(iConfig.getParameter<edm::InputTag>("GroomedBasicJetsTag")),
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
  useOnePassKtAxes( iConfig.exists("UseOnePassKtAxes") ? iConfig.getParameter<bool>("UseOnePassKtAxes") : true ),
  useGroomedJetSubstr( iConfig.exists("UseGroomedJetSubstructure") ? iConfig.getParameter<bool>("UseGroomedJetSubstructure") : false ),
  useUncorrMassForMassDrop( iConfig.exists("UseUncorrMassForMassDrop") ? iConfig.getParameter<bool>("UseUncorrMassForMassDrop") : true ),
  bdiscriminator(iConfig.getParameter<std::string>("Bdiscriminator")),
  doJetFlavor(iConfig.getParameter<bool>("DoJetFlavor")),
  jetFlavorPdgIds(iConfig.getParameter<std::vector<int> >("JetFlavorPdgIds")),
  useAltGSPbDef( iConfig.exists("UseAltGSPbDef") ? iConfig.getParameter<bool>("UseAltGSPbDef") : false ),
  findGluonSplitting( iConfig.exists("FindGluonSplitting") ? iConfig.getParameter<bool>("FindGluonSplitting") : false ),
  findMatrixElement( iConfig.exists("FindMatrixElement") ? iConfig.getParameter<bool>("FindMatrixElement") : false ),
  nsubjettinessCalculator(( useOnePassKtAxes ? Njettiness::onepass_kt_axes : Njettiness::kt_axes ), NsubParameters(1.0, jetRadius, jetRadius))

{
    //now do what ever initialization is needed
    int pvBins=50;
    double pvMin=-0.5, pvMax=49.5;
    int trackBins=200;
    double trackMin=-0.5, trackMax=199.5;
    int ptBins=250;
    double ptMin=0., ptMax=1000.;
    int dRBins=100;
    double dRMin=0., dRMax=5.;
    int etaBins=160;
    double etaMin=-4., etaMax=4.;
    int massBins=200;
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

    h2_JetPt_JetPtOverBosonPt = fs->make<TH2D>("h2_JetPt_JetPtOverBosonPt",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{boson}",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetPtOverGenJetPt = fs->make<TH2D>("h2_JetPt_JetPtOverGenJetPt",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{genjet}",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetPtOverGenJetPt_BosonMatched = fs->make<TH2D>("h2_JetPt_JetPtOverGenJetPt_BosonMatched",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{genjet}",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetMass = fs->make<TH2D>("h2_JetPt_JetMass",";p_{T} [GeV];m_{jet} [GeV]",ptBins,ptMin,ptMax,massBins,massMin,massMax);
    h2_JetPt_JetMass_BosonMatched = fs->make<TH2D>("h2_JetPt_JetMass_BosonMatched",";p_{T} [GeV];m_{jet} [GeV]",ptBins,ptMin,ptMax,massBins,massMin,massMax);
    h2_JetPt_dRsubjets_BosonMatched = fs->make<TH2D>("h2_JetPt_dRsubjets_BosonMatched",";p_{T} [GeV];#DeltaR(subjet_{1},subjet_{2})",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);
    h2_JetPt_dRsubjets_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_dRsubjets_BosonMatched_JetMass",";p_{T} [GeV];#DeltaR(subjet_{1},subjet_{2})",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);

    h2_JetPt_mindRSubjet1Bhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_mindRSubjet1Bhadron_BosonMatched",";p_{T} [GeV];min #DeltaR(subjet_{1},B hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRSubjet2Bhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_mindRSubjet2Bhadron_BosonMatched",";p_{T} [GeV];min #DeltaR(subjet_{2},B hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass",";p_{T} [GeV];min #DeltaR(subjet_{1},B hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass",";p_{T} [GeV];min #DeltaR(subjet_{2},B hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);

    h2_JetPt_SameMatchedBhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched",";p_{T} [GeV];Same matched B hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass",";p_{T} [GeV];Same matched B hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);

    h1_JetCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetCSVDiscr_BosonMatched_JetMass",";Jet CSV Discr;",100,0.,1.);
    h1_SubJetMinCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_SubJetMinCSVDiscr_BosonMatched_JetMass",";SubJet min CSV Discr;",100,0.,1.);
    h1_JetDoubleBDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetDoubleBDiscr_BosonMatched_JetMass",";Jet DoubleB Discr;",100,0.,10.);

    h2_JetPt_JetCSVL_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_JetCSVL_BosonMatched_JetMass",";p_{T} [GeV];Jet CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);
    h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",";p_{T} [GeV];SubJet min CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);

    for(unsigned i=0; i<=(jetPtBins+1); ++i)
    {
      std::string suffix, title;

      if(i==0)
      {
        suffix = Form("%.0ftoInf",(jetPtMin + jetPtBinWidth*i));
        title = Form("p_{T}>%.0f GeV",(jetPtMin + jetPtBinWidth*i));

        h2_nPV_JetMass_Pt[suffix]         = fs->make<TH2D>(("h2_nPV_JetMass_Pt" + suffix).c_str(),(title + ";nPV;m_{jet} [GeV]").c_str(),pvBins,pvMin,pvMax,massBins,massMin,massMax);
        h2_nPV_tau1_Pt[suffix]            = fs->make<TH2D>(("h2_nPV_tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2_Pt[suffix]            = fs->make<TH2D>(("h2_nPV_tau2_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_nPV_tau2tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}/#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_nPV_MassDrop_Pt" + suffix).c_str(),(title + ";nPV;m_{subjet1}/m_{jet}").c_str(),pvBins,pvMin,pvMax,massDropBins,massDropMin,massDropMax);

        h2_JetMass_nTracks_Pt[suffix]         = fs->make<TH2D>(("h2_JetMass_nTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_nSelectedTracks_Pt[suffix] = fs->make<TH2D>(("h2_JetMass_nSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nSelectedTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_tau2tau1_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#tau_{2}/#tau_{1}").c_str(),massBins,massMin,massMax,tauBins,tauMin,tauMax);
        h2_JetMass_SubJetMinCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);

        h2_nTracks_tau2tau1_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);
        h2_nSelectedTracks_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nSelectedTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);

        h2_nTracks_SubJetMinCSVL_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";nTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        h2_nSelectedTracks_SubJetMinCSVL_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";nSelectedTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        h2_tau2tau1_SubJetMinCSVL_Pt[suffix]        = fs->make<TH2D>(("h2_tau2tau1_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";#tau_{2}/#tau_{1};SubJet min CSV Discr").c_str(),tauBins,tauMin,tauMax,100,0.,1.);
      }
      else if(i==(jetPtBins+1))
      {
        suffix = Form("%.0ftoInf",(jetPtMin + jetPtBinWidth*(i-1)));
        title = Form("p_{T}>%.0f GeV",(jetPtMin + jetPtBinWidth*(i-1)));

        h2_nPV_JetMass_Pt[suffix]         = fs->make<TH2D>(("h2_nPV_JetMass_Pt" + suffix).c_str(),(title + ";nPV;m_{jet} [GeV]").c_str(),pvBins,pvMin,pvMax,massBins,massMin,massMax);
        h2_nPV_tau1_Pt[suffix]            = fs->make<TH2D>(("h2_nPV_tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2_Pt[suffix]            = fs->make<TH2D>(("h2_nPV_tau2_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_nPV_tau2tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}/#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_nPV_MassDrop_Pt" + suffix).c_str(),(title + ";nPV;m_{subjet1}/m_{jet}").c_str(),pvBins,pvMin,pvMax,massDropBins,massDropMin,massDropMax);

        h2_JetMass_nTracks_Pt[suffix]         = fs->make<TH2D>(("h2_JetMass_nTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_nSelectedTracks_Pt[suffix] = fs->make<TH2D>(("h2_JetMass_nSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nSelectedTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_tau2tau1_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#tau_{2}/#tau_{1}").c_str(),massBins,massMin,massMax,tauBins,tauMin,tauMax);
        h2_JetMass_SubJetMinCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);

        h2_nTracks_tau2tau1_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);
        h2_nSelectedTracks_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nSelectedTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);

        h2_nTracks_SubJetMinCSVL_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";nTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        h2_nSelectedTracks_SubJetMinCSVL_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";nSelectedTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        h2_tau2tau1_SubJetMinCSVL_Pt[suffix]        = fs->make<TH2D>(("h2_tau2tau1_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";#tau_{2}/#tau_{1};SubJet min CSV Discr").c_str(),tauBins,tauMin,tauMax,100,0.,1.);
      }
      else
      {
        suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*(i-1)),(jetPtMin + jetPtBinWidth*i));
        title = Form("%.0f<p_{T}<%.0f GeV",(jetPtMin + jetPtBinWidth*(i-1)),(jetPtMin + jetPtBinWidth*i));

        h2_nPV_JetMass_Pt[suffix]         = fs->make<TH2D>(("h2_nPV_JetMass_Pt" + suffix).c_str(),(title + ";nPV;m_{jet} [GeV]").c_str(),pvBins,pvMin,pvMax,massBins,massMin,massMax);
        h2_nPV_tau1_Pt[suffix]            = fs->make<TH2D>(("h2_nPV_tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2_Pt[suffix]            = fs->make<TH2D>(("h2_nPV_tau2_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_nPV_tau2tau1_Pt" + suffix).c_str(),(title + ";nPV;#tau_{2}/#tau_{1}").c_str(),pvBins,pvMin,pvMax,tauBins,tauMin,tauMax);
        h2_nPV_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_nPV_MassDrop_Pt" + suffix).c_str(),(title + ";nPV;m_{subjet1}/m_{jet}").c_str(),pvBins,pvMin,pvMax,massDropBins,massDropMin,massDropMax);

        h2_JetMass_nTracks_Pt[suffix]         = fs->make<TH2D>(("h2_JetMass_nTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_nSelectedTracks_Pt[suffix] = fs->make<TH2D>(("h2_JetMass_nSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nSelectedTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_tau2tau1_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#tau_{2}/#tau_{1}").c_str(),massBins,massMin,massMax,tauBins,tauMin,tauMax);
        h2_JetMass_SubJetMinCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);

        h2_nTracks_tau2tau1_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);
        h2_nSelectedTracks_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nSelectedTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);

        h2_nTracks_SubJetMinCSVL_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";nTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        h2_nSelectedTracks_SubJetMinCSVL_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";nSelectedTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        h2_tau2tau1_SubJetMinCSVL_Pt[suffix]        = fs->make<TH2D>(("h2_tau2tau1_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";#tau_{2}/#tau_{1};SubJet min CSV Discr").c_str(),tauBins,tauMin,tauMax,100,0.,1.);
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

    edm::Handle<PatJetCollection> groomedBasicJets;
    if( useSubJets )
    {
      iEvent.getByLabel(groomedBasicJetsTag,groomedBasicJets);
    }

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
    if( findGluonSplitting || findMatrixElement )
    {
      bool bFoundS3Quark = false;
      bool bFoundS2Quark = false;
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
      {
        if( abs(it->pdgId()) == 5 && it->status()==3 ) bFoundS3Quark = true;
        if( abs(it->pdgId()) == 5 && it->status()==2 ) bFoundS2Quark = true;
      }
      // if no status 3 b quark but status 2
      if( (!bFoundS3Quark) && bFoundS2Quark) isGluonSplitting = true;
      // if status 3 b quark found
      if( bFoundS3Quark ) isMatrixElement = true;
    }

    // skip event if it does not pass the selection
    if( (findGluonSplitting && !isGluonSplitting) || (findMatrixElement && !isMatrixElement) ) return;


    if( doBosonMatching )
    {
      int bPrimeCount = 0;
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
      {
        // count b' in the list of GenParticles
        if( abs(it->pdgId()) == 7 && it->status() == 3 ) bPrimeCount++;
        // only take status=3 GenParticles that appear after b' and b'bar have appeared in the list
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


      bool isRightFlavor = false;
      // check jet flavor
      if( doJetFlavor )
      {
        int jetFlavor = it->partonFlavour();

        if( useAltGSPbDef ) // use alternative gluon splitting b-jet definition based on the number of B hadrons inside the jet cone
        {
          int nMatchedBHadrons = 0;
          for(reco::GenParticleCollection::const_iterator gpIt = genParticles->begin(); gpIt != genParticles->end(); ++gpIt)
          {
            int id = abs(gpIt->pdgId());
            // skip GenParticle if not B hadron
            if ( !((id/100)%10 == 5 || (id/1000)%10 == 5) ) continue;

            if( reco::deltaR( it->p4(), gpIt->p4() ) < jetRadius ) ++nMatchedBHadrons;
          }
          //std::cout << "nMatchedBHadrons: " << nMatchedBHadrons << std::endl;
          if( nMatchedBHadrons>=2 ) jetFlavor = 85; // custom jet flavor code for gluon splitting b jets
        }

        for(std::vector<int>::const_iterator pdgIdIt = jetFlavorPdgIds.begin(); pdgIdIt != jetFlavorPdgIds.end(); ++pdgIdIt)
        {
          if( abs(jetFlavor) == abs(*pdgIdIt))
          {
            isRightFlavor = true;
            break;
          }
        }
      }
      else
        isRightFlavor = true;

      // skip the jet if it does not have the right flavor
      if( !isRightFlavor ) continue;


      h1_JetPt->Fill(jetPt, eventWeight);
      h1_JetEta->Fill(it->eta(), eventWeight);


      double jetMass = it->mass();
      PatJetCollection::const_iterator groomedJetMatch;
      bool groomedJetMatchFound = false;
      if( useGroomedJets )
      {
        double dR = jetRadius;
        for(PatJetCollection::const_iterator gjIt = groomedJets->begin(); gjIt != groomedJets->end(); ++gjIt)
        {
          double dR_temp = reco::deltaR( it->p4(), gjIt->p4() );
          if( dR_temp < dR )
          {
            groomedJetMatchFound = true;
            dR = dR_temp;
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

      // skip the jet if it is not matched to a boson
      if( !isBosonMatched ) continue;


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


      // vector of pointers to subjets
      std::vector<const pat::Jet*> subjets;
      if( useSubJets )
      {
        double dR = jetRadius;
        PatJetCollection::const_iterator groomedBasicJetMatch;
        for(PatJetCollection::const_iterator gbjIt = groomedBasicJets->begin(); gbjIt != groomedBasicJets->end(); ++gbjIt)
        {
          double dR_temp = reco::deltaR( it->p4(), gbjIt->p4() );
          if( dR_temp < dR )
          {
            dR = dR_temp;
            groomedBasicJetMatch = gbjIt;
          }
        }
        //std::cout << "number of subjets: " << groomedBasicJetMatch->numberOfDaughters() << std::endl;
        //std::cout << "jet pt: " << it->correctedJet("Uncorrected").p4().pt() << " eta=" << it->correctedJet("Uncorrected").p4().eta() << " phi=" << it->correctedJet("Uncorrected").p4().phi() << " nd=" << it->numberOfDaughters() << std::endl;
        //std::cout << "groomed jet pt: " << groomedJetMatch->correctedJet("Uncorrected").p4().pt() << " eta=" << groomedJetMatch->correctedJet("Uncorrected").p4().eta() << " phi=" << groomedJetMatch->correctedJet("Uncorrected").p4().phi() << " nd=" << groomedJetMatch->numberOfDaughters() << std::endl;
        //std::cout << "groomed basic jet pt: " << groomedBasicJetMatch->p4().pt() << " eta=" << groomedBasicJetMatch->p4().eta() << " phi=" << groomedBasicJetMatch->p4().phi() << std::endl;
        for(unsigned d=0; d<groomedBasicJetMatch->numberOfDaughters(); ++d)
        {
          //std::cout << "subjet " << d << ": pt=" << groomedBasicJetMatch->daughter(d)->p4().pt()  << " eta=" << groomedBasicJetMatch->daughter(d)->p4().eta()  << " phi=" << groomedBasicJetMatch->daughter(d)->p4().phi() << std::endl;
          const reco::Candidate *subjet =  groomedBasicJetMatch->daughter(d);
          const pat::Jet *patsubjet = dynamic_cast<const pat::Jet*>(subjet);
          subjets.push_back(patsubjet);
        }


        if( subjets.size()<2 )
          edm::LogError("TooFewSubjets") << "Less than two subjets found.";
        else
        {
          if( subJetMode=="Kt" || subJetMode=="Pruned" )
          {
            if( subjets.size()>2 )
            {
              edm::LogError("TooManySubjets") << "More than two subjets found. Will take the two subjets closest to the jet axis.";
              std::sort(subjets.begin(), subjets.end(), ordering_dR(&(*it)));
              subjets.erase(subjets.begin()+2,subjets.end());
              //for(unsigned i=0; i<subjets.size(); ++i)
                //std::cout << "dR(jet,subjet) for subjet" << i << ": " << reco::deltaR( it->p4(), subjets.at(i)->p4() ) << std::endl;
            }
            // sort subjets by uncorrected Pt
            std::sort(subjets.begin(), subjets.end(), ordering_Pt("Uncorrected"));
            //for(unsigned i=0; i<subjets.size(); ++i)
              //std::cout << "Uncorrected Pt for subjet" << i << ": " << subjets.at(i)->correctedJet("Uncorrected").pt() << std::endl;
          }
          else if( subJetMode=="Filtered" )
          {
            // sort subjets by uncorrected Pt
            std::sort(subjets.begin(), subjets.end(), ordering_Pt("Uncorrected"));
            //for(unsigned i=0; i<subjets.size(); ++i)
              //std::cout << "Uncorrected Pt for subjet" << i << ": " << subjets.at(i)->correctedJet("Uncorrected").pt() << std::endl;
          }
          else
            edm::LogError("IllegalSubJetMode") << "Allowed subjet modes are Kt, Pruned, and Filtered.";
        }
      }


      double mindRsubjet1 = 999., mindRsubjet2 = 999.;
      reco::GenParticleCollection::const_iterator bHadronMatchSubjet1 = genParticles->end(), bHadronMatchSubjet2 = genParticles->end();
      // find the closest B hadron for each of the two subjets
      if( subjets.size()>1 )
      {
        for(reco::GenParticleCollection::const_iterator gpIt = genParticles->begin(); gpIt != genParticles->end(); ++gpIt)
        {
          int id = abs(gpIt->pdgId());
          // skip GenParticle if not B hadron
          if ( !((id/100)%10 == 5 || (id/1000)%10 == 5) ) continue;

          double dRsubjet1 = reco::deltaR( subjets.at(0)->p4(), gpIt->p4() );
          double dRsubjet2 = reco::deltaR( subjets.at(1)->p4(), gpIt->p4() );
          if( dRsubjet1 < mindRsubjet1 )
          {
            mindRsubjet1 = dRsubjet1;
            bHadronMatchSubjet1 = gpIt;
          }
          if( dRsubjet2 < mindRsubjet2 )
          {
            mindRsubjet2 = dRsubjet2;
            bHadronMatchSubjet2 = gpIt;
          }
        }
      }

      h2_JetPt_mindRSubjet1Bhadron_BosonMatched->Fill(jetPt, (mindRsubjet1<999. ? mindRsubjet1 : -99.), eventWeight);
      h2_JetPt_mindRSubjet2Bhadron_BosonMatched->Fill(jetPt, (mindRsubjet2<999. ? mindRsubjet2 : -99.), eventWeight);


      if( subjets.size()>1 )
      {
        h2_JetPt_dRsubjets_BosonMatched->Fill(jetPt, reco::deltaR( subjets.at(0)->p4(), subjets.at(1)->p4() ), eventWeight);
        if( bHadronMatchSubjet1 != genParticles->end() && bHadronMatchSubjet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadron_BosonMatched->Fill(jetPt, ( bHadronMatchSubjet1 == bHadronMatchSubjet2 ? 1. : 0. ), eventWeight);
      }


      // skip the jet if it does not pass the invariant mass cut
      if( !(jetMass > jetMassMin && jetMass < jetMassMax) ) continue;


      h1_JetPt_BosonMatched_JetMass->Fill(jetPt, eventWeight);
      h1_JetEta_BosonMatched_JetMass->Fill(it->eta(), eventWeight);

      h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass->Fill(jetPt, (mindRsubjet1<999. ? mindRsubjet1 : -99.), eventWeight);
      h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass->Fill(jetPt, (mindRsubjet2<999. ? mindRsubjet2 : -99.), eventWeight);

      if( subjets.size()>1 )
      {
        h2_JetPt_dRsubjets_BosonMatched_JetMass->Fill(jetPt, reco::deltaR( subjets.at(0)->p4(), subjets.at(1)->p4() ), eventWeight);
        if( bHadronMatchSubjet1 != genParticles->end() && bHadronMatchSubjet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass->Fill(jetPt, ( bHadronMatchSubjet1 == bHadronMatchSubjet2 ? 1. : 0. ), eventWeight);
      }


      // get b-tag discriminators
      double jet_CSV_discr = it->bDiscriminator((bdiscriminator).c_str());
      double subJet1_CSV_discr = -999., subJet2_CSV_discr = -999.;
      if( subjets.size()>1 )
      {
        subJet1_CSV_discr = subjets.at(0)->bDiscriminator(bdiscriminator.c_str());
        subJet2_CSV_discr = subjets.at(1)->bDiscriminator(bdiscriminator.c_str());
      }
      double minSubJet_CSV_discr = std::min(subJet1_CSV_discr, subJet2_CSV_discr);
      double jet_DoubleB_discr = it->bDiscriminator("doubleSecondaryVertexHighEffBJetTags");

      h1_JetCSVDiscr_BosonMatched_JetMass->Fill( jet_CSV_discr, eventWeight);
      h1_SubJetMinCSVDiscr_BosonMatched_JetMass->Fill( minSubJet_CSV_discr, eventWeight);
      h1_JetDoubleBDiscr_BosonMatched_JetMass->Fill( jet_DoubleB_discr, eventWeight);

      h2_JetPt_JetCSVL_BosonMatched_JetMass->Fill(jetPt, jet_CSV_discr, eventWeight);;
      h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass->Fill(jetPt, minSubJet_CSV_discr, eventWeight);;

      if( jet_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_CSVL->Fill(jetPt, eventWeight);
      if( jet_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_CSVM->Fill(jetPt, eventWeight);
      if( subJet1_CSV_discr>0.244 && subJet2_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_SubJetCSVL->Fill(jetPt, eventWeight);
      if( subJet1_CSV_discr>0.679 && subJet2_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_SubJetCSVM->Fill(jetPt, eventWeight);
      if( jet_DoubleB_discr>0. ) h1_JetPt_BosonMatched_JetMass_DoubleB->Fill(jetPt, eventWeight);

      PatJetCollection::const_iterator substructJet = it;
      if( useGroomedJetSubstr && groomedJetMatchFound ) substructJet = groomedJetMatch;
      // N-subjettiness
      std::vector<fastjet::PseudoJet> fjConstituents;
      std::vector<edm::Ptr<reco::PFCandidate> > constituents = substructJet->getPFConstituents();
      std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator m;
      for ( m = constituents.begin(); m != constituents.end(); ++m )
      {
        reco::PFCandidatePtr constit = *m;
        if (constit->pt() == 0)
        {
          edm::LogWarning("NullTransverseMomentum") << "dropping input candidate with pt=0";
          continue;
        }
        fjConstituents.push_back(fastjet::PseudoJet(constit->px(),constit->py(),constit->pz(),constit->energy()));
        fjConstituents.back().set_user_index(m - constituents.begin());
      }

      double tau1 = nsubjettinessCalculator.getTau(1,fjConstituents);
      double tau2 = nsubjettinessCalculator.getTau(2,fjConstituents);
      double tau2overtau1 = (tau1>0 ? tau2/tau1 : -10.);

      double nTracks = it->associatedTracks().size();
      double nSelectedTracks = (it->hasTagInfo("impactParameter") ? it->tagInfoTrackIP("impactParameter")->selectedTracks().size() : -99.);

      // fill various 2D histograms
      suffix = Form("%.0ftoInf",jetPtMin);
      h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
      h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
      h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);

      h2_JetMass_nTracks_Pt[suffix]->Fill(jetMass, nTracks, eventWeight);
      h2_JetMass_nSelectedTracks_Pt[suffix]->Fill(jetMass, nSelectedTracks, eventWeight);
      h2_JetMass_tau2tau1_Pt[suffix]->Fill(jetMass, tau2overtau1, eventWeight);
      h2_JetMass_SubJetMinCSVL_Pt[suffix]->Fill(jetMass, minSubJet_CSV_discr, eventWeight);

      h2_nTracks_tau2tau1_Pt[suffix]->Fill(nTracks, tau2overtau1, eventWeight);
      h2_nSelectedTracks_tau2tau1_Pt[suffix]->Fill(nSelectedTracks, tau2overtau1, eventWeight);

      h2_nTracks_SubJetMinCSVL_Pt[suffix]->Fill(nTracks, minSubJet_CSV_discr, eventWeight);
      h2_nSelectedTracks_SubJetMinCSVL_Pt[suffix]->Fill(nSelectedTracks, minSubJet_CSV_discr, eventWeight);
      h2_tau2tau1_SubJetMinCSVL_Pt[suffix]->Fill(tau2overtau1, minSubJet_CSV_discr, eventWeight);
      for(unsigned i=0; i<jetPtBins; ++i)
      {
        if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
        {
          suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
          h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
          h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
          h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);

          h2_JetMass_nTracks_Pt[suffix]->Fill(jetMass, nTracks, eventWeight);
          h2_JetMass_nSelectedTracks_Pt[suffix]->Fill(jetMass, nSelectedTracks, eventWeight);
          h2_JetMass_tau2tau1_Pt[suffix]->Fill(jetMass, tau2overtau1, eventWeight);
          h2_JetMass_SubJetMinCSVL_Pt[suffix]->Fill(jetMass, minSubJet_CSV_discr, eventWeight);

          h2_nTracks_tau2tau1_Pt[suffix]->Fill(nTracks, tau2overtau1, eventWeight);
          h2_nSelectedTracks_tau2tau1_Pt[suffix]->Fill(nSelectedTracks, tau2overtau1, eventWeight);

          h2_nTracks_SubJetMinCSVL_Pt[suffix]->Fill(nTracks, minSubJet_CSV_discr, eventWeight);
          h2_nSelectedTracks_SubJetMinCSVL_Pt[suffix]->Fill(nSelectedTracks, minSubJet_CSV_discr, eventWeight);
          h2_tau2tau1_SubJetMinCSVL_Pt[suffix]->Fill(tau2overtau1, minSubJet_CSV_discr, eventWeight);
        }
      }
      if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
      {
        suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
        h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
        h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
        h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);

        h2_JetMass_nTracks_Pt[suffix]->Fill(jetMass, nTracks, eventWeight);
        h2_JetMass_nSelectedTracks_Pt[suffix]->Fill(jetMass, nSelectedTracks, eventWeight);
        h2_JetMass_tau2tau1_Pt[suffix]->Fill(jetMass, tau2overtau1, eventWeight);
        h2_JetMass_SubJetMinCSVL_Pt[suffix]->Fill(jetMass, minSubJet_CSV_discr, eventWeight);

        h2_nTracks_tau2tau1_Pt[suffix]->Fill(nTracks, tau2overtau1, eventWeight);
        h2_nSelectedTracks_tau2tau1_Pt[suffix]->Fill(nSelectedTracks, tau2overtau1, eventWeight);

        h2_nTracks_SubJetMinCSVL_Pt[suffix]->Fill(nTracks, minSubJet_CSV_discr, eventWeight);
        h2_nSelectedTracks_SubJetMinCSVL_Pt[suffix]->Fill(nSelectedTracks, minSubJet_CSV_discr, eventWeight);
        h2_tau2tau1_SubJetMinCSVL_Pt[suffix]->Fill(tau2overtau1, minSubJet_CSV_discr, eventWeight);
      }

      if( tau2/tau1<nsubjCut )
      {
        if( jet_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_Nsubj_CSVL->Fill(jetPt, eventWeight);
        if( jet_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_Nsubj_CSVM->Fill(jetPt, eventWeight);
        if( subJet1_CSV_discr>0.244 && subJet2_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVL->Fill(jetPt, eventWeight);
        if( subJet1_CSV_discr>0.679 && subJet2_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_Nsubj_SubJetCSVM->Fill(jetPt, eventWeight);
        if( jet_DoubleB_discr>0. ) h1_JetPt_BosonMatched_JetMass_Nsubj_DoubleB->Fill(jetPt, eventWeight);
      }

      // mass drop
      if( useMassDrop )
      {
        double fatJetMass = jetMass;
        double subjetMass = ( subjets.size()>1 ? std::max( subjets.at(0)->mass(), subjets.at(1)->mass() ) : 0. );
        if( useUncorrMassForMassDrop && groomedJetMatchFound && subjets.size()>1 )
        {
          fatJetMass = groomedJetMatch->correctedJet("Uncorrected").mass();
          subjetMass = std::max( subjets.at(0)->correctedJet("Uncorrected").mass(), subjets.at(1)->correctedJet("Uncorrected").mass() );
        }
        double massDrop = ( jetMass>0. ? subjetMass/fatJetMass : -10.);
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

    return;
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
