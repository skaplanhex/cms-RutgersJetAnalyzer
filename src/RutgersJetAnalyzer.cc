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
// $Id: RutgersJetAnalyzer.cc,v 1.23 2013/07/28 19:03:29 ferencek Exp $
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
    const bool             calculateMassDrop;
    const double           jetPtMin;
    const unsigned         jetPtBins;
    const double           jetPtBinWidth;
    const double           jetAbsEtaMax;
    const double           jetMassMin;
    const double           jetMassMax;
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

    TH2D *h2_JetPt_dRmatchedBhadrons_GSPbJets;

    TH1D *h1_JetPt;
    TH1D *h1_JetPt_BosonMatched;
    TH1D *h1_JetPt_BosonMatched_JetMass;
    TH1D *h1_JetPt_BosonMatched_JetMass_CSVL;
    TH1D *h1_JetPt_BosonMatched_JetMass_CSVM;
    TH1D *h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL;
    TH1D *h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVL;
    TH1D *h1_JetPt_BosonMatched_JetMass_SubJetMinCSVM;
    TH1D *h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVM;
    TH1D *h1_JetPt_BosonMatched_JetMass_DoubleB;
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
    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched_SubJetMinCSVL;
    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass_SubJetMinCSVL;

    TH1D *h1_JetCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_SubJetMinCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_SubJetMaxCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_JetDoubleBDiscr_BosonMatched_JetMass;

    TH2D *h2_JetPt_JetCSVL_BosonMatched_JetMass;
    TH2D *h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass;
    TH2D *h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass;

    std::map<std::string, TH2D*> h2_nPV_JetMass_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_MassDrop_Pt;

    std::map<std::string, TH2D*> h2_JetMass_nTracks_Pt;
    std::map<std::string, TH2D*> h2_JetMass_nSelectedTracks_Pt;
    std::map<std::string, TH2D*> h2_JetMass_tau2tau1_Pt;
    std::map<std::string, TH2D*> h2_JetMass_MassDrop_Pt;
    std::map<std::string, TH2D*> h2_JetMass_SubJetMinCSVL_Pt;
    std::map<std::string, TH2D*> h2_JetMass_SubJetMaxCSVL_Pt;
    std::map<std::string, TH2D*> h2_JetMass_TrackJetWidth_Pt;
    std::map<std::string, TH2D*> h2_JetMass_SelectedTrackJetWidth_Pt;
    std::map<std::string, TH2D*> h2_JetMass_maxdRTracks_Pt;
    std::map<std::string, TH2D*> h2_JetMass_maxdRSelectedTracks_Pt;

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
  calculateMassDrop(iConfig.getParameter<bool>("CalculateMassDrop")),
  jetPtMin(iConfig.getParameter<double>("JetPtMin")),
  jetPtBins(iConfig.getParameter<unsigned>("JetPtBins")),
  jetPtBinWidth(iConfig.getParameter<double>("JetPtBinWidth")),
  jetAbsEtaMax(iConfig.getParameter<double>("JetAbsEtaMax")),
  jetMassMin(iConfig.getParameter<double>("JetMassMin")),
  jetMassMax(iConfig.getParameter<double>("JetMassMax")),
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
    h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVL = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVL",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_SubJetMinCSVM = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_SubJetMinCSVM",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVM = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVM",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass_DoubleB = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass_DoubleB",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonDecayProdMatched = fs->make<TH1D>("h1_JetPt_BosonDecayProdMatched",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonDecayProdMatched_JetMass  = fs->make<TH1D>("h1_JetPt_BosonDecayProdMatched_JetMass",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetEta = fs->make<TH1D>("h1_JetEta",";#eta;",etaBins,etaMin,etaMax);
    h1_JetEta_BosonMatched = fs->make<TH1D>("h1_JetEta_BosonMatched",";#eta;",etaBins,etaMin,etaMax);
    h1_JetEta_BosonMatched_JetMass = fs->make<TH1D>("h1_JetEta_BosonMatched_JetMass",";#eta;",etaBins,etaMin,etaMax);

    h2_JetPt_dRmatchedBhadrons_GSPbJets = fs->make<TH2D>("h2_JetPt_dRmatchedBhadrons_GSPbJets",";p_{T} [GeV];#DeltaR",ptBins,ptMin,ptMax,dRBins,0.,1.6);
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
    h2_JetPt_SameMatchedBhadron_BosonMatched_SubJetMinCSVL = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched_SubJetMinCSVL",";p_{T} [GeV];Same matched B hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass",";p_{T} [GeV];Same matched B hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass_SubJetMinCSVL = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass_SubJetMinCSVL",";p_{T} [GeV];Same matched B hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);

    h1_JetCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetCSVDiscr_BosonMatched_JetMass",";Jet CSV Discr;",100,0.,1.);
    h1_SubJetMinCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_SubJetMinCSVDiscr_BosonMatched_JetMass",";SubJet min CSV Discr;",100,0.,1.);
    h1_SubJetMaxCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_SubJetMaxCSVDiscr_BosonMatched_JetMass",";SubJet max CSV Discr;",100,0.,1.);
    h1_JetDoubleBDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetDoubleBDiscr_BosonMatched_JetMass",";Jet DoubleB Discr;",100,0.,10.);

    h2_JetPt_JetCSVL_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_JetCSVL_BosonMatched_JetMass",";p_{T} [GeV];Jet CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);
    h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass",";p_{T} [GeV];SubJet min CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);
    h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass",";p_{T} [GeV];SubJet max CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);

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
        h2_nPV_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_nPV_MassDrop_Pt" + suffix).c_str(),(title + ";nPV;#mu=m_{subjet1}/m_{jet}").c_str(),pvBins,pvMin,pvMax,massDropBins,massDropMin,massDropMax);

        h2_JetMass_nTracks_Pt[suffix]         = fs->make<TH2D>(("h2_JetMass_nTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_nSelectedTracks_Pt[suffix] = fs->make<TH2D>(("h2_JetMass_nSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nSelectedTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        h2_JetMass_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_tau2tau1_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#tau_{2}/#tau_{1}").c_str(),massBins,massMin,massMax,tauBins,tauMin,tauMax);
        h2_JetMass_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_MassDrop_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#mu=m_{subjet1}/m_{jet}").c_str(),massBins,massMin,massMax,massDropBins,massDropMin,massDropMax);
        h2_JetMass_SubJetMinCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_SubJetMaxCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMaxCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet max CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_TrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_TrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_SelectedTrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SelectedTrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_maxdRTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);
        h2_JetMass_maxdRSelectedTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);

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
        h2_JetMass_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_MassDrop_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#mu=m_{subjet1}/m_{jet}").c_str(),massBins,massMin,massMax,massDropBins,massDropMin,massDropMax);
        h2_JetMass_SubJetMinCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_SubJetMaxCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMaxCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet max CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_TrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_TrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_SelectedTrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SelectedTrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_maxdRTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);
        h2_JetMass_maxdRSelectedTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);

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
        h2_JetMass_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_MassDrop_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#mu=m_{subjet1}/m_{jet}").c_str(),massBins,massMin,massMax,massDropBins,massDropMin,massDropMax);
        h2_JetMass_SubJetMinCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_SubJetMaxCSVL_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMaxCSVL_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet max CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_TrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_TrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_SelectedTrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SelectedTrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        h2_JetMass_maxdRTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);
        h2_JetMass_maxdRSelectedTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);

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

    // vector of pointers to status=3 b' or t' decay products
    std::vector<const reco::GenParticle*> resonanceDecayProducts;
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
        if( abs(it->pdgId()) == 5 && it->status()==2 ) { bFoundS2Quark = true; break; } // no need to continue looping after status=2 b quarks has been found
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
      int tPrimeCount = 0;
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
      {
        if( it->status() == 2 ) break; // to speed things up (only works with Pythia6)
        // count b' in the list of GenParticles
        if( abs(it->pdgId()) == 7 && it->status() == 3 ) ++bPrimeCount;
        // count t' in the list of GenParticles
        if( abs(it->pdgId()) == 8 && it->status() == 3 ) ++tPrimeCount;
        // only take status=3 quarks and charged leptons that appear after b' or t' have appeared in the list
        if( (bPrimeCount>0 || tPrimeCount>0 ) && it->status()==3 )
        {
          int dpPdgId = abs(it->pdgId());
          if( (dpPdgId>=1 && dpPdgId<=5) || dpPdgId==11 || dpPdgId==13 || dpPdgId==15  ) resonanceDecayProducts.push_back(&(*it));
        }
      }

      // loop over GenParticles and select bosons isolated from other b' or t' decay products
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
      {
        if( it->status() == 2 ) break; // to speed things up (only works with Pythia6)
        if( abs(it->pdgId()) == abs(bosonPdgId) && it->status() == 3 )
        {
          h1_BosonPt->Fill( it->pt(), eventWeight );
          h1_BosonEta->Fill( it->eta(), eventWeight );

          bool isIsolated = true;
          if( applyBosonIsolation )
          {
            for(std::vector<const reco::GenParticle*>::const_iterator dpIt = resonanceDecayProducts.begin(); dpIt != resonanceDecayProducts.end(); ++dpIt)
            {
              if( &(*it)==(*dpIt) ) continue; // skip the boson itself (should no longer happen since now comparing only to quarks and charged leptons)
              bool isBosonDecayProduct = false;
              if( abs(it->pdgId())==6 ) // special treatment for top quarks
              {
                for(unsigned i=0; i<it->numberOfDaughters(); ++i)
                {
                  if( it->daughter(i)->status()==2 ) continue; // only care about status=3 daughters
                  if( abs(it->daughter(i)->pdgId())==24 ) // if daughter is W
                  {
                    for(unsigned j=0; j<it->daughter(i)->numberOfDaughters(); ++j)
                    {
                      if( it->daughter(i)->daughter(j) == (*dpIt) )
                      {
                        isBosonDecayProduct = true;
                        break;
                      }
                    }
                  }
                  else
                  {
                    if( it->daughter(i) == (*dpIt) )
                    {
                      isBosonDecayProduct = true;
                      break;
                    }
                  }
                }
              }
              else
              {
                for(unsigned i=0; i<it->numberOfDaughters(); ++i)
                {
                  if( it->daughter(i) == (*dpIt) )
                  {
                    isBosonDecayProduct = true;
                    break;
                  }
                }
              }
              if( isBosonDecayProduct ) continue; // skip the boson decay products

              if( reco::deltaR( it->p4(), (*dpIt)->p4() ) < jetRadius ) isIsolated = false;
            }
          }

          if( !isIsolated ) continue;

          h1_BosonPt_Isolated->Fill( it->pt(), eventWeight );
          h1_BosonEta_Isolated->Fill( it->eta(), eventWeight );

          if( doBosonDecayProdSelection )
          {
            bool decayProductsFound = false;

            if( abs(it->pdgId())==6 ) // special treatment for top quarks
            {
              for(unsigned i=0; i<it->numberOfDaughters(); ++i)
              {
                if( it->daughter(i)->status()==2 ) continue; // only care about status=3 daughters
                if( abs(it->daughter(i)->pdgId())<=5 ) decayProducts[&(*it)].push_back(it->daughter(i)); // pick up a quark from the top decay

                if( abs(it->daughter(i)->pdgId())==24 ) // if top decay product is W
                {
                  for(unsigned j=0; j<it->daughter(i)->numberOfDaughters(); ++j)
                  {
                    if( it->daughter(i)->daughter(j)->status()==2 ) continue; // only care about status=3 daughters
                    for(std::vector<int>::const_iterator pdgIdIt = bosonDecayProdPdgIds.begin(); pdgIdIt != bosonDecayProdPdgIds.end(); ++pdgIdIt)
                    {
                      if( abs(it->daughter(i)->daughter(j)->pdgId()) == abs(*pdgIdIt) )
                      {
                        decayProductsFound = true;
                        decayProducts[&(*it)].push_back(it->daughter(i)->daughter(j));
                      }
                    }
                  }
                }
              }
            }
            else
            {
              for(unsigned i=0; i<it->numberOfDaughters(); ++i)
              {
                if( it->daughter(i)->status()==2 ) continue; // only care about status=3 daughters
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
            }

            if( decayProductsFound )
            {
              if( abs(it->pdgId())==6 ) // special treatment for top quarks
              {
                if( decayProducts[&(*it)].size()>3 ) edm::LogError("TooManyDecayProducts") << "More than three boson decay products found.";
                else if( decayProducts[&(*it)].size()<3 ) edm::LogError("TooFewDecayProducts") << "Less than three boson decay products found.";
                else
                {
                  bosons.push_back(&(*it));
                  h1_BosonPt_DecaySel->Fill( it->pt(), eventWeight );
                  h1_BosonEta_DecaySel->Fill( it->eta(), eventWeight );

                  double dRmax = -99.;
                  for(unsigned i=0; i<decayProducts[&(*it)].size(); ++i)
                  {
                    for(unsigned j=i+1; j<decayProducts[&(*it)].size(); ++j)
                    {
                      double dRtemp = reco::deltaR( decayProducts[&(*it)].at(i)->p4(), decayProducts[&(*it)].at(j)->p4() );
                      if( dRtemp>dRmax ) dRmax = dRtemp;
                    }
                  }
                  h2_BosonPt_dRdecay->Fill( it->pt(), dRmax, eventWeight );
                }
              }
              else
              {
                if( decayProducts[&(*it)].size()>2 ) edm::LogError("TooManyDecayProducts") << "More than two boson decay products found.";
                else if( decayProducts[&(*it)].size()<2 ) edm::LogError("TooFewDecayProducts") << "Less than two boson decay products found.";
                else
                {
                  bosons.push_back(&(*it));
                  h1_BosonPt_DecaySel->Fill( it->pt(), eventWeight );
                  h1_BosonEta_DecaySel->Fill( it->eta(), eventWeight );

                  h2_BosonPt_dRdecay->Fill( it->pt(), reco::deltaR( decayProducts[&(*it)].at(0)->p4(), decayProducts[&(*it)].at(1)->p4() ), eventWeight );
                }
              }
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

      if( abs(bosonPdgId)==6 && decayProducts[*bosonIt].size()>2 ) // special treatment for top quarks
      {
        for(PatJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it)
        {
          if( reco::deltaR( decayProducts[*bosonIt].at(0)->p4(), it->p4() ) < jetRadius &&
              reco::deltaR( decayProducts[*bosonIt].at(1)->p4(), it->p4() ) < jetRadius &&
              reco::deltaR( decayProducts[*bosonIt].at(2)->p4(), it->p4() ) < jetRadius )
          {
            h1_BosonPt_DecayProdMatched->Fill( (*bosonIt)->pt(), eventWeight );
            break;
          }
        }
      }
      else if( decayProducts[*bosonIt].size()>1 )
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
          std::vector<const reco::GenParticle*> matchedBHadrons;
          for(reco::GenParticleCollection::const_iterator gpIt = genParticles->begin(); gpIt != genParticles->end(); ++gpIt)
          {
            int id = abs(gpIt->pdgId());
            // skip GenParticle if not B hadron
            if ( !((id/100)%10 == 5 || (id/1000)%10 == 5) ) continue;

            // check if any of daughters is also B hadron
            bool hasBHadronDaughter = false;
            for(unsigned i=0; i<gpIt->numberOfDaughters(); ++i)
            {
              int dId = abs(gpIt->daughter(i)->pdgId());
              if ( (dId/100)%10 == 5 || (dId/1000)%10 == 5 ) { hasBHadronDaughter = true; break; }
            }
            if( hasBHadronDaughter ) continue; // skip excited B hadrons that have other B hadrons as daughters

            if( reco::deltaR( it->p4(), gpIt->p4() ) < jetRadius ) { ++nMatchedBHadrons; matchedBHadrons.push_back(&(*gpIt)); }
          }
          //std::cout << "nMatchedBHadrons: " << nMatchedBHadrons << std::endl;
          if( nMatchedBHadrons>=2 )
          {
            jetFlavor = 85; // custom jet flavor code for gluon splitting b jets
            h2_JetPt_dRmatchedBhadrons_GSPbJets->Fill( jetPt, reco::deltaR( matchedBHadrons.at(0)->p4(), matchedBHadrons.at(1)->p4() ), eventWeight );
          }
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
        // matching based on decay products only used for making some plots
        for(std::vector<const reco::GenParticle*>::const_iterator bosonIt = bosons.begin(); bosonIt != bosons.end(); ++bosonIt)
        {
          if( abs(bosonPdgId)==6 && decayProducts[*bosonIt].size()>2 ) // special treatment for top quarks
          {
            if( reco::deltaR( decayProducts[*bosonIt].at(0)->p4(), it->p4() ) < jetRadius &&
                reco::deltaR( decayProducts[*bosonIt].at(1)->p4(), it->p4() ) < jetRadius &&
                reco::deltaR( decayProducts[*bosonIt].at(2)->p4(), it->p4() ) < jetRadius )
            {
              h1_JetPt_BosonDecayProdMatched->Fill( jetPt, eventWeight );
              if( jetMass > jetMassMin && jetMass < jetMassMax )
                h1_JetPt_BosonDecayProdMatched_JetMass->Fill( jetPt, eventWeight );
              break;
            }
          }
          else if( decayProducts[*bosonIt].size()>1 )
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

      // get b-tag discriminators
      double jet_CSV_discr = it->bDiscriminator((bdiscriminator).c_str());
      double subJet1_CSV_discr = -999., subJet2_CSV_discr = -999.;
      if( subjets.size()>1 )
      {
        subJet1_CSV_discr = subjets.at(0)->bDiscriminator(bdiscriminator.c_str());
        subJet2_CSV_discr = subjets.at(1)->bDiscriminator(bdiscriminator.c_str());
      }
      double minSubJet_CSV_discr = std::min(subJet1_CSV_discr, subJet2_CSV_discr);
      double maxSubJet_CSV_discr = std::max(subJet1_CSV_discr, subJet2_CSV_discr);
      double jet_DoubleB_discr = it->bDiscriminator("doubleSecondaryVertexHighEffBJetTags");

      if( minSubJet_CSV_discr>0.244 )
      {
        if( bHadronMatchSubjet1 != genParticles->end() && bHadronMatchSubjet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadron_BosonMatched_SubJetMinCSVL->Fill(jetPt, ( bHadronMatchSubjet1 == bHadronMatchSubjet2 ? 1. : 0. ), eventWeight);
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

        if( subjets.at(0)->tagInfoSecondaryVertex("secondaryVertex")->nVertices() > 0
            && subjets.at(1)->tagInfoSecondaryVertex("secondaryVertex")->nVertices() > 0 )
        {
          std::string category = "DoubleSecondaryVertex";

          std::cout << category << ": ----------- START ------------" << std::endl;

          std::cout << category << ": Run, lumi, event: " << iEvent.id().run() << ", "
                                                          << iEvent.luminosityBlock() << ", "
                                                          << iEvent.id().event() << std::endl;
          std::cout << category << ": Fat jet pt, eta, phi, mass, jes: " << it->pt() << ", "
                                                                         << it->eta() << ", "
                                                                         << it->phi() << ", "
                                                                         << jetMass << ", "
                                                                         << it->pt()/it->correctedJet("Uncorrected").pt() << std::endl;
          std::cout << category << ": SubJet1 pt, eta, phi, mass, jes: " << subjets.at(0)->pt() << ", "
                                                                         << subjets.at(0)->eta() << ", "
                                                                         << subjets.at(0)->phi() << ", "
                                                                         << subjets.at(0)->mass() << ", "
                                                                         << subjets.at(0)->pt()/subjets.at(0)->correctedJet("Uncorrected").pt() << std::endl;
          std::cout << category << ": SubJet2 pt, eta, phi, mass, jes: " << subjets.at(1)->pt() << ", "
                                                                         << subjets.at(1)->eta() << ", "
                                                                         << subjets.at(1)->phi() << ", "
                                                                         << subjets.at(1)->mass() << ", "
                                                                         << subjets.at(1)->pt()/subjets.at(1)->correctedJet("Uncorrected").pt() << std::endl;
          std::cout << category << ": dR(fat jet, subjet1): "<< reco::deltaR( it->p4(), subjets.at(0)->p4() ) << std::endl;
          std::cout << category << ": dR(fat jet, subjet2): "<< reco::deltaR( it->p4(), subjets.at(1)->p4() ) << std::endl;
          std::cout << category << ": dR(subjet1, subjet2): "<< reco::deltaR( subjets.at(0)->p4(), subjets.at(1)->p4() ) << std::endl;

          std::cout << category << ": ------------ END -------------" << std::endl;
        }
      }

      h1_JetCSVDiscr_BosonMatched_JetMass->Fill( jet_CSV_discr, eventWeight);
      h1_SubJetMinCSVDiscr_BosonMatched_JetMass->Fill( minSubJet_CSV_discr, eventWeight);
      h1_SubJetMaxCSVDiscr_BosonMatched_JetMass->Fill( maxSubJet_CSV_discr, eventWeight);
      h1_JetDoubleBDiscr_BosonMatched_JetMass->Fill( jet_DoubleB_discr, eventWeight);

      h2_JetPt_JetCSVL_BosonMatched_JetMass->Fill(jetPt, jet_CSV_discr, eventWeight);;
      h2_JetPt_SubJetMinCSVL_BosonMatched_JetMass->Fill(jetPt, minSubJet_CSV_discr, eventWeight);
      h2_JetPt_SubJetMaxCSVL_BosonMatched_JetMass->Fill(jetPt, maxSubJet_CSV_discr, eventWeight);

      if( jet_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_CSVL->Fill(jetPt, eventWeight);
      if( jet_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_CSVM->Fill(jetPt, eventWeight);
      if( minSubJet_CSV_discr>0.244 )
      {
        h1_JetPt_BosonMatched_JetMass_SubJetMinCSVL->Fill(jetPt, eventWeight);
        if( bHadronMatchSubjet1 != genParticles->end() && bHadronMatchSubjet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass_SubJetMinCSVL->Fill(jetPt, ( bHadronMatchSubjet1 == bHadronMatchSubjet2 ? 1. : 0. ), eventWeight);
      }
      if( maxSubJet_CSV_discr>0.244 ) h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVL->Fill(jetPt, eventWeight);
      if( minSubJet_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_SubJetMinCSVM->Fill(jetPt, eventWeight);
      if( maxSubJet_CSV_discr>0.679 ) h1_JetPt_BosonMatched_JetMass_SubJetMaxCSVM->Fill(jetPt, eventWeight);
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

      int nTracks = it->associatedTracks().size();
      int nSelectedTracks = ( it->hasTagInfo("impactParameter") ? it->tagInfoTrackIP("impactParameter")->selectedTracks().size() : -99 );

      double trackJetWidth = 0.;
      double trackPtSum = 0.;
      double selectedTrackJetWidth = 0.;
      double selectedTrackPtSum = 0.;
      double maxdRTracks = -99.;
      double maxdRSelectedTracks = -99.;

      for(int i=0; i<nTracks; ++i)
      {
        double dRJetTrack = reco::deltaR( it->eta(), it->phi(), it->associatedTracks().at(i)->eta(), it->associatedTracks().at(i)->phi() );
        double trackPt = it->associatedTracks().at(i)->pt();
        trackJetWidth+=(dRJetTrack*trackPt);
        trackPtSum+=trackPt;

        for(int j=0; j<nTracks; ++j)
        {
          double dRTrkTrk = reco::deltaR( it->associatedTracks().at(i)->eta(), it->associatedTracks().at(i)->phi(), it->associatedTracks().at(j)->eta(), it->associatedTracks().at(j)->phi() );
          if( dRTrkTrk > maxdRTracks ) maxdRTracks = dRTrkTrk;
        }
      }
      trackJetWidth/=trackPtSum;

      for(int i=0; i<nSelectedTracks; ++i)
      {
        double dRJetTrack = reco::deltaR( it->eta(), it->phi(), it->tagInfoTrackIP("impactParameter")->selectedTracks().at(i)->eta(), it->tagInfoTrackIP("impactParameter")->selectedTracks().at(i)->phi() );
        double trackPt = it->tagInfoTrackIP("impactParameter")->selectedTracks().at(i)->pt();
        selectedTrackJetWidth+=(dRJetTrack*trackPt);
        selectedTrackPtSum+=trackPt;

        for(int j=0; j<nSelectedTracks; ++j)
        {
          double dRTrkTrk = reco::deltaR( it->tagInfoTrackIP("impactParameter")->selectedTracks().at(i)->eta(), it->tagInfoTrackIP("impactParameter")->selectedTracks().at(i)->phi(), it->tagInfoTrackIP("impactParameter")->selectedTracks().at(j)->eta(), it->tagInfoTrackIP("impactParameter")->selectedTracks().at(j)->phi() );
          if( dRTrkTrk > maxdRSelectedTracks ) maxdRSelectedTracks = dRTrkTrk;
        }
      }
      selectedTrackJetWidth/=selectedTrackPtSum;

      // fill various 2D histograms
      suffix = Form("%.0ftoInf",jetPtMin);
      h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
      h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
      h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);

      h2_JetMass_nTracks_Pt[suffix]->Fill(jetMass, nTracks, eventWeight);
      h2_JetMass_nSelectedTracks_Pt[suffix]->Fill(jetMass, nSelectedTracks, eventWeight);
      h2_JetMass_tau2tau1_Pt[suffix]->Fill(jetMass, tau2overtau1, eventWeight);
      h2_JetMass_SubJetMinCSVL_Pt[suffix]->Fill(jetMass, minSubJet_CSV_discr, eventWeight);
      h2_JetMass_SubJetMaxCSVL_Pt[suffix]->Fill(jetMass, maxSubJet_CSV_discr, eventWeight);
      h2_JetMass_TrackJetWidth_Pt[suffix]->Fill(jetMass, trackJetWidth, eventWeight);
      h2_JetMass_SelectedTrackJetWidth_Pt[suffix]->Fill(jetMass, selectedTrackJetWidth, eventWeight);
      h2_JetMass_maxdRTracks_Pt[suffix]->Fill(jetMass, maxdRTracks, eventWeight);
      h2_JetMass_maxdRSelectedTracks_Pt[suffix]->Fill(jetMass, maxdRSelectedTracks, eventWeight);

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
          h2_JetMass_SubJetMaxCSVL_Pt[suffix]->Fill(jetMass, maxSubJet_CSV_discr, eventWeight);
          h2_JetMass_TrackJetWidth_Pt[suffix]->Fill(jetMass, trackJetWidth, eventWeight);
          h2_JetMass_SelectedTrackJetWidth_Pt[suffix]->Fill(jetMass, selectedTrackJetWidth, eventWeight);
          h2_JetMass_maxdRTracks_Pt[suffix]->Fill(jetMass, maxdRTracks, eventWeight);
          h2_JetMass_maxdRSelectedTracks_Pt[suffix]->Fill(jetMass, maxdRSelectedTracks, eventWeight);

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
        h2_JetMass_SubJetMaxCSVL_Pt[suffix]->Fill(jetMass, maxSubJet_CSV_discr, eventWeight);
        h2_JetMass_TrackJetWidth_Pt[suffix]->Fill(jetMass, trackJetWidth, eventWeight);
        h2_JetMass_SelectedTrackJetWidth_Pt[suffix]->Fill(jetMass, selectedTrackJetWidth, eventWeight);
        h2_JetMass_maxdRTracks_Pt[suffix]->Fill(jetMass, maxdRTracks, eventWeight);
        h2_JetMass_maxdRSelectedTracks_Pt[suffix]->Fill(jetMass, maxdRSelectedTracks, eventWeight);

        h2_nTracks_tau2tau1_Pt[suffix]->Fill(nTracks, tau2overtau1, eventWeight);
        h2_nSelectedTracks_tau2tau1_Pt[suffix]->Fill(nSelectedTracks, tau2overtau1, eventWeight);

        h2_nTracks_SubJetMinCSVL_Pt[suffix]->Fill(nTracks, minSubJet_CSV_discr, eventWeight);
        h2_nSelectedTracks_SubJetMinCSVL_Pt[suffix]->Fill(nSelectedTracks, minSubJet_CSV_discr, eventWeight);
        h2_tau2tau1_SubJetMinCSVL_Pt[suffix]->Fill(tau2overtau1, minSubJet_CSV_discr, eventWeight);
      }

      // mass drop
      if( calculateMassDrop )
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
        h2_JetMass_MassDrop_Pt[suffix]->Fill(jetMass, massDrop, eventWeight);
        for(unsigned i=0; i<jetPtBins; ++i)
        {
          if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
          {
            suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
            h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
            h2_JetMass_MassDrop_Pt[suffix]->Fill(jetMass, massDrop, eventWeight);
          }
        }
        if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
        {
          suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
          h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
          h2_JetMass_MassDrop_Pt[suffix]->Fill(jetMass, massDrop, eventWeight);
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
