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
// $Id: RutgersJetAnalyzer.cc,v 1.25 2013/08/09 22:47:09 ferencek Exp $
//
//


// system include files
#include <memory>

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
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

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
    const bool             useSubJets;
    const edm::InputTag    groomedBasicJetsTag;
    const bool             useAK5Jets;
    const edm::InputTag    ak5JetsTag;
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
    const bool             useUncorrMassForMassDrop;
    const std::string      bdiscriminator;
    const bool             doJetFlavor;
    const std::vector<int> jetFlavorPdgIds;
    const bool             useAltGSPbDef;
    const bool             findGluonSplitting;
    const bool             findMatrixElement;
    const bool             eventDisplayPrintout;

    edm::Service<TFileService> fs;

    edm::ESHandle<TransientTrackBuilder> trackBuilder;

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

    TH1D *h1_SubJetPt_BosonMatched_JetMass;
    TProfile *p1_SubJetPt_TotalTracks_BosonMatched_JetMass;
    TProfile *p1_SubJetPt_SharedTracks_BosonMatched_JetMass;
    TProfile *p1_SubJetPt_SharedTracksRatio_BosonMatched_JetMass;
    TProfile *p1_SubJetPt_TotalVertexTracks_BosonMatched_JetMass;
    TProfile *p1_SubJetPt_SharedVertexTracks_BosonMatched_JetMass;
    TProfile *p1_SubJetPt_SharedVertexTracksRatio_BosonMatched_JetMass;
    TProfile *p1_JetPt_TotalTracks_BosonMatched_JetMass;
    TProfile *p1_JetPt_SharedTracks_BosonMatched_JetMass;
    TProfile *p1_JetPt_SharedTracksRatio_BosonMatched_JetMass;
    TProfile *p1_JetPt_TotalVertexTracks_BosonMatched_JetMass;
    TProfile *p1_JetPt_SharedVertexTracks_BosonMatched_JetMass;
    TProfile *p1_JetPt_SharedVertexTracksRatio_BosonMatched_JetMass;
    TProfile *p1_dRsubjets_TotalTracks_BosonMatched_JetMass;
    TProfile *p1_dRsubjets_SharedTracks_BosonMatched_JetMass;
    TProfile *p1_dRsubjets_SharedTracksRatio_BosonMatched_JetMass;

    TH2D *h2_JetPt_JetPtOverBosonPt;
    TH2D *h2_JetPt_JetPtOverGenJetPt;
    TH2D *h2_JetPt_JetPtOverGenJetPt_BosonMatched;
    TH2D *h2_JetPt_JetMass;
    TH2D *h2_JetPt_JetMass_BosonMatched;
    TH2D *h2_JetPt_dRsubjets_BosonMatched;
    TH2D *h2_JetPt_dRsubjets_BosonMatched_JetMass;
    TH2D *h2_JetPt_dRak5jets_BosonMatched;
    TH2D *h2_JetPt_dRak5jets_BosonMatched_JetMass;

    TH2D *h2_JetPt_mindRjetBhadron_BosonMatched;
    TH2D *h2_JetPt_mindRSubjet1Bhadron_BosonMatched;
    TH2D *h2_JetPt_mindRSubjet2Bhadron_BosonMatched;
    TH2D *h2_JetPt_mindRak5jetBhadron_BosonMatched;
    TH2D *h2_JetPt_mindRak5jet1Bhadron_BosonMatched;
    TH2D *h2_JetPt_mindRak5jet2Bhadron_BosonMatched;
    TH2D *h2_JetPt_mindRjetBhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_mindRak5jetBhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_mindRak5jet1Bhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_mindRak5jet2Bhadron_BosonMatched_JetMass;

    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched;
    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched_SubJetMinCSVL;
    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass;
    TH2D *h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass_SubJetMinCSVL;
    TH2D *h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched;
    TH2D *h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched_JetMass;

    TH1D *h1_JetCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_SubJetMinCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_SubJetMaxCSVDiscr_BosonMatched_JetMass;
    TH1D *h1_JetDoubleBDiscr_BosonMatched_JetMass;

    TH2D *h2_JetPt_JetCSV_BosonMatched_JetMass;
    TH2D *h2_JetPt_SubJetMinCSV_BosonMatched_JetMass;
    TH2D *h2_JetPt_SubJetMaxCSV_BosonMatched_JetMass;

    TH2D *h2_JetPt_JetJP_BosonMatched_JetMass;
    TH2D *h2_JetPt_SubJetMinJP_BosonMatched_JetMass;
    TH2D *h2_JetPt_SubJetMaxJP_BosonMatched_JetMass;

    TH2D *h2_JetPt_AK5JetCSV_BosonMatched_JetMass;
    TH2D *h2_JetPt_AK5JetMinCSV_BosonMatched_JetMass;

    TH2D *h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0to0p2;
    TH2D *h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p2to0p4;
    TH2D *h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p4to0p6;
    TH2D *h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p6to0p8;

    std::map<std::string, TH2D*> h2_nPV_JetMass_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2_Pt;
    std::map<std::string, TH2D*> h2_nPV_tau2tau1_Pt;
    std::map<std::string, TH2D*> h2_nPV_MassDrop_Pt;

    //std::map<std::string, TH2D*> h2_JetMass_nTracks_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_nSelectedTracks_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_tau2tau1_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_MassDrop_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_SubJetMinCSV_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_SubJetMaxCSV_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_TrackJetWidth_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_SelectedTrackJetWidth_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_maxdRTracks_Pt;
    //std::map<std::string, TH2D*> h2_JetMass_maxdRSelectedTracks_Pt;

    //std::map<std::string, TH2D*> h2_nTracks_tau2tau1_Pt;
    //std::map<std::string, TH2D*> h2_nSelectedTracks_tau2tau1_Pt;

    //std::map<std::string, TH2D*> h2_nTracks_SubJetMinCSV_Pt;
    //std::map<std::string, TH2D*> h2_nSelectedTracks_SubJetMinCSV_Pt;
    //std::map<std::string, TH2D*> h2_tau2tau1_SubJetMinCSV_Pt;
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
  useSubJets(iConfig.getParameter<bool>("UseSubJets")),
  groomedBasicJetsTag(iConfig.getParameter<edm::InputTag>("GroomedBasicJetsTag")),
  useAK5Jets( iConfig.exists("UseAK5Jets") ? iConfig.getParameter<bool>("UseAK5Jets") : false ),
  ak5JetsTag( iConfig.exists("AK5JetsTag") ? iConfig.getParameter<edm::InputTag>("AK5JetsTag") : edm::InputTag("selectedPatJetsAK5PF") ),
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
  useUncorrMassForMassDrop( iConfig.exists("UseUncorrMassForMassDrop") ? iConfig.getParameter<bool>("UseUncorrMassForMassDrop") : true ),
  bdiscriminator(iConfig.getParameter<std::string>("Bdiscriminator")),
  doJetFlavor(iConfig.getParameter<bool>("DoJetFlavor")),
  jetFlavorPdgIds(iConfig.getParameter<std::vector<int> >("JetFlavorPdgIds")),
  useAltGSPbDef( iConfig.exists("UseAltGSPbDef") ? iConfig.getParameter<bool>("UseAltGSPbDef") : false ),
  findGluonSplitting( iConfig.exists("FindGluonSplitting") ? iConfig.getParameter<bool>("FindGluonSplitting") : false ),
  findMatrixElement( iConfig.exists("FindMatrixElement") ? iConfig.getParameter<bool>("FindMatrixElement") : false ),
  eventDisplayPrintout( iConfig.exists("EventDisplayPrintout") ? iConfig.getParameter<bool>("EventDisplayPrintout") : false )

{
    //now do what ever initialization is needed
    int pvBins=50;
    double pvMin=-0.5, pvMax=49.5;
    //int trackBins=200;
    //double trackMin=-0.5, trackMax=199.5;
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

    h1_SubJetPt_BosonMatched_JetMass = fs->make<TH1D>("h1_SubJetPt_BosonMatched_JetMass",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    p1_SubJetPt_TotalTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_SubJetPt_TotalTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_SubJetPt_SharedTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_SubJetPt_SharedTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_SubJetPt_SharedTracksRatio_BosonMatched_JetMass = fs->make<TProfile>("p1_SubJetPt_SharedTracksRatio_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_SubJetPt_TotalVertexTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_SubJetPt_TotalVertexTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_SubJetPt_SharedVertexTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_SubJetPt_SharedVertexTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_SubJetPt_SharedVertexTracksRatio_BosonMatched_JetMass = fs->make<TProfile>("p1_SubJetPt_SharedVertexTracksRatio_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_JetPt_TotalTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_JetPt_TotalTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_JetPt_SharedTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_JetPt_SharedTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_JetPt_SharedTracksRatio_BosonMatched_JetMass = fs->make<TProfile>("p1_JetPt_SharedTracksRatio_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_JetPt_TotalVertexTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_JetPt_TotalVertexTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_JetPt_SharedVertexTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_JetPt_SharedVertexTracks_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_JetPt_SharedVertexTracksRatio_BosonMatched_JetMass = fs->make<TProfile>("p1_JetPt_SharedVertexTracksRatio_BosonMatched_JetMass",";p_{T} [GeV];",20,ptMin,ptMax);
    p1_dRsubjets_TotalTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_dRsubjets_TotalTracks_BosonMatched_JetMass",";#DeltaR(subjet_{1},subjet_{2});",8,0.,0.8);
    p1_dRsubjets_SharedTracks_BosonMatched_JetMass = fs->make<TProfile>("p1_dRsubjets_SharedTracks_BosonMatched_JetMass",";#DeltaR(subjet_{1},subjet_{2});",8,0.,0.8);
    p1_dRsubjets_SharedTracksRatio_BosonMatched_JetMass = fs->make<TProfile>("p1_dRsubjets_SharedTracksRatio_BosonMatched_JetMass",";#DeltaR(subjet_{1},subjet_{2});",8,0.,0.8);

    h2_JetPt_dRmatchedBhadrons_GSPbJets = fs->make<TH2D>("h2_JetPt_dRmatchedBhadrons_GSPbJets",";p_{T} [GeV];#DeltaR",ptBins,ptMin,ptMax,dRBins,0.,1.6);
    h2_JetPt_JetPtOverBosonPt = fs->make<TH2D>("h2_JetPt_JetPtOverBosonPt",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{boson}",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetPtOverGenJetPt = fs->make<TH2D>("h2_JetPt_JetPtOverGenJetPt",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{genjet}",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetPtOverGenJetPt_BosonMatched = fs->make<TH2D>("h2_JetPt_JetPtOverGenJetPt_BosonMatched",";p_{T} [GeV];p_{T}^{jet}/p_{T}^{genjet}",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_JetMass = fs->make<TH2D>("h2_JetPt_JetMass",";p_{T} [GeV];m_{jet} [GeV]",ptBins,ptMin,ptMax,massBins,massMin,massMax);
    h2_JetPt_JetMass_BosonMatched = fs->make<TH2D>("h2_JetPt_JetMass_BosonMatched",";p_{T} [GeV];m_{jet} [GeV]",ptBins,ptMin,ptMax,massBins,massMin,massMax);
    h2_JetPt_dRsubjets_BosonMatched = fs->make<TH2D>("h2_JetPt_dRsubjets_BosonMatched",";p_{T} [GeV];#DeltaR(subjet_{1},subjet_{2})",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);
    h2_JetPt_dRsubjets_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_dRsubjets_BosonMatched_JetMass",";p_{T} [GeV];#DeltaR(subjet_{1},subjet_{2})",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);
    h2_JetPt_dRak5jets_BosonMatched = fs->make<TH2D>("h2_JetPt_dRak5jets_BosonMatched",";p_{T} [GeV];#DeltaR(AK5 jet_{1},AK5 jet_{2})",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);
    h2_JetPt_dRak5jets_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_dRak5jets_BosonMatched_JetMass",";p_{T} [GeV];#DeltaR(AK5 jet_{1},AK5 jet_{2})",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);

    h2_JetPt_mindRjetBhadron_BosonMatched     = fs->make<TH2D>("h2_JetPt_mindRjetBhadron_BosonMatched",    ";p_{T} [GeV];min #DeltaR(fat jet,b hadron)",    ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRSubjet1Bhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_mindRSubjet1Bhadron_BosonMatched",";p_{T} [GeV];min #DeltaR(subjet_{1},b hadron)", ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRSubjet2Bhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_mindRSubjet2Bhadron_BosonMatched",";p_{T} [GeV];min #DeltaR(subjet_{2},b hadron)", ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRak5jetBhadron_BosonMatched  = fs->make<TH2D>("h2_JetPt_mindRak5jetBhadron_BosonMatched", ";p_{T} [GeV];min #DeltaR(AK5 jet,b hadron)",    ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRak5jet1Bhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_mindRak5jet1Bhadron_BosonMatched",";p_{T} [GeV];min #DeltaR(AK5 jet_{1},b hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRak5jet2Bhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_mindRak5jet2Bhadron_BosonMatched",";p_{T} [GeV];min #DeltaR(AK5 jet_{2},b hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRjetBhadron_BosonMatched_JetMass     = fs->make<TH2D>("h2_JetPt_mindRjetBhadron_BosonMatched_JetMass",    ";p_{T} [GeV];min #DeltaR(fat jet,b hadron)",    ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass",";p_{T} [GeV];min #DeltaR(subjet_{1},b hadron)", ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass",";p_{T} [GeV];min #DeltaR(subjet_{2},b hadron)", ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRak5jetBhadron_BosonMatched_JetMass  = fs->make<TH2D>("h2_JetPt_mindRak5jetBhadron_BosonMatched_JetMass", ";p_{T} [GeV];min #DeltaR(AK5 jet,b hadron)",    ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRak5jet1Bhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_mindRak5jet1Bhadron_BosonMatched_JetMass",";p_{T} [GeV];min #DeltaR(AK5 jet_{1},b hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);
    h2_JetPt_mindRak5jet2Bhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_mindRak5jet2Bhadron_BosonMatched_JetMass",";p_{T} [GeV];min #DeltaR(AK5 jet_{2},b hadron)",ptBins,ptMin,ptMax,dRBins,0.,2.);

    h2_JetPt_SameMatchedBhadron_BosonMatched = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched",";p_{T} [GeV];Same matched b hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadron_BosonMatched_SubJetMinCSVL = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched_SubJetMinCSVL",";p_{T} [GeV];Same matched b hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass",";p_{T} [GeV];Same matched b hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass_SubJetMinCSVL = fs->make<TH2D>("h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass_SubJetMinCSVL",";p_{T} [GeV];Same matched b hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched = fs->make<TH2D>("h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched",";p_{T} [GeV];Same matched b hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);
    h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched_JetMass",";p_{T} [GeV];Same matched b hadron",ptBins,ptMin,ptMax,2,-0.5,1.5);

    h1_JetCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetCSVDiscr_BosonMatched_JetMass",";Jet CSV Discr;",100,0.,1.);
    h1_SubJetMinCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_SubJetMinCSVDiscr_BosonMatched_JetMass",";SubJet min CSV Discr;",100,0.,1.);
    h1_SubJetMaxCSVDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_SubJetMaxCSVDiscr_BosonMatched_JetMass",";SubJet max CSV Discr;",100,0.,1.);
    h1_JetDoubleBDiscr_BosonMatched_JetMass = fs->make<TH1D>("h1_JetDoubleBDiscr_BosonMatched_JetMass",";Jet DoubleB Discr;",100,0.,10.);

    h2_JetPt_JetCSV_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_JetCSV_BosonMatched_JetMass",";p_{T} [GeV];Jet CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);
    h2_JetPt_SubJetMinCSV_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SubJetMinCSV_BosonMatched_JetMass",";p_{T} [GeV];SubJet min CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);
    h2_JetPt_SubJetMaxCSV_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SubJetMaxCSV_BosonMatched_JetMass",";p_{T} [GeV];SubJet max CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);

    h2_JetPt_JetJP_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_JetJP_BosonMatched_JetMass",";p_{T} [GeV];Jet JP Discr",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_SubJetMinJP_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SubJetMinJP_BosonMatched_JetMass",";p_{T} [GeV];SubJet min JP Discr",ptBins,ptMin,ptMax,100,0.,2.);
    h2_JetPt_SubJetMaxJP_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_SubJetMaxJP_BosonMatched_JetMass",";p_{T} [GeV];SubJet max JP Discr",ptBins,ptMin,ptMax,100,0.,2.);

    h2_JetPt_AK5JetCSV_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_AK5JetCSV_BosonMatched_JetMass",";p_{T} [GeV];AK5 Jet CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);
    h2_JetPt_AK5JetMinCSV_BosonMatched_JetMass = fs->make<TH2D>("h2_JetPt_AK5JetMinCSV_BosonMatched_JetMass",";p_{T} [GeV];AK5 Jets min CSV Discr",ptBins,ptMin,ptMax,100,0.,1.);

    h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0to0p2   = fs->make<TH2D>("h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0to0p2",  ";Subjet_{1} Discr;Subjet_{2} Discr",100,0.,1.,100,0.,1.);
    h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p2to0p4 = fs->make<TH2D>("h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p2to0p4",";Subjet_{1} Discr;Subjet_{2} Discr",100,0.,1.,100,0.,1.);
    h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p4to0p6 = fs->make<TH2D>("h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p4to0p6",";Subjet_{1} Discr;Subjet_{2} Discr",100,0.,1.,100,0.,1.);
    h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p6to0p8 = fs->make<TH2D>("h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p6to0p8",";Subjet_{1} Discr;Subjet_{2} Discr",100,0.,1.,100,0.,1.);

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

        //h2_JetMass_nTracks_Pt[suffix]         = fs->make<TH2D>(("h2_JetMass_nTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        //h2_JetMass_nSelectedTracks_Pt[suffix] = fs->make<TH2D>(("h2_JetMass_nSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nSelectedTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        //h2_JetMass_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_tau2tau1_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#tau_{2}/#tau_{1}").c_str(),massBins,massMin,massMax,tauBins,tauMin,tauMax);
        //h2_JetMass_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_MassDrop_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#mu=m_{subjet1}/m_{jet}").c_str(),massBins,massMin,massMax,massDropBins,massDropMin,massDropMax);
        //h2_JetMass_SubJetMinCSV_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_SubJetMaxCSV_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMaxCSV_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet max CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_TrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_TrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_SelectedTrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SelectedTrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_maxdRTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);
        //h2_JetMass_maxdRSelectedTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);

        //h2_nTracks_tau2tau1_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);
        //h2_nSelectedTracks_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nSelectedTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);

        //h2_nTracks_SubJetMinCSV_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";nTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        //h2_nSelectedTracks_SubJetMinCSV_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";nSelectedTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        //h2_tau2tau1_SubJetMinCSV_Pt[suffix]        = fs->make<TH2D>(("h2_tau2tau1_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";#tau_{2}/#tau_{1};SubJet min CSV Discr").c_str(),tauBins,tauMin,tauMax,100,0.,1.);
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

        //h2_JetMass_nTracks_Pt[suffix]         = fs->make<TH2D>(("h2_JetMass_nTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        //h2_JetMass_nSelectedTracks_Pt[suffix] = fs->make<TH2D>(("h2_JetMass_nSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nSelectedTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        //h2_JetMass_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_tau2tau1_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#tau_{2}/#tau_{1}").c_str(),massBins,massMin,massMax,tauBins,tauMin,tauMax);
        //h2_JetMass_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_MassDrop_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#mu=m_{subjet1}/m_{jet}").c_str(),massBins,massMin,massMax,massDropBins,massDropMin,massDropMax);
        //h2_JetMass_SubJetMinCSV_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_SubJetMaxCSV_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMaxCSV_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet max CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_TrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_TrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_SelectedTrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SelectedTrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_maxdRTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);
        //h2_JetMass_maxdRSelectedTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);

        //h2_nTracks_tau2tau1_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);
        //h2_nSelectedTracks_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nSelectedTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);

        //h2_nTracks_SubJetMinCSV_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";nTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        //h2_nSelectedTracks_SubJetMinCSV_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";nSelectedTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        //h2_tau2tau1_SubJetMinCSV_Pt[suffix]        = fs->make<TH2D>(("h2_tau2tau1_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";#tau_{2}/#tau_{1};SubJet min CSV Discr").c_str(),tauBins,tauMin,tauMax,100,0.,1.);
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

        //h2_JetMass_nTracks_Pt[suffix]         = fs->make<TH2D>(("h2_JetMass_nTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        //h2_JetMass_nSelectedTracks_Pt[suffix] = fs->make<TH2D>(("h2_JetMass_nSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];nSelectedTracks").c_str(),massBins,massMin,massMax,trackBins,trackMin,trackMax);
        //h2_JetMass_tau2tau1_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_tau2tau1_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#tau_{2}/#tau_{1}").c_str(),massBins,massMin,massMax,tauBins,tauMin,tauMax);
        //h2_JetMass_MassDrop_Pt[suffix]        = fs->make<TH2D>(("h2_JetMass_MassDrop_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];#mu=m_{subjet1}/m_{jet}").c_str(),massBins,massMin,massMax,massDropBins,massDropMin,massDropMax);
        //h2_JetMass_SubJetMinCSV_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet min CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_SubJetMaxCSV_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SubJetMaxCSV_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];SubJet max CSV Discr").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_TrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_TrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_SelectedTrackJetWidth_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_SelectedTrackJetWidth_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track track-jet width").c_str(),massBins,massMin,massMax,100,0.,1.);
        //h2_JetMass_maxdRTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);
        //h2_JetMass_maxdRSelectedTracks_Pt[suffix]   = fs->make<TH2D>(("h2_JetMass_maxdRSelectedTracks_Pt" + suffix).c_str(),(title + ";m_{jet} [GeV];Selected track max#DeltaR(trk,trk)").c_str(),massBins,massMin,massMax,100,0.,2.);

        //h2_nTracks_tau2tau1_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);
        //h2_nSelectedTracks_tau2tau1_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_tau2tau1_Pt" + suffix).c_str(),(title + ";nSelectedTracks;#tau_{2}/#tau_{1}").c_str(),trackBins,trackMin,trackMax,tauBins,tauMin,tauMax);

        //h2_nTracks_SubJetMinCSV_Pt[suffix]         = fs->make<TH2D>(("h2_nTracks_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";nTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        //h2_nSelectedTracks_SubJetMinCSV_Pt[suffix] = fs->make<TH2D>(("h2_nSelectedTracks_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";nSelectedTracks;SubJet min CSV Discr").c_str(),trackBins,trackMin,trackMax,100,0.,1.);
        //h2_tau2tau1_SubJetMinCSV_Pt[suffix]        = fs->make<TH2D>(("h2_tau2tau1_SubJetMinCSV_Pt" + suffix).c_str(),(title + ";#tau_{2}/#tau_{1};SubJet min CSV Discr").c_str(),tauBins,tauMin,tauMax,100,0.,1.);
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

    edm::Handle<PatJetCollection> groomedBasicJets;
    if( useSubJets ) iEvent.getByLabel(groomedBasicJetsTag,groomedBasicJets);

    edm::Handle<PatJetCollection> ak5Jets;
    if( useAK5Jets ) iEvent.getByLabel(ak5JetsTag,ak5Jets);

    edm::Handle<reco::VertexCollection> PVs;
    iEvent.getByLabel(pvTag,PVs);

    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);


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

        if( useAltGSPbDef ) // use alternative gluon splitting b-jet definition based on the number of b hadrons inside the jet cone
        {
          int nMatchedBHadrons = 0;
          std::vector<const reco::GenParticle*> matchedBHadrons;
          for(reco::GenParticleCollection::const_iterator gpIt = genParticles->begin(); gpIt != genParticles->end(); ++gpIt)
          {
            int id = abs(gpIt->pdgId());
            // skip GenParticle if not b hadron
            if ( !((id/100)%10 == 5 || (id/1000)%10 == 5) ) continue;

            // check if any of daughters is also b hadron
            bool hasBHadronDaughter = false;
            for(unsigned i=0; i<gpIt->numberOfDaughters(); ++i)
            {
              int dId = abs(gpIt->daughter(i)->pdgId());
              if ( (dId/100)%10 == 5 || (dId/1000)%10 == 5 ) { hasBHadronDaughter = true; break; }
            }
            if( hasBHadronDaughter ) continue; // skip excited b hadrons that have other b hadrons as daughters

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


      double jetMass = it->userFloat("ca8PFJetsCHSPrunedMass");

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


      bool groomedBasicJetMatchFound = false;
      PatJetCollection::const_iterator groomedBasicJetMatch;
      // try to find matching groomed fat jet
      double dR = jetRadius;
      for(PatJetCollection::const_iterator gbjIt = groomedBasicJets->begin(); gbjIt != groomedBasicJets->end(); ++gbjIt)
      {
        double dR_temp = reco::deltaR( it->p4(), gbjIt->p4() );
        if( dR_temp < dR )
        {
          dR = dR_temp;
          groomedBasicJetMatch = gbjIt;
          groomedBasicJetMatchFound = true;
        }
      }
      // vector of pointers to subjets
      std::vector<const pat::Jet*> subjets;
      if( useSubJets )
      {
        //std::cout << "number of subjets: " << groomedBasicJetMatch->numberOfDaughters() << std::endl;
        //std::cout << "jet pt: " << it->correctedJet("Uncorrected").p4().pt() << " eta=" << it->correctedJet("Uncorrected").p4().eta() << " phi=" << it->correctedJet("Uncorrected").p4().phi() << " nd=" << it->numberOfDaughters() << std::endl;
        //std::cout << "groomed basic jet pt: " << groomedBasicJetMatch->p4().pt() << " eta=" << groomedBasicJetMatch->p4().eta() << " phi=" << groomedBasicJetMatch->p4().phi() << std::endl;
        if ( groomedBasicJetMatchFound )
        {
          for(unsigned d=0; d<groomedBasicJetMatch->numberOfDaughters(); ++d)
          {
            //std::cout << "subjet " << d << ": pt=" << groomedBasicJetMatch->daughter(d)->p4().pt()  << " eta=" << groomedBasicJetMatch->daughter(d)->p4().eta()  << " phi=" << groomedBasicJetMatch->daughter(d)->p4().phi() << std::endl;
            const reco::Candidate *subjet =  groomedBasicJetMatch->daughter(d);
            const pat::Jet *patsubjet = dynamic_cast<const pat::Jet*>(subjet);
            subjets.push_back(patsubjet);
          }
        }

        if( subjets.size()<2 )
          edm::LogError("TooFewSubjets") << "Less than two subjets (" << subjets.size() << ") found.";
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


      // vector of pointers to matching AK5 jets (matching radius is 0.6)
      std::vector<const pat::Jet*> matchedAK5Jets;
      if( useAK5Jets )
      {
        for(PatJetCollection::const_iterator jIt = ak5Jets->begin(); jIt != ak5Jets->end(); ++jIt)
        {
          double dR_temp = reco::deltaR( it->p4(), jIt->p4() );
          if( dR_temp < 0.6 ) matchedAK5Jets.push_back(&(*jIt));
        }
      }
      // sort matched AK5 jets by increasing b-tag discriminator
      std::vector<unsigned> sortedMatchedAK5JetsIdx;
      if( matchedAK5Jets.size()>1 )
      {
        std::multimap<double, unsigned> sortedMatchedAK5Jets;
        for(unsigned i = 0; i<matchedAK5Jets.size(); ++i)
          sortedMatchedAK5Jets.insert(std::make_pair(matchedAK5Jets.at(i)->bDiscriminator(bdiscriminator.c_str()), i));

        for(std::multimap<double, unsigned>::const_iterator it = sortedMatchedAK5Jets.begin(); it != sortedMatchedAK5Jets.end(); ++it)
          sortedMatchedAK5JetsIdx.push_back(it->second);
      }
      else if (matchedAK5Jets.size()==1 )
        sortedMatchedAK5JetsIdx.push_back(0);


      // find the closest b hadron to the fat jet, the two subjets, and the matched AK5 jets
      double mindRfatjet = 999., mindRsubjet1 = 999., mindRsubjet2 = 999., mindRak5jet = 999., mindRak5jet1 = 999., mindRak5jet2 = 999.;
      reco::GenParticleCollection::const_iterator bHadronMatchSubjet1 = genParticles->end(), bHadronMatchSubjet2 = genParticles->end(),
                                                  bHadronMatchAK5jet1 = genParticles->end(), bHadronMatchAK5jet2 = genParticles->end();

      for(reco::GenParticleCollection::const_iterator gpIt = genParticles->begin(); gpIt != genParticles->end(); ++gpIt)
      {
        int id = abs(gpIt->pdgId());
        // skip GenParticle if not b hadron
        if ( !((id/100)%10 == 5 || (id/1000)%10 == 5) ) continue;

        double dRfatjet = reco::deltaR( it->p4(), gpIt->p4() );
        if( dRfatjet < mindRfatjet )
          mindRfatjet = dRfatjet;

        if( subjets.size()>1 )
        {
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

        if( matchedAK5Jets.size()>0 )
        {
          double dRak5jet = reco::deltaR( matchedAK5Jets.at(sortedMatchedAK5JetsIdx.back())->p4(), gpIt->p4() );
          if( dRak5jet < mindRak5jet )
            mindRak5jet = dRak5jet;
        }

        if( matchedAK5Jets.size()>1 )
        {
          double dRak5jet1 = reco::deltaR( matchedAK5Jets.at(sortedMatchedAK5JetsIdx.back())->p4(), gpIt->p4() );
          double dRak5jet2 = reco::deltaR( matchedAK5Jets.at(sortedMatchedAK5JetsIdx.at(matchedAK5Jets.size()-2))->p4(), gpIt->p4() );
          if( dRak5jet1 < mindRak5jet1 )
          {
            mindRak5jet1 = dRak5jet1;
            bHadronMatchAK5jet1 = gpIt;
          }
          if( dRak5jet2 < mindRak5jet2 )
          {
            mindRak5jet2 = dRak5jet2;
            bHadronMatchAK5jet2 = gpIt;
          }
        }
      }

      h2_JetPt_mindRjetBhadron_BosonMatched->    Fill(jetPt, (mindRfatjet<999.  ? mindRfatjet  : -99.), eventWeight);
      h2_JetPt_mindRSubjet1Bhadron_BosonMatched->Fill(jetPt, (mindRsubjet1<999. ? mindRsubjet1 : -99.), eventWeight);
      h2_JetPt_mindRSubjet2Bhadron_BosonMatched->Fill(jetPt, (mindRsubjet2<999. ? mindRsubjet2 : -99.), eventWeight);
      if( useAK5Jets )
      {
        h2_JetPt_mindRak5jetBhadron_BosonMatched-> Fill(jetPt, (mindRak5jet<999.  ? mindRak5jet  : -99.), eventWeight);
        h2_JetPt_mindRak5jet1Bhadron_BosonMatched->Fill(jetPt, (mindRak5jet1<999. ? mindRak5jet1 : -99.), eventWeight);
        h2_JetPt_mindRak5jet2Bhadron_BosonMatched->Fill(jetPt, (mindRak5jet2<999. ? mindRak5jet2 : -99.), eventWeight);
      }

      double dRsubjets = -999.;
      if( subjets.size()>1 )
      {
        dRsubjets = reco::deltaR( subjets.at(0)->p4(), subjets.at(1)->p4() );

        h2_JetPt_dRsubjets_BosonMatched->Fill(jetPt, dRsubjets, eventWeight);
        if( bHadronMatchSubjet1 != genParticles->end() && bHadronMatchSubjet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadron_BosonMatched->Fill(jetPt, ( bHadronMatchSubjet1 == bHadronMatchSubjet2 ? 1. : 0. ), eventWeight);
      }

      if( matchedAK5Jets.size()>1 )
      {
        h2_JetPt_dRak5jets_BosonMatched->Fill(jetPt, reco::deltaR( matchedAK5Jets.at(sortedMatchedAK5JetsIdx.back())->p4(), matchedAK5Jets.at(sortedMatchedAK5JetsIdx.at(matchedAK5Jets.size()-2))->p4() ), eventWeight);
        if( bHadronMatchAK5jet1 != genParticles->end() && bHadronMatchAK5jet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched->Fill(jetPt, ( bHadronMatchAK5jet1 == bHadronMatchAK5jet2 ? 1. : 0. ), eventWeight);
      }

      // get b-tag discriminators
      double jet_CSV_discr = it->bDiscriminator((bdiscriminator).c_str());
      double jet_JP_discr = it->bDiscriminator("jetProbabilityBJetTags");
      double subJet1_CSV_discr = -999., subJet2_CSV_discr = -999.;
      double subJet1_JP_discr = -999., subJet2_JP_discr = -999.;
      double ak5Jet_CSV_discr = -999., ak5Jet1_CSV_discr = -999., ak5Jet2_CSV_discr = -999.;
      if( subjets.size()>1 )
      {
        subJet1_CSV_discr = subjets.at(0)->bDiscriminator(bdiscriminator.c_str());
        subJet2_CSV_discr = subjets.at(1)->bDiscriminator(bdiscriminator.c_str());

        subJet1_JP_discr = subjets.at(0)->bDiscriminator("jetProbabilityBJetTags");
        subJet2_JP_discr = subjets.at(1)->bDiscriminator("jetProbabilityBJetTags");
      }
      if( matchedAK5Jets.size()>0 )
        ak5Jet_CSV_discr = matchedAK5Jets.at(sortedMatchedAK5JetsIdx.back())->bDiscriminator(bdiscriminator.c_str());
      if( matchedAK5Jets.size()>1 )
      {
        ak5Jet1_CSV_discr = matchedAK5Jets.at(sortedMatchedAK5JetsIdx.back())->bDiscriminator(bdiscriminator.c_str());
        ak5Jet2_CSV_discr = matchedAK5Jets.at(sortedMatchedAK5JetsIdx.at(matchedAK5Jets.size()-2))->bDiscriminator(bdiscriminator.c_str());
      }
      double minSubJet_CSV_discr = std::min(subJet1_CSV_discr, subJet2_CSV_discr);
      double maxSubJet_CSV_discr = std::max(subJet1_CSV_discr, subJet2_CSV_discr);
      double minSubJet_JP_discr = std::min(subJet1_JP_discr, subJet2_JP_discr);
      double maxSubJet_JP_discr = std::max(subJet1_JP_discr, subJet2_JP_discr);
      double minAK5Jets_CSV_discr = std::min(ak5Jet1_CSV_discr, ak5Jet2_CSV_discr);
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

      h2_JetPt_mindRjetBhadron_BosonMatched_JetMass->    Fill(jetPt, (mindRfatjet<999.  ? mindRfatjet  : -99.), eventWeight);
      h2_JetPt_mindRSubjet1Bhadron_BosonMatched_JetMass->Fill(jetPt, (mindRsubjet1<999. ? mindRsubjet1 : -99.), eventWeight);
      h2_JetPt_mindRSubjet2Bhadron_BosonMatched_JetMass->Fill(jetPt, (mindRsubjet2<999. ? mindRsubjet2 : -99.), eventWeight);
      if( useAK5Jets )
      {
        h2_JetPt_mindRak5jetBhadron_BosonMatched_JetMass-> Fill(jetPt, (mindRak5jet<999.  ? mindRak5jet  : -99.), eventWeight);
        h2_JetPt_mindRak5jet1Bhadron_BosonMatched_JetMass->Fill(jetPt, (mindRak5jet1<999. ? mindRak5jet1 : -99.), eventWeight);
        h2_JetPt_mindRak5jet2Bhadron_BosonMatched_JetMass->Fill(jetPt, (mindRak5jet2<999. ? mindRak5jet2 : -99.), eventWeight);
      }

//       reco::TaggingVariableList tvlIP = it->tagInfoTrackIP("impactParameter")->taggingVariables();
//
//       std::cout << ">> Fat jet pt, eta, phi, mass: " << it->pt()  << ", "
//                                                      << it->eta() << ", "
//                                                      << it->phi() << ", "
//                                                      << jetMass   << std::endl;
//       std::cout << "   Assoc. tracks: " << it->associatedTracks().size() << std::endl
//                 << "   Assoc. tracks (IPTagInfo): " << it->tagInfoTrackIP("impactParameter")->tracks().size() << std::endl
//                 << "   Assoc. tracks (SVTagInfo): " << it->tagInfoSecondaryVertex("secondaryVertex")->tracks().size() << std::endl
//                 << "   Selected tracks (IPTagInfo): " << it->tagInfoTrackIP("impactParameter")->selectedTracks().size() << std::endl
//                 << "   Selected tracks (SVTagInfo): " << it->tagInfoSecondaryVertex("secondaryVertex")->selectedTracks().size() << std::endl
//                 << "   Vertex tracks (SVTagInfo): " << it->tagInfoSecondaryVertex("secondaryVertex")->vertexTracks().size() << std::endl
//                 << "   TaggingVariableList tracks: " << tvlIP.getList(reco::btau::trackDeltaR,false).size() << std::endl; // the same as selected tracks

      if( subjets.size()>1 )
      {
        h2_JetPt_dRsubjets_BosonMatched_JetMass->Fill(jetPt, dRsubjets, eventWeight);
        if( bHadronMatchSubjet1 != genParticles->end() && bHadronMatchSubjet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadron_BosonMatched_JetMass->Fill(jetPt, ( bHadronMatchSubjet1 == bHadronMatchSubjet2 ? 1. : 0. ), eventWeight);

        int nTotalFat = 0, nSharedFat = 0, nVertexTotalFat = 0, nVertexSharedFat = 0;
        GlobalPoint vertexPosition(PVs->front().x(),PVs->front().y(),PVs->front().z()); // PV position

        for(int i=0; i<2; ++i)
        {
           h1_SubJetPt_BosonMatched_JetMass->Fill(subjets.at(i)->pt(), eventWeight);

           int j = (i==0 ? 1 : 0); // companion subjet index
           int nTracks = ( subjets.at(i)->hasTagInfo("impactParameter") ? subjets.at(i)->tagInfoTrackIP("impactParameter")->selectedTracks().size() : 0 );
           int nTotal = 0, nShared = 0;

           GlobalVector direction_i(subjets.at(i)->px(), subjets.at(i)->py(), subjets.at(i)->pz()); // subjet_i direction
           GlobalVector direction_j(subjets.at(j)->px(), subjets.at(j)->py(), subjets.at(j)->pz()); // subjet_j direction

           for(int t=0; t<nTracks; ++t)
           {
             reco::TransientTrack transientTrack = trackBuilder->build(*(subjets.at(i)->tagInfoTrackIP("impactParameter")->selectedTracks().at(t))); // convert track into transient track
             // calculate decay length wrt subjet_i
             double decayLength_i = 9999.;
             TrajectoryStateOnSurface closest_i = IPTools::closestApproachToJet(transientTrack.impactPointState(), PVs->front(), direction_i, transientTrack.field());
             if (closest_i.isValid())
               decayLength_i = (closest_i.globalPosition() - vertexPosition).mag();
             // calculate distance to subjet_i axis
             double distJetAxis_i = IPTools::jetTrackDistance(transientTrack, direction_i, PVs->front()).second.value();

             if( reco::deltaR( subjets.at(i)->tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->eta(), subjets.at(i)->tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->phi(), subjets.at(i)->eta(), subjets.at(i)->phi() ) < 0.3 &&
                 decayLength_i < 5. &&
                 std::fabs(distJetAxis_i) < 0.07
             )
             {
               ++nTotal;
               ++nTotalFat;

               // calculate decay length wrt subjet_j
               double decayLength_j = 9999.;
               TrajectoryStateOnSurface closest_j = IPTools::closestApproachToJet(transientTrack.impactPointState(), PVs->front(), direction_j, transientTrack.field());
               if (closest_j.isValid())
                 decayLength_j = (closest_j.globalPosition() - vertexPosition).mag();
               // calculate distance to subjet_j axis
               double distJetAxis_j = IPTools::jetTrackDistance(transientTrack, direction_j, PVs->front()).second.value();

               if( reco::deltaR( subjets.at(i)->tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->eta(), subjets.at(i)->tagInfoTrackIP("impactParameter")->selectedTracks().at(t)->phi(), subjets.at(j)->eta(), subjets.at(j)->phi() ) < 0.3 &&
                   decayLength_j < 5. &&
                   std::fabs(distJetAxis_j) < 0.07
               )
               {
                 ++nShared;
                 if(i==0) ++nSharedFat;
               }
             }
           }

           p1_SubJetPt_TotalTracks_BosonMatched_JetMass->Fill(subjets.at(i)->pt(), nTotal, eventWeight);
           p1_SubJetPt_SharedTracks_BosonMatched_JetMass->Fill(subjets.at(i)->pt(), nShared, eventWeight);
           if( nTotal>0 ) p1_SubJetPt_SharedTracksRatio_BosonMatched_JetMass->Fill(subjets.at(i)->pt(), double(nShared)/double(nTotal), eventWeight);

           int nVertexTracks_i = ( subjets.at(i)->hasTagInfo("secondaryVertex") ? subjets.at(i)->tagInfoSecondaryVertex("secondaryVertex")->vertexTracks().size() : 0 );
           int nVertexTracks_j = ( subjets.at(j)->hasTagInfo("secondaryVertex") ? subjets.at(j)->tagInfoSecondaryVertex("secondaryVertex")->vertexTracks().size() : 0 );
           int nVertexTotal = 0, nVertexShared = 0;

           for(int vt_i=0; vt_i<nVertexTracks_i; ++vt_i)
           {
             ++nVertexTotal;
             ++nVertexTotalFat;

             for(int vt_j=0; vt_j<nVertexTracks_j; ++vt_j)
             {
               if( &(*(subjets.at(i)->tagInfoSecondaryVertex("secondaryVertex")->vertexTracks().at(vt_i))) == &(*(subjets.at(j)->tagInfoSecondaryVertex("secondaryVertex")->vertexTracks().at(vt_j))) )
               {
                 ++nVertexShared;
                 if(i==0) ++nVertexSharedFat;
               }
             }
           }

           p1_SubJetPt_TotalVertexTracks_BosonMatched_JetMass->Fill(subjets.at(i)->pt(), nVertexTotal, eventWeight);
           p1_SubJetPt_SharedVertexTracks_BosonMatched_JetMass->Fill(subjets.at(i)->pt(), nVertexShared, eventWeight);
           if( nVertexTotal>0 ) p1_SubJetPt_SharedVertexTracksRatio_BosonMatched_JetMass->Fill(subjets.at(i)->pt(), double(nVertexShared)/double(nVertexTotal), eventWeight);
        }

        p1_JetPt_TotalTracks_BosonMatched_JetMass->Fill(jetPt, nTotalFat-nSharedFat, eventWeight);
        p1_JetPt_SharedTracks_BosonMatched_JetMass->Fill(jetPt, nSharedFat, eventWeight);
        if( nTotalFat>0 ) p1_JetPt_SharedTracksRatio_BosonMatched_JetMass->Fill(jetPt, double(nSharedFat)/double(nTotalFat-nSharedFat), eventWeight);

        p1_JetPt_TotalVertexTracks_BosonMatched_JetMass->Fill(jetPt, nVertexTotalFat-nVertexSharedFat, eventWeight);
        p1_JetPt_SharedVertexTracks_BosonMatched_JetMass->Fill(jetPt, nVertexSharedFat, eventWeight);
        if( nVertexTotalFat>0 ) p1_JetPt_SharedVertexTracksRatio_BosonMatched_JetMass->Fill(jetPt, double(nVertexSharedFat)/double(nVertexTotalFat-nVertexSharedFat), eventWeight);

        p1_dRsubjets_TotalTracks_BosonMatched_JetMass->Fill(dRsubjets, nTotalFat-nSharedFat, eventWeight);
        p1_dRsubjets_SharedTracks_BosonMatched_JetMass->Fill(dRsubjets, nSharedFat, eventWeight);
        if( nTotalFat>0 ) p1_dRsubjets_SharedTracksRatio_BosonMatched_JetMass->Fill(dRsubjets, double(nSharedFat)/double(nTotalFat-nSharedFat), eventWeight);

        if( eventDisplayPrintout )
        {
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
      }

      if( matchedAK5Jets.size()>1 )
      {
        h2_JetPt_dRak5jets_BosonMatched_JetMass->Fill(jetPt, reco::deltaR( matchedAK5Jets.at(sortedMatchedAK5JetsIdx.back())->p4(), matchedAK5Jets.at(sortedMatchedAK5JetsIdx.at(matchedAK5Jets.size()-2))->p4() ), eventWeight);
        if( bHadronMatchAK5jet1 != genParticles->end() && bHadronMatchAK5jet2 != genParticles->end() )
          h2_JetPt_SameMatchedBhadronAK5Jets_BosonMatched_JetMass->Fill(jetPt, ( bHadronMatchAK5jet1 == bHadronMatchAK5jet2 ? 1. : 0. ), eventWeight);
      }

      h1_JetCSVDiscr_BosonMatched_JetMass->Fill( jet_CSV_discr, eventWeight);
      h1_SubJetMinCSVDiscr_BosonMatched_JetMass->Fill( minSubJet_CSV_discr, eventWeight);
      h1_SubJetMaxCSVDiscr_BosonMatched_JetMass->Fill( maxSubJet_CSV_discr, eventWeight);
      h1_JetDoubleBDiscr_BosonMatched_JetMass->Fill( jet_DoubleB_discr, eventWeight);

      h2_JetPt_JetCSV_BosonMatched_JetMass->Fill(jetPt, jet_CSV_discr, eventWeight);;
      h2_JetPt_SubJetMinCSV_BosonMatched_JetMass->Fill(jetPt, minSubJet_CSV_discr, eventWeight);
      h2_JetPt_SubJetMaxCSV_BosonMatched_JetMass->Fill(jetPt, maxSubJet_CSV_discr, eventWeight);

      h2_JetPt_JetJP_BosonMatched_JetMass->Fill(jetPt, jet_JP_discr, eventWeight);;
      h2_JetPt_SubJetMinJP_BosonMatched_JetMass->Fill(jetPt, minSubJet_JP_discr, eventWeight);
      h2_JetPt_SubJetMaxJP_BosonMatched_JetMass->Fill(jetPt, maxSubJet_JP_discr, eventWeight);

      h2_JetPt_AK5JetCSV_BosonMatched_JetMass->Fill(jetPt, ak5Jet_CSV_discr, eventWeight);;
      h2_JetPt_AK5JetMinCSV_BosonMatched_JetMass->Fill(jetPt, minAK5Jets_CSV_discr, eventWeight);

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

      if(dRsubjets>=0. && dRsubjets<0.2)       h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0to0p2->Fill(subJet1_CSV_discr, subJet2_CSV_discr, eventWeight);
      else if(dRsubjets>=0.2 && dRsubjets<0.4) h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p2to0p4->Fill(subJet1_CSV_discr, subJet2_CSV_discr, eventWeight);
      else if(dRsubjets>=0.4 && dRsubjets<0.6) h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p4to0p6->Fill(subJet1_CSV_discr, subJet2_CSV_discr, eventWeight);
      else if(dRsubjets>=0.6 && dRsubjets<0.8) h2_SubJet1CSV_SubJet2CSV_BosonMatched_JetMass_dRsubjets0p6to0p8->Fill(subJet1_CSV_discr, subJet2_CSV_discr, eventWeight);


      double tau1 = it->userFloat("NjettinessCA8:tau1");
      double tau2 = it->userFloat("NjettinessCA8:tau2");
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

      //h2_JetMass_nTracks_Pt[suffix]->Fill(jetMass, nTracks, eventWeight);
      //h2_JetMass_nSelectedTracks_Pt[suffix]->Fill(jetMass, nSelectedTracks, eventWeight);
      //h2_JetMass_tau2tau1_Pt[suffix]->Fill(jetMass, tau2overtau1, eventWeight);
      //h2_JetMass_SubJetMinCSV_Pt[suffix]->Fill(jetMass, minSubJet_CSV_discr, eventWeight);
      //h2_JetMass_SubJetMaxCSV_Pt[suffix]->Fill(jetMass, maxSubJet_CSV_discr, eventWeight);
      //h2_JetMass_TrackJetWidth_Pt[suffix]->Fill(jetMass, trackJetWidth, eventWeight);
      //h2_JetMass_SelectedTrackJetWidth_Pt[suffix]->Fill(jetMass, selectedTrackJetWidth, eventWeight);
      //h2_JetMass_maxdRTracks_Pt[suffix]->Fill(jetMass, maxdRTracks, eventWeight);
      //h2_JetMass_maxdRSelectedTracks_Pt[suffix]->Fill(jetMass, maxdRSelectedTracks, eventWeight);

      //h2_nTracks_tau2tau1_Pt[suffix]->Fill(nTracks, tau2overtau1, eventWeight);
      //h2_nSelectedTracks_tau2tau1_Pt[suffix]->Fill(nSelectedTracks, tau2overtau1, eventWeight);

      //h2_nTracks_SubJetMinCSV_Pt[suffix]->Fill(nTracks, minSubJet_CSV_discr, eventWeight);
      //h2_nSelectedTracks_SubJetMinCSV_Pt[suffix]->Fill(nSelectedTracks, minSubJet_CSV_discr, eventWeight);
      //h2_tau2tau1_SubJetMinCSV_Pt[suffix]->Fill(tau2overtau1, minSubJet_CSV_discr, eventWeight);
      for(unsigned i=0; i<jetPtBins; ++i)
      {
        if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
        {
          suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
          h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
          h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
          h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);

          //h2_JetMass_nTracks_Pt[suffix]->Fill(jetMass, nTracks, eventWeight);
          //h2_JetMass_nSelectedTracks_Pt[suffix]->Fill(jetMass, nSelectedTracks, eventWeight);
          //h2_JetMass_tau2tau1_Pt[suffix]->Fill(jetMass, tau2overtau1, eventWeight);
          //h2_JetMass_SubJetMinCSV_Pt[suffix]->Fill(jetMass, minSubJet_CSV_discr, eventWeight);
          //h2_JetMass_SubJetMaxCSV_Pt[suffix]->Fill(jetMass, maxSubJet_CSV_discr, eventWeight);
          //h2_JetMass_TrackJetWidth_Pt[suffix]->Fill(jetMass, trackJetWidth, eventWeight);
          //h2_JetMass_SelectedTrackJetWidth_Pt[suffix]->Fill(jetMass, selectedTrackJetWidth, eventWeight);
          //h2_JetMass_maxdRTracks_Pt[suffix]->Fill(jetMass, maxdRTracks, eventWeight);
          //h2_JetMass_maxdRSelectedTracks_Pt[suffix]->Fill(jetMass, maxdRSelectedTracks, eventWeight);

          //h2_nTracks_tau2tau1_Pt[suffix]->Fill(nTracks, tau2overtau1, eventWeight);
          //h2_nSelectedTracks_tau2tau1_Pt[suffix]->Fill(nSelectedTracks, tau2overtau1, eventWeight);

          //h2_nTracks_SubJetMinCSV_Pt[suffix]->Fill(nTracks, minSubJet_CSV_discr, eventWeight);
          //h2_nSelectedTracks_SubJetMinCSV_Pt[suffix]->Fill(nSelectedTracks, minSubJet_CSV_discr, eventWeight);
          //h2_tau2tau1_SubJetMinCSV_Pt[suffix]->Fill(tau2overtau1, minSubJet_CSV_discr, eventWeight);
        }
      }
      if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
      {
        suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
        h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
        h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
        h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, tau2overtau1, eventWeight);

        //h2_JetMass_nTracks_Pt[suffix]->Fill(jetMass, nTracks, eventWeight);
        //h2_JetMass_nSelectedTracks_Pt[suffix]->Fill(jetMass, nSelectedTracks, eventWeight);
        //h2_JetMass_tau2tau1_Pt[suffix]->Fill(jetMass, tau2overtau1, eventWeight);
        //h2_JetMass_SubJetMinCSV_Pt[suffix]->Fill(jetMass, minSubJet_CSV_discr, eventWeight);
        //h2_JetMass_SubJetMaxCSV_Pt[suffix]->Fill(jetMass, maxSubJet_CSV_discr, eventWeight);
        //h2_JetMass_TrackJetWidth_Pt[suffix]->Fill(jetMass, trackJetWidth, eventWeight);
        //h2_JetMass_SelectedTrackJetWidth_Pt[suffix]->Fill(jetMass, selectedTrackJetWidth, eventWeight);
        //h2_JetMass_maxdRTracks_Pt[suffix]->Fill(jetMass, maxdRTracks, eventWeight);
        //h2_JetMass_maxdRSelectedTracks_Pt[suffix]->Fill(jetMass, maxdRSelectedTracks, eventWeight);

        //h2_nTracks_tau2tau1_Pt[suffix]->Fill(nTracks, tau2overtau1, eventWeight);
        //h2_nSelectedTracks_tau2tau1_Pt[suffix]->Fill(nSelectedTracks, tau2overtau1, eventWeight);

        //h2_nTracks_SubJetMinCSV_Pt[suffix]->Fill(nTracks, minSubJet_CSV_discr, eventWeight);
        //h2_nSelectedTracks_SubJetMinCSV_Pt[suffix]->Fill(nSelectedTracks, minSubJet_CSV_discr, eventWeight);
        //h2_tau2tau1_SubJetMinCSV_Pt[suffix]->Fill(tau2overtau1, minSubJet_CSV_discr, eventWeight);
      }

      // mass drop
      if( calculateMassDrop )
      {
        double fatJetMass = jetMass;
        double subjetMass = ( subjets.size()>1 ? std::max( subjets.at(0)->mass(), subjets.at(1)->mass() ) : 0. );
        if( useUncorrMassForMassDrop && subjets.size()>1 )
        {
          fatJetMass = groomedBasicJetMatch->correctedJet("Uncorrected").mass();
          subjetMass = std::max( subjets.at(0)->correctedJet("Uncorrected").mass(), subjets.at(1)->correctedJet("Uncorrected").mass() );
        }
        double massDrop = ( jetMass>0. ? subjetMass/fatJetMass : -10.);
        // fill nPV_MassDrop histograms
        suffix = Form("%.0ftoInf",jetPtMin);
        h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
        //h2_JetMass_MassDrop_Pt[suffix]->Fill(jetMass, massDrop, eventWeight);
        for(unsigned i=0; i<jetPtBins; ++i)
        {
          if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
          {
            suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
            h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
            //h2_JetMass_MassDrop_Pt[suffix]->Fill(jetMass, massDrop, eventWeight);
          }
        }
        if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
        {
          suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
          h2_nPV_MassDrop_Pt[suffix]->Fill(nPV, massDrop, eventWeight);
          //h2_JetMass_MassDrop_Pt[suffix]->Fill(jetMass, massDrop, eventWeight);
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
