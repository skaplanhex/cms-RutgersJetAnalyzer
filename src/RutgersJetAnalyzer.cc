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
// $Id: RutgersJetAnalyzer.cc,v 1.7.2.5 2012/11/01 03:48:59 ferencek Exp $
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
    const edm::InputTag pvTag;
    const double	jetRadius;	    // radius for jet clustering
    const bool          doBosonMatching;    // parameter for deciding if matching is on or off
    const double	bosonMatchingRadius;
    const int    	bosonPdgId;
    const bool          doBosonDecayProdSelection;
    const std::vector<int> bosonDecayProdPdgIds;
    const bool          useMassDrop;
    const double	jetPtMin;
    const unsigned      jetPtBins;
    const double        jetPtBinWidth;
    const double	jetAbsEtaMax;
    const double	jetMassMin;
    const double	jetMassMax;
    bool                useUncorrectedJets;

    Njettiness nsubjettinessCalculator;

    edm::Service<TFileService> fs;

    TH1D *h1_nPV;

    TH1D *h1_BosonPt;
    TH1D *h1_BosonEta;
    TH1D *h1_BosonPt_DecaySel;
    TH1D *h1_BosonEta_DecaySel;
    TH1D *h1_BosonPt_Matched;
    TH1D *h1_BosonPt_DecayProdMatched;

    TH2D *h2_BosonPt_dRdecay;

    TH1D *h1_JetPt;
    TH1D *h1_JetPt_BosonMatched;
    TH1D *h1_JetPt_BosonMatched_JetMass;
    TH1D *h1_JetPt_BosonDecayProdMatched;
    TH1D *h1_JetPt_BosonDecayProdMatched_JetMass;
    TH1D *h1_JetEta;
    TH1D *h1_JetEta_BosonMatched;
    TH1D *h1_JetEta_BosonMatched_JetMass;

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
  pvTag(iConfig.getParameter<edm::InputTag>("PvTag")),
  jetRadius(iConfig.getParameter<double>("JetRadius")),
  doBosonMatching(iConfig.getParameter<bool>("DoBosonMatching")),
  bosonMatchingRadius(iConfig.getParameter<double>("BosonMatchingRadius")),
  bosonPdgId(iConfig.getParameter<int>("BosonPdgId")),
  doBosonDecayProdSelection(iConfig.getParameter<bool>("DoBosonDecayProdSelection")),
  bosonDecayProdPdgIds(iConfig.getParameter<std::vector<int> >("BosonDecayProdPdgIds")),
  useMassDrop(iConfig.getParameter<bool>("UseMassDrop")),
  jetPtMin(iConfig.getParameter<double>("JetPtMin")),
  jetPtBins(iConfig.getParameter<unsigned>("JetPtBins")),
  jetPtBinWidth(iConfig.getParameter<double>("JetPtBinWidth")),
  jetAbsEtaMax(iConfig.getParameter<double>("JetAbsEtaMax")),
  jetMassMin(iConfig.getParameter<double>("JetMassMin")),
  jetMassMax(iConfig.getParameter<double>("JetMassMax")),
  useUncorrectedJets(false),
  nsubjettinessCalculator(Njettiness::onepass_kt_axes, NsubParameters(1.0, jetRadius, jetRadius))

{
    //now do what ever initialization is needed
    if ( iConfig.exists("UseUncorrectedJets") )
      useUncorrectedJets = iConfig.getParameter<bool>("UseUncorrectedJets");

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

    h1_BosonPt           = fs->make<TH1D>("h1_BosonPt",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_BosonEta          = fs->make<TH1D>("h1_BosonEta",";#eta;",etaBins,etaMin,etaMax);
    h1_BosonPt_DecaySel  = fs->make<TH1D>("h1_BosonPt_DecaySel",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_BosonEta_DecaySel = fs->make<TH1D>("h1_BosonEta_DecaySel",";#eta;",etaBins,etaMin,etaMax);
    h1_BosonPt_Matched  = fs->make<TH1D>("h1_BosonPt_Matched",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_BosonPt_DecayProdMatched  = fs->make<TH1D>("h1_BosonPt_DecayProdMatched",";p_{T} [GeV];",ptBins,ptMin,ptMax);

    h2_BosonPt_dRdecay = fs->make<TH2D>("h2_BosonPt_dRdecay",";p_{T} [GeV];#DeltaR",ptBins,ptMin,ptMax,dRBins,dRMin,dRMax);

    h1_JetPt = fs->make<TH1D>("h1_JetPt",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched = fs->make<TH1D>("h1_JetPt_BosonMatched",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonMatched_JetMass  = fs->make<TH1D>("h1_JetPt_BosonMatched_JetMass","Jet mass cut;p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonDecayProdMatched = fs->make<TH1D>("h1_JetPt_BosonDecayProdMatched",";p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetPt_BosonDecayProdMatched_JetMass  = fs->make<TH1D>("h1_JetPt_BosonDecayProdMatched_JetMass","Jet mass cut;p_{T} [GeV];",ptBins,ptMin,ptMax);
    h1_JetEta = fs->make<TH1D>("h1_JetEta",";#eta;",etaBins,etaMin,etaMax);
    h1_JetEta_BosonMatched = fs->make<TH1D>("h1_JetEta_BosonMatched",";#eta;",etaBins,etaMin,etaMax);
    h1_JetEta_BosonMatched_JetMass = fs->make<TH1D>("h1_JetEta_BosonMatched_JetMass","Jet mass cut;#eta;",etaBins,etaMin,etaMax);

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

    // vector of pointers to bosons
    std::vector<const reco::GenParticle*> bosons;
    // map to vectors of pointers to boson decay products
    std::map<const reco::GenParticle*,std::vector<const reco::Candidate*> > decayProducts;

    // loop over GenParticles and select bosons
    for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
    {
      if( abs(it->pdgId()) == abs(bosonPdgId) && it->status() == 3 )
      {
        h1_BosonPt->Fill( it->pt(), eventWeight );
	h1_BosonEta->Fill( it->eta(), eventWeight );

	if( doBosonDecayProdSelection )
	{
	  bool decayProductsFound = false;

	  for(unsigned i=0; i<it->numberOfDaughters(); ++i)
	  {
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
      if( useGroomedJets )
      {
        bool matchFound = false;
	for(PatJetCollection::const_iterator gjIt = groomedJets->begin(); gjIt != groomedJets->end(); ++gjIt)
	{
          if( reco::deltaR( it->p4(), gjIt->p4() ) < jetRadius )
	  {
            matchFound = true;
	    jetMass = gjIt->mass();
	    break;
	  }
	}
	if( !matchFound ) edm::LogError("NoMatchingGroomedJet") << "Matching groomed jet not found. Using the original jet mass.";
      }

      bool isMatched = false;

      // perform boson matching
      if( doBosonMatching )
      {
        for(std::vector<const reco::GenParticle*>::const_iterator bosonIt = bosons.begin(); bosonIt != bosons.end(); ++bosonIt)
	{
          if( reco::deltaR( (*bosonIt)->p4(), it->p4() ) < bosonMatchingRadius )
	  {
            isMatched = true;
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
        isMatched = true;

      // skip the jet if it is not matched to a boson
      if( !isMatched ) continue;

      h1_JetPt_BosonMatched->Fill(jetPt, eventWeight);
      h1_JetEta_BosonMatched->Fill(it->eta(), eventWeight);

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

      if( !useMassDrop )
      {
        std::vector<fastjet::PseudoJet> fjConstituents;
        std::vector<edm::Ptr<reco::PFCandidate> > constituents = it->getPFConstituents();
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

        // fill nPV_tau histograms
        suffix = Form("%.0ftoInf",jetPtMin);
        h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
        h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
        h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
        for(unsigned i=0; i<jetPtBins; ++i)
        {
          if( jetPt>(jetPtMin + jetPtBinWidth*i) && jetPt<=(jetPtMin + jetPtBinWidth*(i+1)) )
          {
            suffix = Form("%.0fto%.0f",(jetPtMin + jetPtBinWidth*i),(jetPtMin + jetPtBinWidth*(i+1)));
            h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
            h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
            h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
          }
        }
        if( jetPt>(jetPtMin+jetPtBinWidth*jetPtBins))
        {
          suffix = Form("%.0ftoInf",(jetPtMin+jetPtBinWidth*jetPtBins));
          h2_nPV_tau1_Pt[suffix]->Fill(nPV, tau1, eventWeight);
          h2_nPV_tau2_Pt[suffix]->Fill(nPV, tau2, eventWeight);
          h2_nPV_tau2tau1_Pt[suffix]->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
        }
      }
      else
      {
        double massDrop = std::max( it->daughter(0)->mass(), it->daughter(1)->mass() ) / jetMass;
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
