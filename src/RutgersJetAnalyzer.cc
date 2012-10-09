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
// $Id: RutgersJetAnalyzer.cc,v 1.7.2.2 2012/10/08 19:23:14 ferencek Exp $
//
//


// system include files
#include <memory>

// FastJet include files
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
// N-subjettiness include files
#include "RutgersSandbox/RutgersJetAnalyzer/plugins/NjettinessPlugin.hh"
#include "RutgersSandbox/RutgersJetAnalyzer/plugins/Nsubjettiness.hh"
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
    typedef boost::shared_ptr<fastjet::ClusterSequence>        ClusterSequencePtr;
    typedef boost::shared_ptr<fastjet::JetDefinition::Plugin>  PluginPtr;
    typedef boost::shared_ptr<fastjet::JetDefinition>          JetDefPtr;
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
    const edm::InputTag genParticleTag;
    const edm::InputTag jetsTag;
    const edm::InputTag groomedJetsTag;
    const edm::InputTag pvTag;
    const double        fjInputPtMin;       // minimum pT of input constituents clustered by FastJet
    const double        fjJetPtMin;         // minimum jet pT returned by FastJet
    const double	fjRadius;	    // radius for jet clustering
    const bool          doWMatching;	    // parameter for deciding if matching on or off
    const double	wMatchingRadius;
    const double	leptonMatchingRadius;
    const double	jetPtMin;
    const double	jetAbsEtaMax;
    const double	jetMassMin;
    const double	jetMassMax;
    const bool          useGroomedJets;
    const bool          useEventWeight;
    bool                useUncorrectedJets;
    JetDefPtr           jetDefinitionAK;  // Anti-kT jet definition
    JetDefPtr           jetDefinitionKT;  // kT jet definition
    edm::Service<TFileService> fs;

    TH1D *h1_nPV;

    TH1D *h1_JetPt;
    TH1D *h1_JetPt_JetMass;
    TH1D *h1_JetEta;
    TH1D *h1_JetEta_JetMass;

    TH2D *h2_nPV_JetMass_Pt300toInf;
    TH2D *h2_nPV_JetMass_Pt300to500;
    TH2D *h2_nPV_JetMass_Pt500to700;
    TH2D *h2_nPV_JetMass_Pt700to900;
    TH2D *h2_nPV_JetMass_Pt900toInf;

    TH2D *h2_nPV_tau1_Pt300toInf;
    TH2D *h2_nPV_tau1_Pt300to500;
    TH2D *h2_nPV_tau1_Pt500to700;
    TH2D *h2_nPV_tau1_Pt700to900;
    TH2D *h2_nPV_tau1_Pt900toInf;

    TH2D *h2_nPV_tau2_Pt300toInf;
    TH2D *h2_nPV_tau2_Pt300to500;
    TH2D *h2_nPV_tau2_Pt500to700;
    TH2D *h2_nPV_tau2_Pt700to900;
    TH2D *h2_nPV_tau2_Pt900toInf;

    TH2D *h2_nPV_tau2tau1_Pt300toInf;
    TH2D *h2_nPV_tau2tau1_Pt300to500;
    TH2D *h2_nPV_tau2tau1_Pt500to700;
    TH2D *h2_nPV_tau2tau1_Pt700to900;
    TH2D *h2_nPV_tau2tau1_Pt900toInf;
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

  genParticleTag(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
  jetsTag(iConfig.getParameter<edm::InputTag>("JetsTag")),
  groomedJetsTag(iConfig.getParameter<edm::InputTag>("GroomedJetsTag")),
  pvTag(iConfig.getParameter<edm::InputTag>("PvTag")),
  fjInputPtMin(iConfig.getParameter<double>("FJInputPtMin")),
  fjJetPtMin(iConfig.getParameter<double>("FJJetPtMin")),
  fjRadius(iConfig.getParameter<double>("FJRadius")),
  doWMatching(iConfig.getParameter<bool>("DoWMatching")),
  wMatchingRadius(iConfig.getParameter<double>("WMatchingRadius")),
  leptonMatchingRadius(iConfig.getParameter<double>("LeptonMatchingRadius")),
  jetPtMin(iConfig.getParameter<double>("JetPtMin")),
  jetAbsEtaMax(iConfig.getParameter<double>("JetAbsEtaMax")),
  jetMassMin(iConfig.getParameter<double>("JetMassMin")),
  jetMassMax(iConfig.getParameter<double>("JetMassMax")),
  useGroomedJets(iConfig.getParameter<bool>("UseGroomedJets")),
  useEventWeight(iConfig.getParameter<bool>("UseEventWeight")),
  useUncorrectedJets(false)

{
    //now do what ever initialization is needed
    jetDefinitionAK = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, fjRadius) );
    jetDefinitionKT = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, fjRadius) );

    if ( iConfig.exists("UseUncorrectedJets") )
      useUncorrectedJets = iConfig.getParameter<bool>("UseUncorrectedJets");

    h1_nPV = fs->make<TH1D>("h1_nPV","PV Multiplicity;nPV;",51,-0.5,50.5);

    h1_JetPt         = fs->make<TH1D>("h1_JetPt",";p_{T} [GeV];",1000,0.,1000.);
    h1_JetPt_JetMass = fs->make<TH1D>("h1_JetPt_JetMass","Jet mass cut;p_{T} [GeV];",1000,0.,1000.);
    h1_JetEta         = fs->make<TH1D>("h1_JetEta",";#eta;",120,-3.,3.);
    h1_JetEta_JetMass = fs->make<TH1D>("h1_JetEta_JetMass","Jet mass cut;#eta;",120,-3.,3.);

    h2_nPV_JetMass_Pt300toInf = fs->make<TH2D>("h2_nPV_JetMass_Pt300toInf","p_{T}>900 GeV;nPV;m_{jet} [GeV]",51,-0.5,50.5,200,0,200.);
    h2_nPV_JetMass_Pt300to500 = fs->make<TH2D>("h2_nPV_JetMass_Pt300to500","300<p_{T}<500 GeV;nPV;m_{jet} [GeV]",51,-0.5,50.5,200,0,200.);
    h2_nPV_JetMass_Pt500to700 = fs->make<TH2D>("h2_nPV_JetMass_Pt500to700","500<p_{T}<700 GeV;nPV;m_{jet} [GeV]",51,-0.5,50.5,200,0,200.);
    h2_nPV_JetMass_Pt700to900 = fs->make<TH2D>("h2_nPV_JetMass_Pt700to900","700<p_{T}<900 GeV;nPV;m_{jet} [GeV]",51,-0.5,50.5,200,0,200.);
    h2_nPV_JetMass_Pt900toInf = fs->make<TH2D>("h2_nPV_JetMass_Pt900toInf","p_{T}>900 GeV;nPV;m_{jet} [GeV]",51,-0.5,50.5,200,0,200.);

    h2_nPV_tau1_Pt300toInf = fs->make<TH2D>("h2_nPV_tau1_Pt300toInf","p_{T}>900 GeV;nPV;#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau1_Pt300to500 = fs->make<TH2D>("h2_nPV_tau1_Pt300to500","300<p_{T}<500 GeV;nPV;#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau1_Pt500to700 = fs->make<TH2D>("h2_nPV_tau1_Pt500to700","500<p_{T}<700 GeV;nPV;#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau1_Pt700to900 = fs->make<TH2D>("h2_nPV_tau1_Pt700to900","700<p_{T}<900 GeV;nPV;#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau1_Pt900toInf = fs->make<TH2D>("h2_nPV_tau1_Pt900toInf","p_{T}>900 GeV;nPV;#tau_{1}",51,-0.5,50.5,100,0,1.);

    h2_nPV_tau2_Pt300toInf = fs->make<TH2D>("h2_nPV_tau2_Pt300toInf","p_{T}>900 GeV;nPV;#tau_{2}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2_Pt300to500 = fs->make<TH2D>("h2_nPV_tau2_Pt300to500","300<p_{T}<500 GeV;nPV;#tau_{2}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2_Pt500to700 = fs->make<TH2D>("h2_nPV_tau2_Pt500to700","500<p_{T}<700 GeV;nPV;#tau_{2}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2_Pt700to900 = fs->make<TH2D>("h2_nPV_tau2_Pt700to900","700<p_{T}<900 GeV;nPV;#tau_{2}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2_Pt900toInf = fs->make<TH2D>("h2_nPV_tau2_Pt900toInf","p_{T}>900 GeV;nPV;#tau_{2}",51,-0.5,50.5,100,0,1.);

    h2_nPV_tau2tau1_Pt300toInf = fs->make<TH2D>("h2_nPV_tau2tau1_Pt300toInf","p_{T}>300 GeV;nPV;#tau_{2}/#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2tau1_Pt300to500 = fs->make<TH2D>("h2_nPV_tau2tau1_Pt300to500","300<p_{T}<500 GeV;nPV;#tau_{2}/#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2tau1_Pt500to700 = fs->make<TH2D>("h2_nPV_tau2tau1_Pt500to700","500<p_{T}<700 GeV;nPV;#tau_{2}/#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2tau1_Pt700to900 = fs->make<TH2D>("h2_nPV_tau2tau1_Pt700to900","700<p_{T}<900 GeV;nPV;#tau_{2}/#tau_{1}",51,-0.5,50.5,100,0,1.);
    h2_nPV_tau2tau1_Pt900toInf = fs->make<TH2D>("h2_nPV_tau2tau1_Pt900toInf","p_{T}>900 GeV;nPV;#tau_{2}/#tau_{1}",51,-0.5,50.5,100,0,1.);
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

    // vectors of pointers to Ws and status=3 charged leptons
    std::vector<const reco::GenParticle*> Ws, st3ChargedLeptons;

    for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it)
    {
      if( abs(it->pdgId()) == 24 ) Ws.push_back(&(*it));
      if( it->status() == 3 && ( abs(it->pdgId()) == 11 || abs(it->pdgId()) == 13 || abs(it->pdgId()) == 15 ) ) st3ChargedLeptons.push_back(&(*it));
    }

    // loop over jets
    for(PatJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it)
    {
      double jetPt = it->pt();
      // skip the jet if it does not pass pT and eta cuts
      if( !(jetPt > jetPtMin && fabs(it->eta()) < jetAbsEtaMax) ) continue;

      bool isMatched = false;

      // perform W matching with lepton veto
      if( doWMatching )
      {
        for(std::vector<const reco::GenParticle*>::const_iterator itW = Ws.begin(); itW != Ws.end(); ++itW)
	{
          if( reco::deltaR( (*itW)->p4(), it->p4() ) < wMatchingRadius )
	  {
            isMatched = true;
	    break;
	  }
	}

	for(std::vector<const reco::GenParticle*>::const_iterator itL = st3ChargedLeptons.begin(); itL != st3ChargedLeptons.end(); ++itL)
	{
          if( reco::deltaR( (*itL)->p4(), it->p4() ) < leptonMatchingRadius )
	  {
            isMatched = false;
	    break;
	  }
	}
      }
      else
        isMatched = true;

      // skip the jet if it is not matched to a W boson
      if( !isMatched ) continue;

      h1_JetPt->Fill(jetPt, eventWeight);
      h1_JetEta->Fill(it->eta(), eventWeight);

      double jetMass = it->mass();
      if( useGroomedJets )
      {
        bool matchFound = false;
	for(PatJetCollection::const_iterator itGJ = groomedJets->begin(); itGJ != groomedJets->end(); ++itGJ)
	{
          if( reco::deltaR( it->p4(), itGJ->p4() ) < fjRadius )
	  {
            matchFound = true;
	    jetMass = itGJ->mass();
	    break;
	  }
	}
	if( !matchFound ) edm::LogError("NoMatchingGroomedJet") << "Matching groomed jet not found. Using the original jet mass.";
      }

      // fill nPV_JetMass histograms
      h2_nPV_JetMass_Pt300toInf->Fill(nPV, jetMass, eventWeight);
      if( jetPt>300 && jetPt<=500 )      h2_nPV_JetMass_Pt300to500->Fill(nPV, jetMass, eventWeight);
      else if( jetPt>500 && jetPt<=700 ) h2_nPV_JetMass_Pt500to700->Fill(nPV, jetMass, eventWeight);
      else if( jetPt>700 && jetPt<=900 ) h2_nPV_JetMass_Pt700to900->Fill(nPV, jetMass, eventWeight);
      else                               h2_nPV_JetMass_Pt900toInf->Fill(nPV, jetMass, eventWeight);

      // skip the jet if it does not pass the invariant mass cut
      if( !(jetMass > jetMassMin && jetMass < jetMassMax) ) continue;

      h1_JetPt_JetMass->Fill(jetPt, eventWeight);
      h1_JetEta_JetMass->Fill(it->eta(), eventWeight);

      std::vector<fastjet::PseudoJet> fjInputs;
      std::vector<edm::Ptr<reco::PFCandidate> > constituents = it->getPFConstituents();
      std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator m;
      for ( m = constituents.begin(); m != constituents.end(); ++m )
      {
	reco::PFCandidatePtr constit = *m;
	if (constit->pt() < fjInputPtMin) continue;
	if (constit->pt() == 0)
	{
          edm::LogError("NullTransverseMomentum") << "dropping input candidate with pt=0";
	  continue;
	}
	fjInputs.push_back(fastjet::PseudoJet(constit->px(),constit->py(),constit->pz(),constit->energy()));
	fjInputs.back().set_user_index(m - constituents.begin());
      }

      ClusterSequencePtr reclusterSeq = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *jetDefinitionAK ));
      std::vector<fastjet::PseudoJet> reclusteredJets = fastjet::sorted_by_pt(reclusterSeq->inclusive_jets(fjJetPtMin));

      double beta = 1.0, R0 = fjRadius, Rcut = fjRadius;
      fastjet::Nsubjettiness nSub1OnePass_reco(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
      double tau1 = nSub1OnePass_reco(reclusteredJets[0]);
      fastjet::Nsubjettiness nSub2OnePass_reco(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
      double tau2 = nSub2OnePass_reco(reclusteredJets[0]);

      // fill nPV_tau histograms
      h2_nPV_tau1_Pt300toInf->Fill(nPV, tau1, eventWeight);
      h2_nPV_tau2_Pt300toInf->Fill(nPV, tau2, eventWeight);
      h2_nPV_tau2tau1_Pt300toInf->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
      if( jetPt>300 && jetPt<=500 )
      {
        h2_nPV_tau1_Pt300to500->Fill(nPV, tau1, eventWeight);
        h2_nPV_tau2_Pt300to500->Fill(nPV, tau2, eventWeight);
        h2_nPV_tau2tau1_Pt300to500->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
      }
      else if( jetPt>500 && jetPt<=700 )
      {
        h2_nPV_tau1_Pt500to700->Fill(nPV, tau1, eventWeight);
        h2_nPV_tau2_Pt500to700->Fill(nPV, tau2, eventWeight);
	h2_nPV_tau2tau1_Pt500to700->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
      }
      else if( jetPt>700 && jetPt<=900 )
      {
        h2_nPV_tau1_Pt700to900->Fill(nPV, tau1, eventWeight);
	h2_nPV_tau2_Pt700to900->Fill(nPV, tau2, eventWeight);
	h2_nPV_tau2tau1_Pt700to900->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
      }
      else
      {
        h2_nPV_tau1_Pt900toInf->Fill(nPV, tau1, eventWeight);
        h2_nPV_tau2_Pt900toInf->Fill(nPV, tau2, eventWeight);
	h2_nPV_tau2tau1_Pt900toInf->Fill(nPV, (tau1>0 ? tau2/tau1 : -10.), eventWeight);
      }
    }

    //############################################################################################
    //## This section demonstrates how to cluster jets from scratch
    //##
    //    edm::Handle<reco::CandidateView> inputsHandle;
    //    iEvent.getByLabel(inputsTag, inputsHandle);
    // 
    //    std::cout << "Number of GenParticles: " << inputsHandle->size() << std::endl;
    // 
    //    std::vector<edm::Ptr<reco::Candidate> > inputs; // inputs for jet clustering
    // 
    //    for (size_t i = 0; i < inputsHandle->size(); ++i) inputs.push_back(inputsHandle->ptrAt(i));
    // 
    //    std::vector<fastjet::PseudoJet> fjInputs; // FastJet inputs for jet clustering
    // 
    //    // convert candidates to fastjet::PseudoJets
    //    std::vector<edm::Ptr<reco::Candidate> >::const_iterator inBegin = inputs.begin(),
    //      inEnd = inputs.end(), i = inBegin;
    //    for (; i != inEnd; ++i )
    //    {
    //      reco::CandidatePtr input = *i;
    //      if (input->pt() < inputPtMin) continue;
    //      if (input->pt() == 0) {
    //        edm::LogError("NullTransverseMomentum") << "dropping input candidate with pt=0";
    //        continue;
    //      }
    // 
    //      fjInputs.push_back(fastjet::PseudoJet(input->px(),input->py(),input->pz(), input->energy()));
    //      fjInputs.back().set_user_index(i - inBegin);
    //    }
    // 
    //    // define clustering sequence
    //    ClusterSequencePtr clusterSeqAK = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *jetDefinitionAK ) );
    // 
    //    // get clustered jets
    //    std::vector<fastjet::PseudoJet> jetsAK = fastjet::sorted_by_pt(clusterSeqAK->inclusive_jets(jetPtMin));
    // 
    //    std::cout << "Number of GenJets: " << jetsAK.size() << std::endl;
    //##
    //############################################################################################
    
    //############################################################################################
    //## This section demonstrates how to initialize N-subjettiness
    //##
    //    double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
    //    double R0 = 1.2; // Characteristic jet radius for normalization
    //    double Rcut = 1.4; // maximum R particles can be from axis to be included in jet
    //    fastjet::Nsubjettiness nSub2OnePass(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
    //    JetDefPtr nsubOnepass_jetDef = JetDefPtr( new fastjet::JetDefinition(new fastjet::NjettinessPlugin(3, Njettiness::onepass_kt_axes, beta, R0, Rcut)) );
    //##
    //############################################################################################
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
