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
// $Id$
//
//


// system include files
#include <memory>
#include <boost/shared_ptr.hpp>

// FastJet include files
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

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


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      const edm::InputTag genJetsTag;
      const edm::InputTag inputsTag;
      const double        inputPtMin;       // minimum pT of input constituents
      const double        jetPtMin;         // minimum jet pT
      JetDefPtr           jetDefinitionAK;  // Anti-kT jet definition
      JetDefPtr           jetDefinitionKT;  // kT jet definition
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
RutgersJetAnalyzer::RutgersJetAnalyzer(const edm::ParameterSet& iConfig):
   genJetsTag(iConfig.getParameter<edm::InputTag>("GenJetsTag")),
   inputsTag(iConfig.getParameter<edm::InputTag>("InputsTag")),
   inputPtMin(iConfig.getParameter<double>("InputPtMin")),
   jetPtMin(iConfig.getParameter<double>("JetPtMin"))
{
   //now do what ever initialization is needed
   jetDefinitionAK = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.2) );
   jetDefinitionKT = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, 1.4) );
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
   edm::Handle<reco::GenJetCollection> genJets;
   iEvent.getByLabel(genJetsTag, genJets);

   if ( genJets->size() >= 2 )
   {
     // loop over two leading GenJets
     for ( size_t j = 0; j < 2; ++j )
     {
       //std::cout << "Jet " << i << " has " << genJets->at(i).getJetConstituents().size() << " constituents" << std::endl;

       std::vector<edm::Ptr<reco::Candidate> > inputs; // inputs for subjet clustering
       inputs = genJets->at(j).getJetConstituents();

       std::vector<fastjet::PseudoJet> fjInputs; // FastJet inputs for jet clustering

       // convert candidates to fastjet::PseudoJets
       std::vector<edm::Ptr<reco::Candidate> >::const_iterator inBegin = inputs.begin(),
         inEnd = inputs.end(), i = inBegin;
       for ( ; i != inEnd; ++i )
       {
         reco::CandidatePtr input = *i;
         if (input->pt() < inputPtMin) continue;
         if (input->pt() == 0) {
           edm::LogError("NullTransverseMomentum") << "dropping input candidate with pt=0";
           continue;
         }

         fjInputs.push_back(fastjet::PseudoJet(input->px(),input->py(),input->pz(), input->energy()));
         fjInputs.back().set_user_index(i - inBegin);
       }
       // define clustering sequence
       ClusterSequencePtr clusterSeqKT = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *jetDefinitionKT ) );
       // kT jets
       std::vector<fastjet::PseudoJet> kTJets = fastjet::sorted_by_pt(clusterSeqKT->inclusive_jets(jetPtMin));

       // get 2 subjets only for jets that have at least 2 constituents
       if ( fjInputs.size() >= 2 && kTJets.size() >=1 ) {

         std::vector<fastjet::PseudoJet> subJets = fastjet::sorted_by_pt(kTJets[0].exclusive_subjets(2));

         std::cout << "dR(subjet1,subjet2) for Jet" << j << ": " << subJets[0].delta_R(subJets[1]) << std::endl;
       }
     } // end for ()
   } // end i f()

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
