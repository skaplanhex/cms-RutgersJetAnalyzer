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
// $Id: RutgersJetAnalyzer.cc,v 1.3 2012/07/30 18:16:55 mzientek Exp $
//
//


// system include files
#include <memory>

// FastJet include files
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

// user include files
#include <boost/shared_ptr.hpp>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RutgersSandbox/RutgersJetAnalyzer/plugins/NjettinessPlugin.hh"
#include "RutgersSandbox/RutgersJetAnalyzer/plugins/Nsubjettiness.hh"
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
      edm::Service<TFileService> fs;
      TH1D *hT1;
      TH1D *hT2;
      TH1D *hT2byT1;
      TH1D *hpt_nocut;
      TH1D *hpt_cut;
      TH1D *heta_nocut;
      TH1D *heta_cut;
      TH1D *hmass_nocut;
      TH1D *hmass_cut;
  
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

  //hT1->Write();
  //hT2->Write();
  //hT2byT1->Write();

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
       reco::GenJet genjet = genJets->at(j);
       fastjet::PseudoJet *genjetfast = new fastjet::PseudoJet(genjet.px(),genjet.py(),genjet.pz(),genjet.energy());
         if (genjetfast->pt() > 30) {
       hpt_nocut->Fill(genjetfast->pt());
       heta_nocut->Fill(genjetfast->eta());
       hmass_nocut->Fill(genjetfast->m());
       if (115 < genjetfast->m() && genjetfast->m() < 125)
	 hpt_cut->Fill(genjetfast->pt());
       if (genjetfast->pt() > 200) {
	 hmass_cut->Fill(genjetfast->m());
	 if (115 < genjetfast->m() && genjetfast->m() < 125)
	   heta_cut->Fill(genjetfast->eta());
       }
         }
	 if(genjetfast->pt() > 200 && abs(genjetfast->eta())<2.5 && 115 < genjetfast->m() && genjetfast->m()<125){
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
	 long double T1sum=0;
	 long double d0=0;
	 double R1=1.2;
	 long double T2sum=0;
	 for(int k=0; k < (int)fjInputs.size(); k++){
	   long double dR = genjetfast->delta_R(fjInputs[k]);
	   long double dR1 = subJets[0].delta_R(fjInputs[k]);
	   long double dR2 = subJets[1].delta_R(fjInputs[k]);
	   if(dR1 < dR2)
	     T2sum += dR1*fjInputs[k].pt();
	   else
	     T2sum += dR2*fjInputs[k].pt();
	   T1sum += dR*fjInputs[k].pt();
	   d0 += fjInputs[k].pt()*R1;
	 }
	 double T1=T1sum/d0;
	 double T2=T2sum/d0;
	 double ratio = T2/T1;
	 //std:: cout << T1 << " " << T2 << " " << ratio << std::endl;
	 hT1->Fill(T1);
	 hT2->Fill(T2);
	 hT2byT1->Fill(ratio);
       }
	 } //end boosted & higgs mass conditions
         //std::cout << "dR(subjet1,subjet2) for Jet" << j << ": " << subJets[0].delta_R(subJets[1]) << std::endl;
     } // end for ()
} // end if()
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
  hT1 = fs->make<TH1D>("hT1" , "T_1" , 100 , 0 , 1 );
  hT2 = fs->make<TH1D>("hT2" , "T_2" , 100 , 0 , 1 );
  hT2byT1 = fs->make<TH1D>("hT2/T1" , "T_2/T_1" , 100 , 0 , 1 );
  hpt_nocut = fs->make<TH1D>("hpt_nocut","Jet pT with no cuts",250,0,450);
  hpt_cut = fs->make<TH1D>("hpt_cut","Jet pT with inv mass cut",250,0,450);
  heta_nocut = fs->make<TH1D>("heta_nocut","Jet eta with no cuts",450,-5,5);
  heta_cut = fs->make<TH1D>("heta_cut","Jet eta with cuts",200,-5,5);
  hmass_nocut = fs->make<TH1D>("hmass_nocut","Jet mass with no cuts",300,0,300);
  hmass_cut = fs->make<TH1D>("hmass_cut","Jet mass with cuts",300,0,300);
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
