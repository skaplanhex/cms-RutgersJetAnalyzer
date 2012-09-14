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
// $Id: RutgersJetAnalyzer.cc,v 1.2.2.5 2012/09/13 21:31:20 mzientek Exp $
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
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TH1.h"
#include "TH2.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RutgersSandbox/RutgersJetAnalyzer/plugins/Nsubjettiness.hh"
#include "RutgersSandbox/RutgersJetAnalyzer/plugins/NjettinessPlugin.hh"

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
      const edm::InputTag genParticleTag;
      const edm::InputTag inputsTag;
      const edm::InputTag jetsTag;
      const edm::InputTag pvtag;
      const double        inputPtMin;       // minimum pT of input constituents
      const double        jetPtMin;         // minimum jet pT
      const int           matching;	    // parameter for deciding if matching on or off
      const double	  radius;	    // radius for jet clustering
      JetDefPtr           jetDefinitionAK;  // Anti-kT jet definition
      JetDefPtr           jetDefinitionKT;  // kT jet definition
      edm::Service<TFileService> fs;
      //TH1D *hT1;
      //TH1D *hT2;
      //TH1D *hT2byT1;
      //TH1D *hpt_nocut;
      //TH1D *hpt_cut;
      //TH1D *heta_nocut;
      //TH1D *heta_cut;
      //TH1D *hmass_nocut;
      //TH1D *hmass_cut;
      //TH1D *hpt_nocut_cluster;
      //TH1D *heta_nocut_cluster;
      //TH1D *hmass_nocut_cluster;
      /*TH1D *h_tau1onepass_gen;
      TH1D *h_tau2onepass_gen;
      TH1D *h_T2T1onepass_gen; 
      TH2D *h_T2_T1_gen; 
      TH2D *h_tau1_pt_gen;
      TH2D *h_tau2_pt_gen;
      TH2D *h_T2T1_pt_gen;
      TH2D *h_T2T1_pt_noM_gen;*/
      /*TH1D *h_invmass_all;
      TH1D *h_invmass_06;
      TH1D *h_invmass_04;
      TH1D *h_invmass_02;
      TH1D *h_pt_300_449;
      TH1D *h_m_300_449;
      TH1D *h_pt_450_599;
      TH1D *h_m_450_599;*/
      TH1D *h_tau1_reco;
      TH1D *h_tau2_reco;
      TH1D *h_T2T1_reco;
      TH2D *h_tau1_pt_reco;
      TH2D *h_tau2_pt_reco;
      TH2D *h_T2T1_pt_mcut_reco;
      TH2D *h_T2T1_pt_nomcut_reco;
      TH2D *h_T2T1_pv_mcut_reco;
      TH2D *h_T2T1_pv_nomcut_reco;
      TH2D *h_T2T1_pv_mcut_reco2;
      TH2D *h_T2T1_pv_nomcut_reco2;

};

//
// constants, enums and typedefs
//
double ptcut = 300.0;
double etacut = 1.3;
double massmin = 65.0;
double massmax = 95.0;
double fixDeltaR = 0.3;
//
// static data member definitions
//

//
// constructors and destructor
//
RutgersJetAnalyzer::RutgersJetAnalyzer(const edm::ParameterSet& iConfig):
   genJetsTag(iConfig.getParameter<edm::InputTag>("GenJetsTag")),
   genParticleTag(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
   inputsTag(iConfig.getParameter<edm::InputTag>("InputsTag")),
   jetsTag(iConfig.getParameter<edm::InputTag>("JetsTag")),
   pvtag(iConfig.getParameter<edm::InputTag>("PvTag")),
   inputPtMin(iConfig.getParameter<double>("InputPtMin")),
   jetPtMin(iConfig.getParameter<double>("JetPtMin")),
   matching(iConfig.getParameter<double>("Matching")),
   radius(iConfig.getParameter<double>("Radius"))

{
   //now do what ever initialization is needed
   jetDefinitionAK = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, radius) );
   jetDefinitionKT = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, radius) );
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
    //Loop over gen particles to find Ws (pdgId = +-24), Zs (pdgId = 23), Hs (pdgId = 25) and b-quarks (pdgId = +-5)
    edm::Handle<reco::GenParticleCollection> genHandle;
    iEvent.getByLabel(genParticleTag, genHandle);
    //std::cout << "Number of GenParticles: " << genHandle->size() << std::endl;
    //std::cout << "Run #: " << iEvent.id().run() << std::endl;
    //std::cout << "Event #: " << iEvent.id().event() << std::endl;
    //std::cout << "Lumi #: " << iEvent.luminosityBlock() << std::endl;
    std::vector <fastjet::PseudoJet> allWs; 
    for( reco::GenParticleCollection::const_iterator  iter = genHandle->begin(); iter != genHandle->end(); ++iter )
     {
       if( abs(iter->pdgId()) == 24 )
	{
	  fastjet::PseudoJet temp = fastjet::PseudoJet(iter->px(),iter->py(),iter->pz(),iter->energy());
	  allWs.push_back(temp);
	  //std::cout<<"Found W with eta: "<<iter->eta()<<" & with phi: "<<iter->phi()<<std::endl; 
	}
     }
   int numW = allWs.size();
   /*std::vector <fastjet::PseudoJet> allZs;
   for( reco::GenParticleCollection::const_iterator  iter3 = genHandle->begin(); iter3 != genHandle->end(); ++iter3 )
     {
       if( abs(iter3->pdgId()) == 23 )
	{
	  fastjet::PseudoJet temp3 = fastjet::PseudoJet(iter3->px(),iter3->py(),iter3->pz(),iter3->energy());
	  allZs.push_back(temp3);
	}
     }
   int numZ = allZs.size();
   std::vector <fastjet::PseudoJet> allHs;
   for( reco::GenParticleCollection::const_iterator  iter4 = genHandle->begin(); iter4 != genHandle->end(); ++iter4 )
     {
       if( abs(iter4->pdgId()) == 25 )
	{
	  fastjet::PseudoJet temp4 = fastjet::PseudoJet(iter4->px(),iter4->py(),iter4->pz(),iter4->energy());
	  allHs.push_back(temp4);
	}
     }
   int numH = allHs.size();
   std::vector <fastjet::PseudoJet> allbs;
   for( reco::GenParticleCollection::const_iterator iter2 = genHandle->begin(); iter2 != genHandle->end(); ++iter2 )
     {
       if( abs(iter2->pdgId()) == 5 )
	{
	  fastjet::PseudoJet temp2 = fastjet::PseudoJet(iter2->px(),iter2->py(),iter2->pz(),iter2->energy());
	  allbs.push_back(temp2);
	}
     }
   int numb = allbs.size();*/

   //RECONSTRUCTED JET LEVEL (PAT JETS)
   edm::Handle< std::vector<reco::Vertex> > PVs;
   iEvent.getByLabel(pvtag,PVs);
   int numPVs = (int) PVs->size();

   //convert patjets to fastjet PseudoJets
   edm::Handle< std::vector<pat::Jet> > patjets;
   iEvent.getByLabel(jetsTag,patjets);
   double numJets = patjets->size();
   std::vector<double> dR_r(numJets,99); //vector for remembering distance to closest W
   std::vector<int> matchedindex_r(numJets,-1); //matched if index is not -1
   for (int i=0; i < numJets; i++)
     {
	pat::Jet ijet = patjets->at(i);
	std::vector<edm::Ptr<reco::PFCandidate> > recoConstit = patjets->at(i).getPFConstituents();
	double patpx=ijet.px(), patpy=ijet.py(), patpz=ijet.pz(), patE=ijet.energy();
	fastjet::PseudoJet tempjet = fastjet::PseudoJet(patpx,patpy,patpz,patE);
        if (matching==1)
          {
	    for (int k=0; k<numW; k++) //match to W
	      {
	        if (tempjet.delta_R(allWs[k]) <= fixDeltaR && tempjet.delta_R(allWs[k]) <= dR_r[i])
	          {
		    dR_r[i] = tempjet.delta_R(allWs[k]);
		    matchedindex_r[i] = 24; //W jet
		    /*for (int l=0; l<numZ; l++) //look for closer Z
		      {
			if (tempjet.delta_R(allZs[l]) <= fixDeltaR && tempjet.delta_R(allZs[l]) <= dR_r[i])
			  {
			    dR_r[i] = tempjet.delta_R(allZs[l]);
			    matchedindex_r[i] = 23; //Z jet
			    for (int m=0; m<numH; m++) //look for closer H
			      {
				if (tempjet.delta_R(allHs[m]) <= fixDeltaR && tempjet.delta_R(allHs[m]) <= dR_r[i])
				  {
				    dR_r[i] = tempjet.delta_R(allHs[m]);
				    matchedindex_r[i] = 25; //H jet
				  }
			      }
			  }
		      }*/
	          }
	      }
          }
        if (matching==0) matchedindex_r[i] = 7; //no matching applied for QCD

        if (!(tempjet.pt()>ptcut && fabs(tempjet.eta())<etacut && matchedindex_r[i]!=-1)) continue;

	std::vector<fastjet::PseudoJet> fjInputs;
	std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator inBegin = recoConstit.begin(),
        inEnd = recoConstit.end(), m = inBegin;
	for ( ; m != inEnd; ++m )
	{
	  reco::PFCandidatePtr constit = *m;
	  if (constit->pt() < inputPtMin) continue;
	  if (constit->pt() == 0)
	  	{
		  edm::LogError("NullTransverseMomentum") << "dropping input candidate with pt=0";
		  continue;
		}
	  fjInputs.push_back(fastjet::PseudoJet(constit->px(),constit->py(),constit->pz(),constit->energy()));
	  fjInputs.back().set_user_index(m - inBegin);
	}
	
	ClusterSequencePtr reclusterSeq = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *jetDefinitionAK ));
	std::vector<fastjet::PseudoJet> reclusteredJets = fastjet::sorted_by_pt(reclusterSeq->inclusive_jets(jetPtMin));
	std::vector<fastjet::PseudoJet> subJets = fastjet::sorted_by_pt(reclusteredJets[0].exclusive_subjets(2));
	double beta = 1.0;
	double R0 = radius;
	double Rcut = R0;
	fastjet::Nsubjettiness nSub1OnePass_reco(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	double tau1onepass_reco = nSub1OnePass_reco(reclusteredJets[0]);
	fastjet::Nsubjettiness nSub2OnePass_reco(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	double tau2onepass_reco = nSub2OnePass_reco(reclusteredJets[0]);

	std::vector<double> dR_sub1(2,99); //vector for remembering distance to closest b
	std::vector<int> matchedindex_sub1(2,-1); //subjet matched if index is not -1
	/*for (int h=0; h<numb; h++) //match subjets to b-quarks
	 {
		if (subJets1[0].delta_R(allbs[h]) <= fixDeltaR && subJets1[0].delta_R(allbs[h]) <= dR_sub1[0])
		  {
		    dR_sub1[0] = subJets1[0].delta_R(allbs[h]);
		    matchedindex_sub1[0] = 5;
		  }
		if (subJets1[1].delta_R(allbs[h]) <= fixDeltaR && subJets1[1].delta_R(allbs[h]) <= dR_sub1[1])
		  {
		    dR_sub1[1] = subJets1[1].delta_R(allbs[h]);
		    matchedindex_sub1[1] = 5;
		  }
	 }*/

        if (tau1onepass_reco != 0 || tau2onepass_reco !=0)
	  {
	    h_T2T1_pt_nomcut_reco->Fill(reclusteredJets[0].pt(),(tau2onepass_reco/tau1onepass_reco));
	
	    if (tempjet.m()>massmin && tempjet.m()<massmax)
	      {
		h_tau1_reco->Fill(tau1onepass_reco);
		h_tau2_reco->Fill(tau2onepass_reco);
		h_T2T1_reco->Fill(tau2onepass_reco/tau1onepass_reco);
		h_tau1_pt_reco->Fill(reclusteredJets[0].pt(),tau1onepass_reco);
		h_tau2_pt_reco->Fill(reclusteredJets[0].pt(),tau2onepass_reco);
		h_T2T1_pt_mcut_reco->Fill(reclusteredJets[0].pt(),(tau2onepass_reco/tau1onepass_reco));
	      }
	    if (tempjet.pt()>=700 && tempjet.pt()<900 && tempjet.m()>massmin && tempjet.m()<massmax)
	      h_T2T1_pv_mcut_reco->Fill(numPVs,(tau2onepass_reco/tau1onepass_reco));

	    if (tempjet.pt()>=700 && tempjet.pt()<900)
	      h_T2T1_pv_nomcut_reco->Fill(numPVs,(tau2onepass_reco/tau1onepass_reco));

	    if (tempjet.pt()>=300 && tempjet.pt()<500 && tempjet.m()>massmin && tempjet.m()<massmax)
	      h_T2T1_pv_mcut_reco2->Fill(numPVs,(tau2onepass_reco/tau1onepass_reco));

	    if (tempjet.pt()>=300 && tempjet.pt()<500)
	      h_T2T1_pv_nomcut_reco2->Fill(numPVs,(tau2onepass_reco/tau1onepass_reco));
          }
      }
   /*
   //GENERATOR JET LEVEL
   edm::Handle<reco::GenJetCollection> genJets;
   iEvent.getByLabel(genJetsTag, genJets);
   double nJets = genJets->size();
   std::vector<double> dR(nJets,99); //vector for remembering distance to closest W
   std::vector<int> matchedindex(nJets,-1); //matched if index is not -1
   //if ( genJets->size() >= 2 )
   {
     // loop over two leading GenJets
     for ( size_t j = 0; j < genJets->size(); ++j )
     {
       //std::cout << "Jet " << i << " has " << genJets->at(i).getJetConstituents().size() << " constituents" << std::endl;
       std::vector<edm::Ptr<reco::Candidate> > inputs; // inputs for subjet clustering
       inputs = genJets->at(j).getJetConstituents();
       reco::GenJet genjet = genJets->at(j);
       fastjet::PseudoJet *genjetfast = new fastjet::PseudoJet(genjet.px(),genjet.py(),genjet.pz(),genjet.energy());
       //double rap = genjetfast->rap();
       //double phi = genjetfast->phi();
        if (matching==1)
          {
	    for (int k=0; k<numW; k++) //match to W
	      {
	        if (tempjet.delta_R(allWs[k]) <= fixDeltaR && tempjet.delta_R(allWs[k]) <= dR_r[j])
	          {
		    dR_r[j] = tempjet.delta_R(allWs[k]);
		    matchedindex[j] = 24; //W jet
		    for (int l=0; l<numZ; l++) //look for closer Z
		      {
			if (tempjet.delta_R(allZs[l]) <= fixDeltaR && tempjet.delta_R(allZs[l]) <= dR_r[j])
			  {
			    dR_r[j] = tempjet.delta_R(allZs[l]);
			    matchedindex[j] = 23; //Z jet
			    for (int m=0; m<numH; m++) //look for closer H
			      {
				if (tempjet.delta_R(allHs[m]) <= fixDeltaR && tempjet.delta_R(allHs[m]) <= dR_r[j])
				  {
				    dR_r[j] = tempjet.delta_R(allHs[m]);
				    matchedindex[j] = 25; //H jet
				  }
			      }
			  }
		      }
	          }
	      }
          }
       if (matching==0) matchedindex[j] = 7; //no real matching done for QCD, just change index so will still be analyzed
       if (genjetfast->pt() > 30 && (matchedindex[j]!=-1) ) 
	{
          hpt_nocut->Fill(genjetfast->pt());
          heta_nocut->Fill(genjetfast->eta());
          hmass_nocut->Fill(genjetfast->m());
        }
       if (genjetfast->pt()>ptcut && fabs(genjetfast->eta())<etacut && (matchedindex[j]!=-1))
	{ */
	  /*if (genjetfast->pt()>=300 && genjetfast->pt()<450)
	    {
	      h_pt_300_449->Fill(genjetfast->pt());
	      h_m_300_449->Fill(genjetfast->m());
	    }
	  if (genjetfast->pt()>=450 && genjetfast->pt()<600)
	    {
	      h_pt_450_599->Fill(genjetfast->pt());
	      h_m_450_599->Fill(genjetfast->m());
	    }*/
	 /* std::vector<fastjet::PseudoJet> noMfj;
	  std::vector<edm::Ptr<reco::Candidate> >::const_iterator inBegin = inputs.begin(),
            inEnd = inputs.end(), i = inBegin;
	  for ( ; i != inEnd; ++i )
	    {
	      reco::CandidatePtr input3 = *i;
	      if (input3->pt() < inputPtMin) continue;
	      if (input3->pt() == 0)
	      	{
		  edm::LogError("NullTransverseMomentum") << "dropping input candidate with pt=0";
		  continue;
		}
	      noMfj.push_back(fastjet::PseudoJet(input3->px(),input3->py(),input3->pz(),input3->energy()));
	      noMfj.back().set_user_index(i - inBegin);
	    }
	  ClusterSequencePtr clusterSeqKT1 = ClusterSequencePtr( new fastjet::ClusterSequence( noMfj, *jetDefinitionKT ) );
	  std::vector<fastjet::PseudoJet> kTJets1 = fastjet::sorted_by_pt(clusterSeqKT1->inclusive_jets(jetPtMin)); // kT jets
	  double beta = 1.0;
	  double R0 = radius;
	  double Rcut = R0+0.2;
	  fastjet::Nsubjettiness nSub1OnePass_gen1(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	  double tau1onepass_gen1 = nSub1OnePass_gen1(kTJets1[0]);
	  fastjet::Nsubjettiness nSub2OnePass_gen1(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	  double tau2onepass_gen1 = nSub2OnePass_gen1(kTJets1[0]);
	  if (tau1onepass_gen1 != 0 || tau2onepass_gen1 !=0)
	   {
	     h_T2T1_pt_noM_gen->Fill(kTJets1[0].pt(),(tau2onepass_gen1/tau1onepass_gen1));
	   }
	}
       if(genjetfast->pt()>ptcut && fabs(genjetfast->eta())<etacut && massmin<genjetfast->m() && genjetfast->m()<massmax && (matchedindex[j]!=-1) )
	{
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
       std::vector<fastjet::PseudoJet> kTJets = fastjet::sorted_by_pt(clusterSeqKT->inclusive_jets(jetPtMin)); // kT jets
       //Nsubjettiness plugin calculates Tau values:
       double beta = 1.0;
       double R0 = radius;
       double Rcut = R0+0.2;
       fastjet::Nsubjettiness nSub1OnePass_gen(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
       double tau1onepass_gen = nSub1OnePass_gen(kTJets[0]);
       fastjet::Nsubjettiness nSub2OnePass_gen(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
       double tau2onepass_gen = nSub2OnePass_gen(kTJets[0]);
	 if (tau1onepass_gen != 0 || tau2onepass_gen !=0)
	   {
	     h_tau1onepass_gen->Fill(tau1onepass_gen);
	     h_tau2onepass_gen->Fill(tau2onepass_gen);
	     h_T2T1onepass_gen->Fill(tau2onepass_gen/tau1onepass_gen);
	     h_T2_T1_gen->Fill(tau1onepass_gen,tau2onepass_gen);
	     h_tau1_pt_gen->Fill((kTJets[0].pt()),tau1onepass_gen);
	     h_tau2_pt_gen->Fill((kTJets[0].pt()),tau2onepass_gen);
	     h_T2T1_pt_gen->Fill((kTJets[0].pt()),(tau2onepass_gen/tau1onepass_gen));
	     /*h_invmass_all->Fill(kTJets[0].m());
	     if ((tau2onepass_gen/tau1onepass_gen)<0.6) h_invmass_06->Fill(kTJets[0].m());
	     if ((tau2onepass_gen/tau1onepass_gen)<0.4) h_invmass_04->Fill(kTJets[0].m());
	     if ((tau2onepass_gen/tau1onepass_gen)<0.2) h_invmass_02->Fill(kTJets[0].m());*/
	  // }

       //JetDefPtr nsubOnepass_gen_jetDef = JetDefPtr( new fastjet::JetDefinition(new fastjet::NjettinessPlugin(2, Njettiness::onepass_kt_axes, beta, R0, Rcut)) ); //make two subjets from plugin --- DON'T USE THESE THOUGH!

       /*// get 2 subjets only for jets that have at least 2 constituents
       if ( fjInputs.size() >= 2 && kTJets.size() >=1 ) 
	{
         std::vector<fastjet::PseudoJet> subJets = fastjet::sorted_by_pt(kTJets[0].exclusive_subjets(2));
	 long double T1sum=0;
	 long double d0=0;
	 double R1=1.2;
	 long double T2sum=0;
	 for(int k=0; k < (int)fjInputs.size(); k++)
	 {
	   long double dR = genjetfast->delta_R(fjInputs[k]);
	   long double dR1 = subJets[0].delta_R(fjInputs[k]);
	   long double dR2 = subJets[1].delta_R(fjInputs[k]);
	   if(dR1 < dR2)
	    {
	     T2sum += dR1*fjInputs[k].pt();
	    }
	   else
	    {
	     T2sum += dR2*fjInputs[k].pt();
	    }
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
        }*//*
	 } //end boosted & higgs mass conditions
         //std::cout << "dR(subjet1,subjet2) for Jet" << j << ": " << subJets[0].delta_R(subJets[1]) << std::endl;
     } // end for ()
} // end if()*/
   //############################################################################################
   //## This section clusters jets from scratch
   //##
 /*   edm::Handle<reco::CandidateView> inputsHandle;
    iEvent.getByLabel(inputsTag, inputsHandle);
 
    //std::cout << "Number of GenParticles: " << inputsHandle->size() << std::endl;
    std::vector<edm::Ptr<reco::Candidate> > inputs2; // inputs for jet clustering

    for (size_t i = 0; i < inputsHandle->size(); ++i) inputs2.push_back(inputsHandle->ptrAt(i));

    std::vector<fastjet::PseudoJet> fjInputs2; // FastJet inputs for jet clustering
 
    // convert candidates to fastjet::PseudoJets
    std::vector<edm::Ptr<reco::Candidate> >::const_iterator inBegin2 = inputs2.begin(),
      inEnd = inputs2.end(), i = inBegin2;
    for (; i != inEnd; ++i )
    {
      reco::CandidatePtr input2 = *i;
      if (input2->pt() < inputPtMin) continue;
      if (input2->pt() == 0) {
        edm::LogError("NullTransverseMomentum") << "dropping input candidate with pt=0";
        continue;
      }

      fjInputs2.push_back(fastjet::PseudoJet(input2->px(),input2->py(),input2->pz(), input2->energy()));
      fjInputs2.back().set_user_index(i - inBegin2);
    }
 
   // define clustering sequence & get clustered jets
   ClusterSequencePtr clusterSeqAK = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs2, *jetDefinitionAK ) );
   std::vector<fastjet::PseudoJet> jetsAK = fastjet::sorted_by_pt(clusterSeqAK->inclusive_jets(jetPtMin));
   //std::cout << "Number of GenJets: " << jetsAK.size() << std::endl;
   if ( jetsAK.size() >= 2 )
   {
     // loop over two leading GenJets
     for ( size_t j = 0; j < 2; ++j )
     { 
       fastjet::PseudoJet *genjet_cluster = new fastjet::PseudoJet(jetsAK[j].px(),jetsAK[j].py(),jetsAK[j].pz(),jetsAK[j].E());
       if (genjet_cluster->pt() > 30)
	{
	  hpt_nocut_cluster->Fill(genjet_cluster->pt());
	  heta_nocut_cluster->Fill(genjet_cluster->eta());
	  hmass_nocut_cluster->Fill(genjet_cluster->m());
	}
     }
   }*/
   //############################################################################################
   //## This section demonstrates how to initialize N-subjettiness
   //##
   /*double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
   double R0 = 1.2; // Characteristic jet radius for normalization
   double Rcut = 1.4; // maximum R particles can be from axis to be included in jet
   if ( jetsAK.size() >= 2){
     for (size_t k=0; k < 2; ++k ){
	if (jetsAK[k].pt() > 200 && fabs(jetsAK[k].eta()) < 2.5 && jetsAK[k].m() > 95 && jetsAK[k].m() < 125) {
	  fastjet::Nsubjettiness nSub1OnePass(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	  double tau1onepass = nSub1OnePass(jetsAK[k]);
	  fastjet::Nsubjettiness nSub2OnePass(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
	  double tau2onepass = nSub2OnePass(jetsAK[k]);
	    if (tau1onepass != 0 || tau2onepass !=0){
	      h_tau1onepass_cluster->Fill(tau1onepass);
	      h_tau2onepass_cluster->Fill(tau2onepass);
	      h_T2T1onepass_cluster->Fill(tau2onepass/tau1onepass);
	    }
	  //JetDefPtr nsubOnepass_jetDef = JetDefPtr( new fastjet::JetDefinition(new fastjet::NjettinessPlugin(2, Njettiness::onepass_kt_axes, beta, R0, Rcut)) );
	  
	}//end if condition on AK jets
     }//end for loop over 2 leading jets
   }//end condition must have 2 jets*/
   //##
   //##
   //############################################################################################
}


// ------------ method called once each job just before starting event loop  ------------
void 
RutgersJetAnalyzer::beginJob()
{
  //hT1 = fs->make<TH1D>("hT1" , "T_1" , 100 , 0 , 1 );
  //hT2 = fs->make<TH1D>("hT2" , "T_2" , 100 , 0 , 1 );
  //hT2byT1 = fs->make<TH1D>("hT2/T1" , "T_2/T_1" , 100 , 0 , 1 );
  //hpt_nocut = fs->make<TH1D>("hpt_nocut","Jet pT with no cuts",250,0,450);
  //hpt_cut = fs->make<TH1D>("hpt_cut","Jet pT with inv mass cut",250,0,450);
  //heta_nocut = fs->make<TH1D>("heta_nocut","Jet eta with no cuts",450,-5,5);
  //heta_cut = fs->make<TH1D>("heta_cut","Jet eta with cuts",200,-5,5);
  //hmass_nocut = fs->make<TH1D>("hmass_nocut","Jet mass with no cuts",300,0,300);
  //hmass_cut = fs->make<TH1D>("hmass_cut","Jet mass with cuts",300,0,300);
  //hpt_nocut_cluster = fs->make<TH1D>("hpt_nocut_cluster","Jet pT with no cuts",250,0,450);
  //heta_nocut_cluster = fs->make<TH1D>("heta_nocut_cluster","Jet eta with no cuts",450,-5,5);
  //hmass_nocut_cluster = fs->make<TH1D>("hmass_nocut_cluster","Jet mass with no cuts",300,0,300);
  /*h_tau1onepass_gen = fs->make<TH1D>("h_tau1onepass_gen","Tau1 OnePass kT CMSSW Gen Jet clustering",100,0,1);  
  h_tau2onepass_gen = fs->make<TH1D>("h_tau2onepass_gen","Tau2 OnePass kT CMSSW Gen Jet clustering",100,0,1);
  h_T2T1onepass_gen = fs->make<TH1D>("h_T2T1onepass_gen","T2/T1 OnePass kT CMSSW Gen Jet clustering",100,0,1);
  h_T2_T1_gen = fs->make<TH2D>("h_T2_T1_gen","Tau2 vs. Tau1 OnePass kT CMSSW Gen Jet clustering",100,0,1,100,0,1);
  h_tau1_pt_gen = fs->make<TH2D>("h_tau1_pt_gen","Tau1 OnePass kT CMSSW Gen Jet clustering vs kT Jet pT",500,0,500,100,0,1);
  h_tau2_pt_gen = fs->make<TH2D>("h_tau2_pt_gen","Tau2 OnePass kT CMSSW Gen Jet clustering vs kT Jet pT",500,0,500,100,0,1);
  h_T2T1_pt_gen = fs->make<TH2D>("h_T2T1_pt_gen","T2/T1 OnePass kT Gen Jet vs kT Jet pT with mass cut",1000,0,1000,100,0,1);
  h_T2T1_pt_noM_gen = fs->make<TH2D>("h_T2T1_noM_gen","T2/T1 OnePass kT Gen Jet vs kT Jet pT no mass cut",1000,0,1000,100,0,1);*/
  //*****inorder for invmass plots to be usefull*** need to remove invmass window requirements!********
  /*h_invmass_all = fs->make<TH1D>("h_invmass_all","Jet Inv Mass for all Tau2/Tau1 ",80,40,120);
  h_invmass_06 = fs->make<TH1D>("h_invmass_06","Jet Inv Mass for Tau2/Tau1 < 0.6",80,40,120);
  h_invmass_04 = fs->make<TH1D>("h_invmass_04","Jet Inv Mass for Tau2/Tau1 < 0.4",80,40,120);
  h_invmass_02 = fs->make<TH1D>("h_invmass_02","Jet Inv Mass for Tau2/Tau1 < 0.2",80,40,120);
  h_pt_300_449 = fs->make<TH1D>("h_pt_300_449","Jet pT for pT 300 - 449",150,300,450);
  h_m_300_449 = fs->make<TH1D>("h_m_300_449","Jet Inv Mass for pT 300 - 449",150,300,450);
  h_pt_450_599 = fs->make<TH1D>("h_pt_450_599","Jet pT for pT 450 - 600",150,450,600);
  h_m_450_599 = fs->make<TH1D>("h_m_450_599","Jet Inv Mass for pT 450 - 600",150,450,600);*/
  h_tau1_reco = fs->make<TH1D>("h_tau1_reco","Tau1 OnePass kT CMSSW reco",100,0,1);  
  h_tau2_reco = fs->make<TH1D>("h_tau2_reco","Tau2 OnePass kT CMSSW reco",100,0,1);  
  h_T2T1_reco = fs->make<TH1D>("h_T2T1_reco","T2/T1 OnePass kT CMSSW reco",100,0,1);  
  h_tau1_pt_reco = fs->make<TH2D>("h_tau1_pt_reco","Tau1 OnePass kT CMSSW reco vs kT Jet pT",500,0,500,100,0,1);
  h_tau2_pt_reco = fs->make<TH2D>("h_tau2_pt_reco","Tau2 OnePass kT CMSSW reco vs kT Jet pT",500,0,500,100,0,1);
  h_T2T1_pt_mcut_reco = fs->make<TH2D>("h_T2T1_pt_mcut_reco","T2/T1 vs. pT with mcut",1000,0,1000,100,0,1);
  h_T2T1_pt_nomcut_reco = fs->make<TH2D>("h_T2T1_pt_nomcut_reco","T2/T1 vs. pT with nomcut",1000,0,1000,100,0,1);
  h_T2T1_pv_mcut_reco = fs->make<TH2D>("h_T2T1_pv_mcut_reco","T2/T1 vs. Num of PVs mcut for pT 700-900 GeV",50,0,50,100,0,1);
  h_T2T1_pv_nomcut_reco = fs->make<TH2D>("h_T2T1_pv_nomcut_reco","T2/T1 vs. Num of PVs nomcut for pT 700-900 GeV",50,0,50,100,0,1);
  h_T2T1_pv_mcut_reco2 = fs->make<TH2D>("h_T2T1_pv_mcut_reco2","T2/T1 vs. Num of PVs mcut for pT 300-500 GeV",50,0,50,100,0,1);
  h_T2T1_pv_nomcut_reco2 = fs->make<TH2D>("h_T2T1_pv_nomcut_reco2","T2/T1 vs. Num of PVs nomcut for pT 300-500 GeV",50,0,50,100,0,1);
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
RutgersJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RutgersJetAnalyzer);
