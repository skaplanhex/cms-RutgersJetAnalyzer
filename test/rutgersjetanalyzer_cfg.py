import FWCore.ParameterSet.Config as cms

###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('reportEvery',
    1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1)"
)

options.register('outputfile',
		 "outputfile.root",
		 VarParsing.multiplicity.singleton,
		 VarParsing.varType.string,
		 "outputfile"
		 )
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 30000)

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

# Options and output report
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

#flat QCD sample -------
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/FE6F3A25-D998-E111-B776-0030487D5E45.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/FCE8334B-C098-E111-A8D4-0025904B12FE.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/FA9CC11F-D498-E111-BAC9-0025901D40CA.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/F8A0ECCC-E098-E111-9407-0030487F1A49.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/F88200CA-DE98-E111-B250-0030487D858D.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/F6901201-DA98-E111-B667-003048D43958.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/F49003A6-BE98-E111-B664-00266CF25E44.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/F2FC706C-DF98-E111-A859-0030487D5EBB.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/F246CC96-D898-E111-8E8C-003048F0E18C.root',
       '/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v5/0000/F087AC48-DA98-E111-95EE-003048C692C0.root'

#WW events -------------
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/36B20258-8A95-E111-A646-0026189438BC.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/38E44AE5-8595-E111-A05F-0030486790BE.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/36EDDBA8-D495-E111-8FF1-0026189438D4.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0AD250E5-9C95-E111-93DA-003048D15E02.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0AD5E4E0-8295-E111-BB9C-002618FDA204.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0AE0567D-9195-E111-9881-002354EF3BDA.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0C0C4F0A-7A95-E111-80ED-001A92810ADE.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0C49E64A-9C95-E111-AA4F-0026189438ED.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0CAF6BD5-9495-E111-8F73-001A928116AE.root',
 	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0CC9D92B-8995-E111-83DD-003048678F8A.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0CF3A916-7C95-E111-9FC7-003048678FF8.root',
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/0E0C9FE3-8595-E111-81B7-003048678FB4.root',

#HH events -------------
	#'file:/eos/uscms/store/user/skaplan/dihiggs_1-9999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_10000-19999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_20000-29999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_30000-39999.root',
        #'file:/eos/uscms/store/user/skaplan/dihiggs_40000-49999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_50000-59999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_60000-69999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_70000-79999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_80000-89999.root',
	#'file:/eos/uscms/store/user/skaplan/dihiggs_90000-100000.root'
#Example file ----------
       # '/store/mc/Summer12/ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigpp/AODSIM/PU_S7_START52_V9-v1/0000/FA40039F-03B4-E111-A081-0030486790BE.root'
    )
)

# Load module that produces GenJetParticles
process.load('RecoJets.Configuration.GenJetParticles_cff')

# Define Anti-kT jets with R=1.2
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak6GenJets = ak5GenJets.clone( rParam = 0.6, src = cms.InputTag("genParticlesForJetsNoNu") )

# Initialize RutgersJetAnalyzer
process.rutgersJetAnalyzer = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenJetsTag = cms.InputTag('ak6GenJets'),
    InputsTag = cms.InputTag('genParticlesForJetsNoNu'),
    GenParticleTag = cms.InputTag('genParticles'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0)
)

# Path definition
process.p = cms.Path(process.genParticlesForJetsNoNu*process.ak6GenJets*process.rutgersJetAnalyzer)
#Output histograms to root file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputfile)
    )
