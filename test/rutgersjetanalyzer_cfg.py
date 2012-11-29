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
options.setDefault('maxEvents', -1)

options.parseArguments()

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
	DUMMY_INPUTFILES	

#'file:/uscmst1b_scratch/lpc1/lpcphys/ferencek/Higgs/CMSSW_5_2_5/src/RutgersSandbox/RutgersJetAnalyzer/test/patTuple_PF2PAT.root'


#flat QCD sample -------
       #'/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/'
#WW events -------------
	#'/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/
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
    )
)

# Load module that produces GenJetParticles
process.load('RecoJets.Configuration.GenJetParticles_cff')

# Define Anti-kT jets with R=0.6
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak6GenJets = ak5GenJets.clone( rParam = 0.6, src = cms.InputTag("genParticlesForJetsNoNu") )

## Produce a collection of good primary vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone(
        minNdof = cms.double(4.0), # this is >= 4
        maxZ = cms.double(24.0),
        maxRho = cms.double(2.0)
    ),
    src = cms.InputTag('offlinePrimaryVertices')
)

# Initialize RutgersJetAnalyzer
process.rutgersJetAnalyzer = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenJetsTag = cms.InputTag('ak6GenJets'),
    GenParticleTag = cms.InputTag('prunedGenParticles'),
    InputsTag = cms.InputTag('genParticlesForJetsNoNu'),
    JetsTag = cms.InputTag('selectedPatJetsAK6PF'),
    GroomedJetsTag = cms.InputTag('selectedPatJetsAK6PrunedPF'),
    PvTag = cms.InputTag('goodOfflinePrimaryVertices'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0),
    Matching = cms.double(1),
    DoWMatching = cms.bool(True),
    WMatchingRadius = cms.double(0.3),
    LeptonMatchingRadius = cms.double(0.4),
    Radius = cms.double(0.6),
    UseGroomedJets = cms.bool(True)
)

#if Matching==0 NO MATCHING DONE - use for QCD	
#if Matching==1 MATCHING IS DONE - use for WW

#specify radius for jet clustering and rParam in defining Anti-kT jets

# Path definition
#process.p = cms.Path(process.genParticlesForJetsNoNu*process.ak6GenJets*process.rutgersJetAnalyzer)
process.p = cms.Path(process.goodOfflinePrimaryVertices*process.rutgersJetAnalyzer)
#Output histograms to root file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(DUMMY_OUTPUTFILE)
	#options.outputfile)
    )
