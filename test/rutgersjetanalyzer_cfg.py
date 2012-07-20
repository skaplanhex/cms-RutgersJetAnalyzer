import FWCore.ParameterSet.Config as cms

###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('reportEvery',
    1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1)"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 10)

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
        '/store/mc/Summer12/ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigpp/AODSIM/PU_S7_START52_V9-v1/0000/FA40039F-03B4-E111-A081-0030486790BE.root'
    )
)

# Load module that produces GenJetParticles
process.load('RecoJets.Configuration.GenJetParticles_cff')

# Define Anti-kT jets with R=1.2
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak12GenJets = ak5GenJets.clone( rParam = 1.2 )

# Initialize RutgersJetAnalyzer
process.rutgersJetAnalyzer = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenJetsTag = cms.InputTag('ak12GenJets'),
    InputsTag = cms.InputTag('genParticlesForJets'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0)
)

# Path definition
process.p = cms.Path(process.genJetParticles*process.ak12GenJets*process.rutgersJetAnalyzer)
