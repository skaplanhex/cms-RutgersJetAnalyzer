###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('outFilename',
                 "outfile.root",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "output file name"
                 )
options.register('reportEvery',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Report every N events (default is N=10)"
                 )
options.register('wantSummary',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print out trigger and timing summary"
                 )
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', -1)

options.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
DUMMY_INPUTFILES
    )
)

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(DUMMY_OUTPUTFILE)
)

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
    GenParticleTag       = cms.InputTag('prunedGenParticles'),
    JetsTag              = cms.InputTag('selectedPatJetsAK6PF'),
    GroomedJetsTag       = cms.InputTag('selectedPatJetsAK6TrimmedPF'),
    PvTag                = cms.InputTag('goodOfflinePrimaryVertices'),
    FJInputPtMin         = cms.double(0.),
    FJJetPtMin           = cms.double(0.),
    FJRadius             = cms.double(0.6),
    DoWMatching          = cms.bool(True),
    WMatchingRadius      = cms.double(0.3),
    LeptonMatchingRadius = cms.double(0.4),
    JetPtMin             = cms.double(300.),
    JetAbsEtaMax         = cms.double(1.5),
    JetMassMin           = cms.double(65.),
    JetMassMax           = cms.double(95.),
    UseGroomedJets       = cms.bool(True),
    UseEventWeight       = cms.bool(False)
)

## Path definition
process.p = cms.Path(
     process.goodOfflinePrimaryVertices*
     process.rutgersJetAnalyzer
)

## Schedule definition
process.schedule = cms.Schedule(process.p)
