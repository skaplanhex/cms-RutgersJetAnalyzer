###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('outfilename',
                 "outfile.root",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "output file name"
                 )
options.register('jetCollection',
                 "selectedPatJetsAK6PF",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Jet collection used in analysis"
                 )
options.register('reportEvery',
                 10,
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
options.setDefault('maxEvents', 10)

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
      'file:///cms/ferencek/store/skaplan/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1_PATTuple/f7377c82d9b827962c4fc87795e9adaf/patTuple_PF2PAT_100_2_6Ow.root'
      )
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

process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outfilename)
)

process.rutgersJetAnalyzerAK6 = cms.EDAnalyzer('RutgersJetAnalyzer',
   GenParticleTag = cms.InputTag('prunedGenParticles'),
   JetsTag = cms.InputTag('selectedPatJetsAK6PF'),
   pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
   InputPtMin = cms.double(0.0),
   JetPtMin = cms.double(0.0)
   )
process.rutgersJetAnalyzerAK7 = cms.EDAnalyzer('RutgersJetAnalyzer',
   GenParticleTag = cms.InputTag('prunedGenParticles'),
   JetsTag = cms.InputTag('selectedPatJetsPFlow'),
   pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
   InputPtMin = cms.double(0.0),
   JetPtMin = cms.double(0.0)
   )
process.rutgersJetAnalyzerAK8 = cms.EDAnalyzer('RutgersJetAnalyzer',
   GenParticleTag = cms.InputTag('prunedGenParticles'),
   JetsTag = cms.InputTag('selectedPatJetsAK8PF'),
   pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
   InputPtMin = cms.double(0.0),
   JetPtMin = cms.double(0.0)
   )
process.rutgersJetAnalyzerAK10 = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenParticleTag = cms.InputTag('prunedGenParticles'),
    JetsTag = cms.InputTag('selectedPatJetsAK10PF'),
    pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0)
    )
## Path definition
process.p = cms.Path(
     process.goodOfflinePrimaryVertices*
     process.rutgersJetAnalyzerAK6*
     process.rutgersJetAnalyzerAK7*
     process.rutgersJetAnalyzerAK8*
     process.rutgersJetAnalyzerAK10
     )

### Add PF2PAT output to the created file
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
##process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
##process.out.outputCommands =  cms.untracked.vstring('drop *')
#process.out.outputCommands = cms.untracked.vstring('drop *',
#'keep recoPFCandidates_particleFlow_*_*',
#*patEventContentNoCleaning )

## Delete predefined Endpath (needed for running with CRAB)
#del process.out
#del process.outpath

## Schedule definition
process.schedule = cms.Schedule(process.p)
#process.schedule.append(process.outpath)