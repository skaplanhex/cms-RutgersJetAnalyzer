###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('runOnData',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag',
    'START53_V27',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global tag to be used"
)
options.register('outFilename',
    'patTuple_PF2PAT_v3.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('useExternalInput', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use external input"
)
options.register('externalInput', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Path of an external list of input files"
)
options.register('dumpPythonCfg', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of the rewritten cfg file"
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
options.register('usePFchs',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PFchs"
)
options.register('useExtPFchs',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use extended PFchs"
)
options.register('useHLTFiltering',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use HLT filtering"
)
options.register('triggers',
    'HLT_PFJet400_*',
    VarParsing.multiplicity.list,
    VarParsing.varType.string,
    "List of triggers for HLT filtering"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 10)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Jet energy corrections
inputJetCorrLabelAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
inputJetCorrLabelAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if not options.usePFchs:
    inputJetCorrLabelAK5 = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    inputJetCorrLabelAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if options.runOnData:
    inputJetCorrLabelAK5[1].append('L2L3Residual')
    inputJetCorrLabelAK7[1].append('L2L3Residual')

import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = options.globalTag + '::All'

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12_DR53X/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FA0D87EB-6CE9-E111-BA2C-0018F3D096F0.root'
    )
)
# If using external input files
if options.useExternalInput:
    process.source.fileNames = cms.untracked.vstring( open(options.externalInput,"r").read().splitlines() )

if options.runOnData:
    process.source.fileNames = ['/store/data/Run2012A/Jet/AOD/13Jul2012-v1/00000/449F09C1-94D9-E111-BEED-00266CFACC38.root']

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Configure PAT to use PF2PAT instead of AOD sources
## this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outFilename),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)

postfix = "PFlow"
jetAlgo="AK7"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
          jetCorrections=inputJetCorrLabelAK7, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

## Top projections in PF2PAT
getattr(process,"pfPileUp"+postfix).checkClosestZVertex = False
if options.useExtPFchs:
    getattr(process,"pfPileUp"+postfix).useJets = True
    getattr(process,"pfPileUp"+postfix).useMuons = True
getattr(process,"pfNoPileUp"+postfix).enable = options.usePFchs
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done to turn off b tagging)
switchJetCollection(process,
    cms.InputTag("pfNoTauPFlow"),
    doJTA        = False,
    doBTagging   = False,
    jetCorrLabel = inputJetCorrLabelAK7,
    doType1MET   = False,
    genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
    doJetID      = False,
    postfix = postfix
)

from PhysicsTools.PatAlgos.tools.coreTools import *
## Remove taus from the PAT sequence
removeSpecificPATObjects(process,names=['Taus'],postfix=postfix)

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

## Filter based on HLT paths
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths = cms.vstring(options.triggers),        # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''),             # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),                          # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True)                           # throw exception on unknown path names
)

## Good primary vertex event filter
process.primaryVertexFilter = cms.EDFilter('VertexSelector',
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string('!isFake & ndof > 4 & abs(z) <= 24 & position.Rho <= 2'),
    filter = cms.bool(True)
)

## Define a sequence of trigger filters
process.trigSeq = cms.Sequence(
    process.hltFilter
)

## We only want to use HLT filtering on data
if not options.useHLTFiltering:
    process.trigSeq.remove(process.hltFilter)

## Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
patEventContentNoCleaning.append('drop *_selectedPatPFParticlesPFlow_*_*')
patEventContentNoCleaning.append('drop *_selectedPatJetsPFlow_*_*')
patEventContentNoCleaning.append('keep GenEventInfoProduct_generator_*_*')
patEventContentNoCleaning.append('keep PileupSummaryInfos_*_*_*')
## GenParticles
patEventContentNoCleaning.append('keep recoGenParticles_genParticles_*_*')
## For PF jets
patEventContentNoCleaning.append('keep *_pfNoElectronPFlow_*_*')
patEventContentNoCleaning.append('keep *_kt6PFJets_rho_*')
## HLT trigger results
patEventContentNoCleaning.append('drop *_TriggerResults_*_*')
patEventContentNoCleaning.append('keep *_TriggerResults_*_HLT')
## For b tagging
patEventContentNoCleaning.append('keep *_offlineBeamSpot_*_*')
patEventContentNoCleaning.append('keep *_generalTracks_*_*')
patEventContentNoCleaning.append('keep *_goodOfflinePrimaryVertices_*_*')
patEventContentNoCleaning.append('keep recoMuons_muons_*_*')
patEventContentNoCleaning.append('keep recoTracks_globalMuons_*_*')
patEventContentNoCleaning.append('keep recoTracks_standAloneMuons_*_*')
patEventContentNoCleaning.append('keep *_gsfElectrons_*_*')
patEventContentNoCleaning.append('keep *_gsfElectronCores_*_*')
patEventContentNoCleaning.append('keep *_electronGsfTracks_*_*')

process.out.outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning)

## Path definition
process.p = cms.Path(
    process.primaryVertexFilter
    * process.trigSeq
    * process.goodOfflinePrimaryVertices
    * getattr(process,"patPF2PATSequence"+postfix)
)

## Define Endpath
process.outpath = cms.EndPath(process.out)

## Schedule definition
process.schedule = cms.Schedule(process.p)
process.schedule.append(process.outpath)

## Rewrite the cfg file
if options.dumpPythonCfg != '':
    open(options.dumpPythonCfg,'w').write(process.dumpPython())
