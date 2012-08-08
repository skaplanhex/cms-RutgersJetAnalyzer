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
    'START52_V9::All',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global tag to be used"
)
options.register('reportEvery',
    10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=10)"
)
options.register('usePFchs',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PFchs"
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

## Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = options.globalTag

## Events to process
process.maxEvents.input = options.maxEvents

## Options and Output Report
process.options.wantSummary = False

## Input files
process.source.fileNames = [
    '/store/mc/Summer12/ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigpp/AODSIM/PU_S7_START52_V9-v1/0000/FA40039F-03B4-E111-A081-0030486790BE.root'
]

## Configure PAT to use PF2PAT instead of AOD sources
## this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
jetAlgo="AK7"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
          jetCorrections=inputJetCorrLabelAK7, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

## Top projections in PF2PAT
getattr(process,"pfNoPileUp"+postfix).enable = options.usePFchs
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

## Load GenJetParticles
process.load('RecoJets.Configuration.GenJetParticles_cff')

## Define AK6 jets (GEN and RECO)
from RecoJets.Configuration.RecoGenJets_cff import ak5GenJetsNoNu
process.ak6GenJetsNoNu = ak5GenJetsNoNu.clone( rParam = 0.6 )
process.ak6PFlow = process.pfJetsPFlow.clone( rParam = 0.6 )
## Add AK6 jets to PAT
addJetCollection(
    process,
    cms.InputTag('ak6PFlow'),
    'AK6', 'PF',
    doJTA=True,
    doBTagging=True,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID = False,
    genJetCollection = cms.InputTag("ak6GenJetsNoNu")
)

## Define a sequence for RECO jets and append it to the PF2PAT sequence
process.jetSeq = cms.Sequence(
    process.ak6PFlow
)
process.PFBRECOPFlow *= process.jetSeq

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

## Path definition
process.p = cms.Path(
    (
    process.genParticlesForJetsNoNu*
    process.ak6GenJetsNoNu+
    process.goodOfflinePrimaryVertices
    )*
    getattr(process,"patPF2PATSequence"+postfix)
)

### Add PF2PAT output to the created file
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
##process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
##process.out.outputCommands =  cms.untracked.vstring('drop *')
#process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   #'keep recoPFCandidates_particleFlow_*_*',
                                                   #*patEventContentNoCleaning )

#process.out.fileName = 'patTuple_PF2PAT.root'

## Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

## Schedule definition
process.schedule = cms.Schedule(process.p)
#process.schedule.append(process.outpath)
