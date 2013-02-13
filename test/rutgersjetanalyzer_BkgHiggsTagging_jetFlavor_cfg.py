import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import copy

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('runOnData',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag',
    'START53_V7F',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global tag to be used"
)
options.register('outFilename',
    'outfile.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('reportEvery',
    100,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=100)"
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
options.register('jetRadius',
    0.8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Distance parameter R for jet clustering (default is 0.8)"
)
options.register('doJTA',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run jet-track association"
)
options.register('useExplicitJTA',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use explicit jet-track association"
)
options.register('doBTagging',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run b tagging"
)
options.register('runOnSignal',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on signal"
)
options.register('objectType',
    'H',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Heavy object type (H or W)"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', -1)

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

## b tagging
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos','inclusiveSecondaryVertexFinderTagInfos','inclusiveSecondaryVertexFinderTagInfosFiltered']
             #,'softMuonTagInfos','softElectronTagInfos']
bTagDiscriminators = ['jetProbabilityBJetTags','trackCountingHighPurBJetTags', 'trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags',
                      'simpleSecondaryVertexHighPurBJetTags','simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags',
                      'combinedSecondaryVertexBJetTags','combinedInclusiveSecondaryVertexBJetTags','doubleSecondaryVertexHighEffBJetTags']
                      #'combinedSecondaryVertexMVABJetTags','softMuonBJetTags','softMuonByPtBJetTags','softMuonByIP3dBJetTags','combinedMVABJetTags']
bTagDiscriminatorsSub = copy.deepcopy(bTagDiscriminators)
bTagDiscriminatorsSub.remove('doubleSecondaryVertexHighEffBJetTags')


process = cms.Process("USER")

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
DUMMY_INPUTFILES
    )
)

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(DUMMY_OUTPUTFILE)
)

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

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
from PhysicsTools.PatAlgos.tools.coreTools import *
## Only keep jets in the default sequence
removeAllPATObjectsBut(process, ['Jets'])

## Remove MC matching when running over data
if options.runOnData:
    removeMCMatching( process, ['All'] )

## GenParticles for GenJets
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu = genParticlesForJetsNoNu

## Define Anti-kT jets (GEN and RECO)
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.akGenJetsNoNu = ak5GenJets.clone(
    rParam = options.jetRadius,
    src = cms.InputTag("genParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
process.akPFJets = ak5PFJets.clone(
    rParam = options.jetRadius,
    doAreaFastjet = cms.bool(True),
    src = cms.InputTag("pfNoElectronPFlow"),
    srcPVs = cms.InputTag("goodOfflinePrimaryVertices")
)

## Anti-kT jets and subjets (GEN and RECO) (two collections are produced)
from RutgersSandbox.RutgersSubJetAlgorithm.ak5GenJetsRU_cfi import ak5GenJetsRU
process.akGenJetsNoNuRU = ak5GenJetsRU.clone(
    rParam = process.akGenJetsNoNu.rParam,
    src = process.akGenJetsNoNu.src,
    srcPVs = process.akGenJetsNoNu.srcPVs
)
from RutgersSandbox.RutgersSubJetAlgorithm.ak5PFJetsRU_cfi import ak5PFJetsRU
process.akPFJetsRU = ak5PFJetsRU.clone(
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet
)

## Anti-kT groomed jets
from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
process.akPFJetsTrimmed = ak5PFJetsTrimmed.clone(
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.akPFJetsFiltered = ak5PFJetsFiltered.clone(
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet,
    writeCompound = cms.bool(False),
    jetCollInstanceName=cms.string("")
)
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.akGenJetsNoNuFiltered = ak5GenJets.clone(
    rParam = options.jetRadius,
    src = cms.InputTag("genParticlesForJetsNoNu"),
    useFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.akPFJetsFilteredCompound = ak5PFJetsFiltered.clone(
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.akPFJetsPruned = ak5PFJetsPruned.clone(
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet,
    writeCompound = cms.bool(False),
    jetCollInstanceName=cms.string("")
)
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.akGenJetsNoNuPruned = ak5GenJets.clone(
    SubJetParameters,
    rParam = options.jetRadius,
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.akPFJetsPrunedCompound = ak5PFJetsPruned.clone(
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

## CA jets and subjets (GEN and RECO) (two collections are produced)
from RutgersSandbox.RutgersSubJetAlgorithm.ak5GenJetsRU_cfi import ak5GenJetsRU
process.caGenJetsNoNuRU = ak5GenJetsRU.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = process.akGenJetsNoNu.rParam,
    src = process.akGenJetsNoNu.src,
    srcPVs = process.akGenJetsNoNu.srcPVs
)
from RutgersSandbox.RutgersSubJetAlgorithm.ak5PFJetsRU_cfi import ak5PFJetsRU
process.caPFJetsRU = ak5PFJetsRU.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet
)

## Define CA jets (GEN and RECO)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.caGenJetsNoNu = ca4GenJets.clone(
    rParam = options.jetRadius,
    src = cms.InputTag("genParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.caPFJets = ca4PFJets.clone(
    rParam = options.jetRadius,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.caPFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = options.jetRadius,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet,
    writeCompound = cms.bool(False),
    jetCollInstanceName=cms.string("")
)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.caGenJetsNoNuFiltered = ca4GenJets.clone(
    rParam = options.jetRadius,
    src = cms.InputTag("genParticlesForJetsNoNu"),
    useFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.caPFJetsFilteredCompound = ak5PFJetsFiltered.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.caGenJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    rParam = options.jetRadius,
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.caPFJetsPrunedCompound = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = process.akPFJets.rParam,
    src = process.akPFJets.src,
    srcPVs = process.akPFJets.srcPVs,
    doAreaFastjet = process.akPFJets.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

## PATify above jets
from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection
switchJetCollection(process,
    cms.InputTag('akPFJets'),
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    genJetCollection=cms.InputTag("akGenJetsNoNu"),
    doJetID=False
)
#addJetCollection(
    #process,
    #cms.InputTag('akPFJetsTrimmed'),
    #'AKTrimmed','PF',
    #doJTA=False,
    #doBTagging=False,
    #btagInfo=bTagInfos,
    #btagdiscriminators=bTagDiscriminators,
    #jetCorrLabel=inputJetCorrLabelAK7,
    #doType1MET=False,
    #doL1Cleaning=False,
    #doL1Counters=False,
    #doJetID=False,
    #genJetCollection=cms.InputTag("akGenJetsNoNu")
#)
#addJetCollection(
    #process,
    #cms.InputTag('akPFJetsFiltered'),
    #'AKFiltered','PF',
    #doJTA=False,
    #doBTagging=False,
    #btagInfo=bTagInfos,
    #btagdiscriminators = bTagDiscriminators,
    #jetCorrLabel=inputJetCorrLabelAK7,
    #doType1MET=False,
    #doL1Cleaning=False,
    #doL1Counters=False,
    #doJetID=False,
    #genJetCollection=cms.InputTag("akGenJetsNoNu")
#)
#addJetCollection(
    #process,
    #cms.InputTag('akPFJetsPruned'),
    #'AKPruned','PF',
    #doJTA=False,
    #doBTagging=False,
    #btagInfo=bTagInfos,
    #btagdiscriminators=bTagDiscriminators,
    #jetCorrLabel=inputJetCorrLabelAK7,
    #doType1MET=False,
    #doL1Cleaning=False,
    #doL1Counters=False,
    #doJetID=False,
    #genJetCollection=cms.InputTag("akGenJetsNoNu")
#)
addJetCollection(
    process,
    cms.InputTag('caPFJets'),
    'CA','PF',
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("caGenJetsNoNu")
)
addJetCollection(
    process,
    cms.InputTag('caPFJetsPruned'),
    'CAPruned','PF',
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("caGenJetsNoNu")
)

## If running H tagging
if options.objectType=='H':
    #addJetCollection(
        #process,
        #cms.InputTag('akPFJets'),
        #'AKJTA', 'PF',
        #doJTA=options.doJTA,
        #doBTagging=options.doBTagging,
        #btagInfo=bTagInfos,
        #btagdiscriminators=bTagDiscriminators,
        #jetCorrLabel=inputJetCorrLabelAK7,
        #doType1MET=False,
        #doL1Cleaning=False,
        #doL1Counters=False,
        #doJetID=False,
        #genJetCollection=cms.InputTag('akGenJetsNoNu')
    #)
    #addJetCollection(
        #process,
        #cms.InputTag('caPFJets'),
        #'CAJTA', 'PF',
        #doJTA=options.doJTA,
        #doBTagging=options.doBTagging,
        #btagInfo=bTagInfos,
        #btagdiscriminators=bTagDiscriminators,
        #jetCorrLabel=inputJetCorrLabelAK7,
        #doType1MET=False,
        #doL1Cleaning=False,
        #doL1Counters=False,
        #doJetID=False,
        #genJetCollection=cms.InputTag('caGenJetsNoNu')
    #)
    #addJetCollection(
        #process,
        #cms.InputTag('akPFJetsRU','SubJets'),
        #'AKSub', 'PF',
        #doJTA=options.doJTA,
        #doBTagging=options.doBTagging,
        #btagInfo=bTagInfos,
        #btagdiscriminators=bTagDiscriminatorsSub,
        #jetCorrLabel=inputJetCorrLabelAK5,
        #doType1MET=False,
        #doL1Cleaning=False,
        #doL1Counters=False,
        #doJetID=False,
        #genJetCollection=cms.InputTag('akGenJetsNoNuRU','SubJets')
    #)
    #addJetCollection(
        #process,
        #cms.InputTag('akPFJetsFilteredCompound','SubJets'),
        #'AKFilteredSub', 'PF',
        #doJTA=options.doJTA,
        #doBTagging=options.doBTagging,
        #btagInfo=bTagInfos,
        #btagdiscriminators=bTagDiscriminatorsSub,
        #jetCorrLabel=inputJetCorrLabelAK5,
        #doType1MET=False,
        #doL1Cleaning=False,
        #doL1Counters=False,
        #doJetID=False,
        #genJetCollection=cms.InputTag('akGenJetsNoNuFiltered','SubJets')
    #)
    #addJetCollection(
        #process,
        #cms.InputTag('akPFJetsPrunedCompound','SubJets'),
        #'AKPrunedSub', 'PF',
        #doJTA=options.doJTA,
        #doBTagging=options.doBTagging,
        #btagInfo=bTagInfos,
        #btagdiscriminators=bTagDiscriminatorsSub,
        #jetCorrLabel=inputJetCorrLabelAK5,
        #doType1MET=False,
        #doL1Cleaning=False,
        #doL1Counters=False,
        #doJetID=False,
        #genJetCollection=cms.InputTag('akGenJetsNoNuPruned','SubJets')
    #)
    #addJetCollection(
        #process,
        #cms.InputTag('caPFJetsRU','SubJets'),
        #'CASub', 'PF',
        #doJTA=options.doJTA,
        #doBTagging=options.doBTagging,
        #btagInfo=bTagInfos,
        #btagdiscriminators=bTagDiscriminatorsSub,
        #jetCorrLabel=inputJetCorrLabelAK5,
        #doType1MET=False,
        #doL1Cleaning=False,
        #doL1Counters=False,
        #doJetID=False,
        #genJetCollection=cms.InputTag('caGenJetsNoNuRU','SubJets')
    #)
    #addJetCollection(
        #process,
        #cms.InputTag('caPFJetsFilteredCompound','SubJets'),
        #'CAFilteredSub', 'PF',
        #doJTA=options.doJTA,
        #doBTagging=options.doBTagging,
        #btagInfo=bTagInfos,
        #btagdiscriminators=bTagDiscriminatorsSub,
        #jetCorrLabel=inputJetCorrLabelAK5,
        #doType1MET=False,
        #doL1Cleaning=False,
        #doL1Counters=False,
        #doJetID=False,
        #genJetCollection=cms.InputTag('caGenJetsNoNuFiltered','SubJets')
    #)
    addJetCollection(
        process,
        cms.InputTag('caPFJetsPrunedCompound','SubJets'),
        'CAPrunedSub', 'PF',
        doJTA=options.doJTA,
        doBTagging=options.doBTagging,
        btagInfo=bTagInfos,
        btagdiscriminators=bTagDiscriminatorsSub,
        jetCorrLabel=inputJetCorrLabelAK5,
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False,
        genJetCollection=cms.InputTag('caGenJetsNoNuPruned','SubJets')
    )

## Define jet sequences
process.genJetSeq = cms.Sequence(
    process.akGenJetsNoNu
    + process.caGenJetsNoNu
)
process.genJetSeqExtra = cms.Sequence(
    #process.akGenJetsNoNuRU
    #+ process.akGenJetsNoNuFiltered
    #+ process.akGenJetsNoNuPruned
    #+ process.caGenJetsNoNuRU
    #+ process.caGenJetsNoNuFiltered
    process.caGenJetsNoNuPruned
)
process.recoJetSeq = cms.Sequence(
    process.akPFJets
    #+ process.akPFJetsTrimmed
    #+ process.akPFJetsFiltered
    #+ process.akPFJetsPruned
    + process.caPFJets
    + process.caPFJetsPruned
)
process.recoJetSeqExtra = cms.Sequence(
    #process.akPFJetsRU
    #+ process.akPFJetsFilteredCompound
    #+ process.akPFJetsPrunedCompound
    #+ process.caPFJetsRU
    #+ process.caPFJetsFilteredCompound
    process.caPFJetsPrunedCompound
)

## If running H tagging
if options.objectType=='H':
    process.genJetSeq = cms.Sequence( process.genJetSeq + process.genJetSeqExtra )
    process.recoJetSeq = cms.Sequence( process.recoJetSeq + process.recoJetSeqExtra )

## Define combined jet+PAT sequence
process.jetPATSequence = cms.Sequence( process.recoJetSeq + process.patDefaultSequence )

## If using explicit jet-track association
if options.useExplicitJTA:
    from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorExplicit
    for m in getattr(process,'jetPATSequence').moduleNames():
        if m.startswith('jetTracksAssociatorAtVertex'):
            print 'Switching ' + m + ' to explicit jet-track association'
            setattr( process, m, ak5JetTracksAssociatorExplicit.clone(jets = getattr(getattr(process,m),'jets')) )

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='jetPATSequence', skipAddedJets=False)

if options.doBTagging:
    # Set the cone size for the jet-track association to the jet radius
    if hasattr( process, 'jetTracksAssociatorAtVertexAKJTAPF' ):
        process.jetTracksAssociatorAtVertexAKJTAPF.coneSize = cms.double(options.jetRadius)
    if hasattr( process, 'secondaryVertexTagInfosAKJTAPF' ):
        process.secondaryVertexTagInfosAKJTAPF.trackSelection.jetDeltaRMax = cms.double(options.jetRadius)
    if hasattr( process, 'jetTracksAssociatorAtVertexCAJTAPF' ):
        process.jetTracksAssociatorAtVertexCAJTAPF.coneSize = cms.double(options.jetRadius)
    if hasattr( process, 'secondaryVertexTagInfosCAJTAPF' ):
        process.secondaryVertexTagInfosCAJTAPF.trackSelection.jetDeltaRMax = cms.double(options.jetRadius)
    # Set the jet-SV dR to the jet radius
    if hasattr( process, 'inclusiveSecondaryVertexFinderTagInfosFilteredAOD' ):
        process.inclusiveSecondaryVertexFinderTagInfosFilteredAOD.extSVDeltaRToJet = cms.double(options.jetRadius)
    if hasattr( process, 'inclusiveSecondaryVertexFinderTagInfosFilteredCAPF' ):
        process.inclusiveSecondaryVertexFinderTagInfosFilteredCAPF.extSVDeltaRToJet = cms.double(options.jetRadius)
    if hasattr( process, 'inclusiveSecondaryVertexFinderTagInfosFilteredAKJTAPF' ):
        process.inclusiveSecondaryVertexFinderTagInfosFilteredAKJTAPF.extSVDeltaRToJet = cms.double(options.jetRadius)
    if hasattr( process, 'inclusiveSecondaryVertexFinderTagInfosFilteredCAJTAPF' ):
        process.inclusiveSecondaryVertexFinderTagInfosFilteredCAJTAPF.extSVDeltaRToJet = cms.double(options.jetRadius)

## Initialize instances of the RutgersJetAnalyzer
process.jetAnalyzerDefaultJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(False),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    UseSubJets                = cms.bool(False),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsRU'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(85.),
    JetMassMax                = cms.double(150.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerTrimmedJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsRU'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerTrimmedJets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    UseGroomedJetSubstructure = cms.bool(True),
    UseSubJets                = cms.bool(False),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsRU'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerFilteredJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKFilteredPF'),
    UseSubJets                = cms.bool(False),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsRU'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(80.),
    JetMassMax                = cms.double(140.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerPrunedJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerPrunedJetMassKtSub = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsRU'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerPrunedJetMassFilteredSub = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsFilteredCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKFilteredSubPF'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerPrunedJetMassJTACone = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsAKJTAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerPrunedJetMassKtAxes = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKPrunedPF'),
    UseSubJets                = cms.bool(False),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseOnePassKtAxes          = cms.bool(False),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerPrunedJets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKPrunedPF'),
    UseGroomedJetSubstructure = cms.bool(True),
    UseSubJets                = cms.bool(False),
    GroomedBasicJetsTag       = cms.InputTag('akPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCAPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJetMassKtSub = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsRU'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCASubPF'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJetMassFilteredSub = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsFilteredCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCAFilteredSubPF'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJetMassJTACone = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAJTAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCAPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJetMass_bQuarksGSP = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCAPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(True),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(True),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJetMass_bQuarksME = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCAPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(True),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(True)
)
process.jetAnalyzerCAPrunedJetMass_cQuarks = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCAPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(True),
    JetFlavorPdgId            = cms.vint32(4),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJetMass_udsQuarks_g = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsCAPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(True),
    JetFlavorPdgId            = cms.vint32(1,2,3,21),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)
process.jetAnalyzerCAPrunedJets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsCAPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    UseGroomedJetSubstructure = cms.bool(True),
    UseSubJets                = cms.bool(False),
    GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKPrunedSubPF'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    NsubjCut                  = cms.double(0.45),
    Bdiscriminator            = cms.string("combinedSecondaryVertexBJetTags"),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgId            = cms.vint32(5),
    FindGluonSplitting        = cms.bool(False),
    FindMatrixElement         = cms.bool(False)
)

## If running over MC, add GenJets
if not options.runOnData:
    process.jetPATSequence = cms.Sequence( (process.genParticlesForJetsNoNu * process.genJetSeq) + process.jetPATSequence )

## Define jet analyzer sequence
process.jetAnalyzerSequence = cms.Sequence(
    process.jetAnalyzerDefaultJetMass
    #+ process.jetAnalyzerTrimmedJetMass
    #+ process.jetAnalyzerTrimmedJets
    #+ process.jetAnalyzerFilteredJetMass
    #+ process.jetAnalyzerPrunedJetMass
    #+ process.jetAnalyzerPrunedJetMassKtAxes
    #+ process.jetAnalyzerPrunedJets
    + process.jetAnalyzerCAPrunedJetMass
    #+ process.jetAnalyzerCAPrunedJets
)

process.jetAnalyzerSequenceExtra = cms.Sequence(
    #process.jetAnalyzerPrunedJetMassKtSub
    #+ process.jetAnalyzerPrunedJetMassFilteredSub
    #+ process.jetAnalyzerPrunedJetMassJTACone
    #+ process.jetAnalyzerCAPrunedJetMassKtSub
    #+ process.jetAnalyzerCAPrunedJetMassFilteredSub
    #+ process.jetAnalyzerCAPrunedJetMassJTACone
    process.jetAnalyzerCAPrunedJetMass_bQuarksGSP
    + process.jetAnalyzerCAPrunedJetMass_bQuarksME
    + process.jetAnalyzerCAPrunedJetMass_cQuarks
    + process.jetAnalyzerCAPrunedJetMass_udsQuarks_g
)

## If running H tagging
if options.objectType=='H':
    process.jetAnalyzerSequence = cms.Sequence( process.jetAnalyzerSequence + process.jetAnalyzerSequenceExtra )

## If running W tagging
if options.objectType=='W':
    for m in getattr(process,'jetAnalyzerSequence').moduleNames():
        setattr( getattr(process,m), 'BosonPdgId', cms.int32(24) )
        setattr( getattr(process,m), 'ApplyBosonIsolation', cms.bool(False) )
        setattr( getattr(process,m), 'BosonDecayProdPdgIds', cms.vint32(1,2,3,4,5,6) )
        setattr( getattr(process,m), 'JetPtMin', cms.double(500.) )
        setattr( getattr(process,m), 'JetPtBins', cms.uint32(2) )
        setattr( getattr(process,m), 'UseSubJets', cms.bool(False) )
        if m.startswith('jetAnalyzerDefaultJet'):
            setattr( getattr(process,m), 'JetMassMin', cms.double(70.) )
            setattr( getattr(process,m), 'JetMassMax', cms.double(125.) )
        if m.startswith('jetAnalyzerTrimmedJet'):
            setattr( getattr(process,m), 'JetMassMin', cms.double(60.) )
            setattr( getattr(process,m), 'JetMassMax', cms.double(100.) )
        if m.startswith('jetAnalyzerFilteredJet'):
            setattr( getattr(process,m), 'JetMassMin', cms.double(70.) )
            setattr( getattr(process,m), 'JetMassMax', cms.double(110.) )
        if m.startswith('jetAnalyzerPrunedJet'):
            setattr( getattr(process,m), 'JetMassMin', cms.double(55.) )
            setattr( getattr(process,m), 'JetMassMax', cms.double(95.) )
        if m.startswith('jetAnalyzerCAPrunedJet'):
            setattr( getattr(process,m), 'JetMassMin', cms.double(55.) )
            setattr( getattr(process,m), 'JetMassMax', cms.double(95.) )

## If running over background samples
if not options.runOnSignal:
    for m in getattr(process,'jetAnalyzerSequence').moduleNames():
        setattr( getattr(process,m), 'UseEventWeight', cms.bool(True) )
        setattr( getattr(process,m), 'DoBosonMatching', cms.bool(False) )

## Path definition
process.p = cms.Path(
    process.jetPATSequence
    * process.jetAnalyzerSequence
)

# Delete predefined output module (needed for running with CRAB)
del process.out
