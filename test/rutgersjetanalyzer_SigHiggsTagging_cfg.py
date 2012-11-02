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
options.register('jetRadius',
    0.8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Distance parameter R for jet clustering (default is 0.5)"
)
options.register('doBTagging',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run b tagging"
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

## b-tagging
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos','inclusiveSecondaryVertexFinderTagInfos','softMuonTagInfos','softElectronTagInfos','inclusiveSecondaryVertexFinderTagInfosFiltered']
bTagDiscriminators = ['jetProbabilityBJetTags','trackCountingHighPurBJetTags', 'trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags',
                      'simpleSecondaryVertexHighPurBJetTags','simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags',
                      'combinedSecondaryVertexBJetTags','combinedSecondaryVertexMVABJetTags','combinedInclusiveSecondaryVertexBJetTags',
                      'softMuonBJetTags','softMuonByPtBJetTags','softMuonByIP3dBJetTags','combinedMVABJetTags','doubleSecondaryVertexHighEffBJetTags']
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
    src = cms.InputTag("pfNoElectronPFlow")
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
#from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
#process.akPFJetsFiltered = ak5PFJetsFiltered.clone(
    #rParam = process.akPFJets.rParam,
    #src = process.akPFJets.src,
    #srcPVs = process.akPFJets.srcPVs,
    #doAreaFastjet = process.akPFJets.doAreaFastjet,
    #writeCompound = cms.bool(False),
    #jetCollInstanceName=cms.string("")
#)
#from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
#process.akPFJetsPruned = ak5PFJetsPruned.clone(
    #rParam = process.akPFJets.rParam,
    #src = process.akPFJets.src,
    #srcPVs = process.akPFJets.srcPVs,
    #doAreaFastjet = process.akPFJets.doAreaFastjet,
    #writeCompound = cms.bool(False),
    #jetCollInstanceName=cms.string("")
#)

### Define CA jets (GEN and RECO)
#from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
#process.caGenJetsNoNu = ca4GenJets.clone(
    #rParam = options.jetRadius,
    #src = cms.InputTag("genParticlesForJetsNoNu")
#)
#from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
#process.caPFJetsPruned = ak5PFJetsPruned.clone(
    #jetAlgorithm = cms.string("CambridgeAachen"),
    #rParam = options.jetRadius,
    #src = process.akPFJets.src,
    #srcPVs = process.akPFJets.srcPVs,
    #doAreaFastjet = process.akPFJets.doAreaFastjet,
#)

## PATify above jets
from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection
switchJetCollection(process,
    cms.InputTag('akPFJets'),
    doJTA=True,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    genJetCollection=cms.InputTag("akGenJetsNoNu"),
    doJetID=False
)
addJetCollection(
    process,
    cms.InputTag('akPFJets'),
    'AK', 'PF',
    doJTA=True,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('akGenJetsNoNu')
)
addJetCollection(
    process,
    cms.InputTag('akPFJetsRU','SubJets'),
    'AKSub', 'PF',
    doJTA=True,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminatorsSub,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('akGenJetsNoNuRU','SubJets')
)
addJetCollection(
    process,
    cms.InputTag('akPFJetsTrimmed'),
    'AKTrimmed','PF',
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("akGenJetsNoNu")
)
#addJetCollection(
    #process,
    #cms.InputTag('akPFJetsFiltered'),
    #'AKFiltered','PF',
    #doJTA=True,
    #doBTagging=options.doBTagging,
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
    #doJTA=True,
    #doBTagging=options.doBTagging,
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
    #cms.InputTag('caPFJetsPruned'),
    #'CAPruned','PF',
    #doJTA=True,
    #doBTagging=options.doBTagging,
    #btagInfo=bTagInfos,
    #btagdiscriminators=bTagDiscriminators,
    #jetCorrLabel=inputJetCorrLabelAK7,
    #doType1MET=False,
    #doL1Cleaning=False,
    #doL1Counters=False,
    #doJetID=False,
    #genJetCollection=cms.InputTag("caGenJetsNoNu")
#)

## Define a sequence for jets
process.jetSeq = cms.Sequence(
    process.genParticlesForJetsNoNu *
    (
    process.akGenJetsNoNu
    + process.akGenJetsNoNuRU
    #+ process.caGenJetsNoNu
    )
    + process.akPFJets
    + process.akPFJetsRU
    + process.akPFJetsTrimmed
    #+ process.akPFJetsFiltered
    #+ process.akPFJetsPruned
    #+ process.caPFJetsPruned
)

## If running over data, remove GenJets
if options.runOnData:
    process.jetSeq.remove(process.genParticlesForJetsNoNu)
    process.jetSeq.remove(process.akGenJetsNoNu)
    process.jetSeq.remove(process.akGenJetsNoNuRU)
    process.jetSeq.remove(process.caGenJetsNoNu)

process.jetPATSequence = cms.Sequence( process.jetSeq + process.patDefaultSequence)

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='jetPATSequence', skipAddedJets=False)

if options.doBTagging:
    # Set the cone size for the jet-track association to the jet radius
    process.jetTracksAssociatorAtVertexAKPF.coneSize = cms.double(options.jetRadius)
    # Set the jet-SV dR to the jet radius
    process.inclusiveSecondaryVertexFinderTagInfosFilteredAOD.extSVDeltaRToJet = cms.double(options.jetRadius)
    process.inclusiveSecondaryVertexFinderTagInfosFilteredAKPF.extSVDeltaRToJet = cms.double(options.jetRadius)

## Initialize instances of the RutgersJetAnalyzer
#process.jetAnalyzerDefaultJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    #UseEventWeight            = cms.bool(False),
    #GenParticleTag            = cms.InputTag('genParticles'),
    #JetsTag                   = cms.InputTag('selectedPatJets'),
    #UseGroomedJets            = cms.bool(False),
    #GroomedJetsTag            = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    #UseSubJets                = cms.bool(False),
    #SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    #PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    #JetRadius                 = cms.double(options.jetRadius),
    #DoBosonMatching           = cms.bool(True),
    #BosonMatchingRadius       = cms.double(0.5),
    #BosonPdgId                = cms.int32(25),
    #DoBosonDecayProdSelection = cms.bool(True),
    #BosonDecayProdPdgIds      = cms.vint32(5),
    #UseMassDrop               = cms.bool(False),
    #JetPtMin                  = cms.double(300.),
    #JetPtBins                 = cms.uint32(3),
    #JetPtBinWidth             = cms.double(200.),
    #JetAbsEtaMax              = cms.double(1.5),
    #JetMassMin                = cms.double(90.),
    #JetMassMax                = cms.double(150.),
    #NsubjCut                  = cms.double(0.45)
#)
process.jetAnalyzerTrimmedJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    UseSubJets                = cms.bool(True),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(80.),
    JetMassMax                = cms.double(130.),
    NsubjCut                  = cms.double(0.45)
)
process.jetAnalyzerTrimmedJetMassJTACone = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJetsAKPF'),
    UseGroomedJets            = cms.bool(True),
    GroomedJetsTag            = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    UseSubJets                = cms.bool(True),
    SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    UseMassDrop               = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(80.),
    JetMassMax                = cms.double(130.),
    NsubjCut                  = cms.double(0.45)
)
#process.jetAnalyzerTrimmedJets = cms.EDAnalyzer('RutgersJetAnalyzer',
    #UseEventWeight            = cms.bool(False),
    #GenParticleTag            = cms.InputTag('genParticles'),
    #JetsTag                   = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    #UseGroomedJets            = cms.bool(False),
    #GroomedJetsTag            = cms.InputTag('selectedPatJetsAKTrimmedPF'),
    #UseSubJets                = cms.bool(False),
    #SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    #PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    #JetRadius                 = cms.double(options.jetRadius),
    #DoBosonMatching           = cms.bool(True),
    #BosonMatchingRadius       = cms.double(0.5),
    #BosonPdgId                = cms.int32(25),
    #DoBosonDecayProdSelection = cms.bool(True),
    #BosonDecayProdPdgIds      = cms.vint32(5),
    #UseMassDrop               = cms.bool(False),
    #JetPtMin                  = cms.double(300.),
    #JetPtBins                 = cms.uint32(3),
    #JetPtBinWidth             = cms.double(200.),
    #JetAbsEtaMax              = cms.double(1.5),
    #JetMassMin                = cms.double(80.),
    #JetMassMax                = cms.double(130.),
    #NsubjCut                  = cms.double(0.45)
#)
#process.jetAnalyzerFilteredJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    #UseEventWeight            = cms.bool(False),
    #GenParticleTag            = cms.InputTag('genParticles'),
    #JetsTag                   = cms.InputTag('selectedPatJets'),
    #UseGroomedJets            = cms.bool(True),
    #GroomedJetsTag            = cms.InputTag('selectedPatJetsAKFilteredPF'),
    #UseSubJets                = cms.bool(False),
    #SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    #PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    #JetRadius                 = cms.double(options.jetRadius),
    #DoBosonMatching           = cms.bool(True),
    #BosonMatchingRadius       = cms.double(0.5),
    #BosonPdgId                = cms.int32(25),
    #DoBosonDecayProdSelection = cms.bool(True),
    #BosonDecayProdPdgIds      = cms.vint32(5),
    #UseMassDrop               = cms.bool(False),
    #JetPtMin                  = cms.double(300.),
    #JetPtBins                 = cms.uint32(3),
    #JetPtBinWidth             = cms.double(200.),
    #JetAbsEtaMax              = cms.double(1.5),
    #JetMassMin                = cms.double(80.),
    #JetMassMax                = cms.double(140.),
    #NsubjCut                  = cms.double(0.45)
#)
#process.jetAnalyzerPrunedJetMass = cms.EDAnalyzer('RutgersJetAnalyzer',
    #UseEventWeight            = cms.bool(False),
    #GenParticleTag            = cms.InputTag('genParticles'),
    #JetsTag                   = cms.InputTag('selectedPatJets'),
    #UseGroomedJets            = cms.bool(True),
    #GroomedJetsTag            = cms.InputTag('selectedPatJetsAKPrunedPF'),
    #UseSubJets                = cms.bool(False),
    #SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    #PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    #JetRadius                 = cms.double(options.jetRadius),
    #DoBosonMatching           = cms.bool(True),
    #BosonMatchingRadius       = cms.double(0.5),
    #BosonPdgId                = cms.int32(25),
    #DoBosonDecayProdSelection = cms.bool(True),
    #BosonDecayProdPdgIds      = cms.vint32(5),
    #UseMassDrop               = cms.bool(False),
    #JetPtMin                  = cms.double(300.),
    #JetPtBins                 = cms.uint32(3),
    #JetPtBinWidth             = cms.double(200.),
    #JetAbsEtaMax              = cms.double(1.5),
    #JetMassMin                = cms.double(75.),
    #JetMassMax                = cms.double(130.),
    #NsubjCut                  = cms.double(0.45)
#)
#process.jetAnalyzerCAPrunedJets = cms.EDAnalyzer('RutgersJetAnalyzer',
    #UseEventWeight            = cms.bool(False),
    #GenParticleTag            = cms.InputTag('genParticles'),
    #JetsTag                   = cms.InputTag('selectedPatJetsCAPrunedPF'),
    #UseGroomedJets            = cms.bool(False),
    #GroomedJetsTag            = cms.InputTag('selectedPatJetsCAPrunedPF'),
    #UseSubJets                = cms.bool(False),
    #SubJetsTag                = cms.InputTag('selectedPatJetsAKSubPF'),
    #PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    #JetRadius                 = cms.double(options.jetRadius),
    #DoBosonMatching           = cms.bool(True),
    #BosonMatchingRadius       = cms.double(0.5),
    #BosonPdgId                = cms.int32(25),
    #DoBosonDecayProdSelection = cms.bool(True),
    #BosonDecayProdPdgIds      = cms.vint32(5),
    #UseMassDrop               = cms.bool(True),
    #JetPtMin                  = cms.double(300.),
    #JetPtBins                 = cms.uint32(3),
    #JetPtBinWidth             = cms.double(200.),
    #JetAbsEtaMax              = cms.double(1.5),
    #JetMassMin                = cms.double(75.),
    #JetMassMax                = cms.double(130.),
    #NsubjCut                  = cms.double(0.45)
#)

## Path definition
process.p = cms.Path(
    process.jetPATSequence
    * (
    #process.jetAnalyzerDefaultJetMass
    process.jetAnalyzerTrimmedJetMass
    + process.jetAnalyzerTrimmedJetMassJTACone
    #+ process.jetAnalyzerTrimmedJets
    #+ process.jetAnalyzerFilteredJetMass
    #+ process.jetAnalyzerPrunedJetMass
    #+ process.jetAnalyzerCAPrunedJets
    )
)

# Delete predefined output module (needed for running with CRAB)
del process.out
