###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
import string

options = VarParsing ('python')

options.register('runOnData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
## Make sure correct global tags are used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions)
options.register('mcGlobalTag', 'START53_V27',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', 'FT53_V21A_AN6',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
)
options.register('outFilename', 'outfile',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('reportEvery', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=10)"
)
options.register('wantSummary', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
options.register('usePFchs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PFchs"
)
options.register('doJTA', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run jet-track association"
)
options.register('useExplicitJTA', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use explicit jet-track association"
)
options.register('doBTagging', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run b tagging"
)
options.register('jetRadius', 0.8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Distance parameter R for jet clustering (default is 0.8)"
)
options.register('doBosonMatching', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do boson matching"
)
options.register('applyBosonIsolation', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply boson isolation"
)
options.register('useEventWeight', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use event weight"
)
options.register('useAK5Jets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use AK5 jets"
)
options.register('runQCDFlavorExtra', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run QCD flavor extra"
)
options.register('runOnWBkg', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on W background"
)
options.register('runOnZBkg', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on Z background"
)
options.register('runOnTopBkg', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on top background"
)
options.register('csvDiscriminator', 'combinedSecondaryVertexBJetTags',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "CSV discriminator name"
)
options.register('useSVClustering', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use SV clustering"
)
options.register('useSVMomentum', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use SV momentum"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

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
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos','inclusiveSecondaryVertexFinderTagInfos']
             #,'inclusiveSecondaryVertexFinderFilteredTagInfos','softMuonTagInfos','secondaryVertexNegativeTagInfos']
bTagDiscriminators = ['jetProbabilityBJetTags','jetBProbabilityBJetTags','combinedSecondaryVertexBJetTags'
                      #,'trackCountingHighPurBJetTags','trackCountingHighEffBJetTags'
                      #,'simpleSecondaryVertexHighPurBJetTags','simpleSecondaryVertexHighEffBJetTags'
                      ,'combinedInclusiveSecondaryVertexBJetTags']
                      #,'simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags'
                      #,'doubleSecondaryVertexHighEffBJetTags']

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

## Geometry and Detector Conditions
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag + '::All'

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

#-------------------------------------
## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:patTuple_PF2PAT_v2.root'
        # 'file:/cms/ferencek/store/ferencek/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_2_NZA.root'
        'file:/cms/ferencek/store/skaplan/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_1_bnp.root'
        #'file:/cms/ferencek/store/skaplan/BprimeBprimeToBHBHinc_M-1500_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_1_djM.root'
        #'file:/cms/ferencek/store/ferencek/TprimeToBWinc_M-1000_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_1_9ZZ.root'
        #'file:/cms/ferencek/store/ferencek/BprimeBprimeToBZBZinc_M-1200_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_1_vSQ.root'
        #'file:/cms/ferencek/store/ferencek/TprimeToTHinc_M-1700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_1_gMe.root'
        #'file:/cms/ferencek/store/ferencek/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_1_mRj.root'
        #'file:/cms/ferencek/store/ferencek/RadionToHHTo4B_M-1500_TuneZ2star_8TeV-nm-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_PATTuple_v2/6950c4b6452599a829ed09f6192c8cf5/patTuple_PF2PAT_v2_1_1_fdQ.root'
    )
)

#-------------------------------------
outFilename = string.replace(options.outFilename,'.root','') + '_mc.root'
if options.runOnData :
    outFilename = string.replace(options.outFilename,'.root','') + '_data.root'

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(outFilename)
)
#-------------------------------------
## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)

#-------------------------------------
## CA8 jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNu = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJetsCHS = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("pfNoElectronPFlow"),
    srcPVs = cms.InputTag("goodOfflinePrimaryVertices"),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(20.)
)
## CA8 filtered jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNuFiltered = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    useFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.ca8PFJetsCHSFiltered = ak5PFJetsFiltered.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = process.ca8PFJetsCHS.src,
    srcPVs = process.ca8PFJetsCHS.srcPVs,
    doAreaFastjet = process.ca8PFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.)
)
## CA8 MassDrop-BDRS filtered jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
## Compared to the above filtered jets, here dynamic filtering radius is used (as in arXiv:0802.2470)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNuMDBDRSFiltered = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    useMassDropTagger = cms.bool(True),
    muCut = cms.double(0.667),
    yCut = cms.double(0.08),
    useFiltering = cms.bool(True),
    useDynamicFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    rFiltFactor = cms.double(0.5),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsMassDropFiltered
process.ca8PFJetsCHSMDBDRSFiltered = ak5PFJetsMassDropFiltered.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = process.ca8PFJetsCHS.src,
    srcPVs = process.ca8PFJetsCHS.srcPVs,
    doAreaFastjet = process.ca8PFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.),
    useDynamicFiltering = cms.bool(True),
    rFiltFactor = cms.double(0.5)
)
## CA8 Kt-BDRS filtered jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
## Compared to the above filtered jets, here dynamic filtering radius is used (as in arXiv:0802.2470)
## However, here the mass drop is replaced by finding two Kt subjets which then set the size of the filtering radius
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNuKtBDRSFiltered = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    useKtPruning = cms.bool(True),
    zcut = cms.double(0.),
    rcut_factor = cms.double(9999.),
    useFiltering = cms.bool(True),
    useDynamicFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    rFiltFactor = cms.double(0.5),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.ca8PFJetsCHSKtBDRSFiltered = ak5PFJetsFiltered.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = process.ca8PFJetsCHS.src,
    srcPVs = process.ca8PFJetsCHS.srcPVs,
    doAreaFastjet = process.ca8PFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.),
    usePruning = cms.bool(True),
    useKtPruning = cms.bool(True),
    zcut = cms.double(0.),
    rcut_factor = cms.double(9999.),
    useDynamicFiltering = cms.bool(True),
    rFiltFactor = cms.double(0.5)
)
## CA8 pruned jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ca8GenJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ca8PFJetsCHSPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = process.ca8PFJetsCHS.src,
    srcPVs = process.ca8PFJetsCHS.srcPVs,
    doAreaFastjet = process.ca8PFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.)
)
## CA8 jets with Kt subjets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
## Kt subjets produced using Kt-based pruning with very loose pruning cuts (pruning is effectively disabled)
process.ca8GenJetsNoNuKtPruned = ca4GenJets.clone(
    SubJetParameters.clone(
        zcut = cms.double(0.),
        rcut_factor = cms.double(9999.)
    ),
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    useKtPruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ca8PFJetsCHSKtPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = process.ca8PFJetsCHS.src,
    srcPVs = process.ca8PFJetsCHS.srcPVs,
    doAreaFastjet = process.ca8PFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.),
    useKtPruning = cms.bool(True),
    zcut = cms.double(0.),
    rcut_factor = cms.double(9999.)
)
## CA8 trimmed jets (Reco only)
from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
process.ca8PFJetsCHSTrimmed = ak5PFJetsTrimmed.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = process.ca8PFJetsCHS.src,
    srcPVs = process.ca8PFJetsCHS.srcPVs,
    doAreaFastjet = process.ca8PFJetsCHS.doAreaFastjet,
    jetPtMin = cms.double(20.)
)

#-------------------------------------
## PATify the above jets
from PhysicsTools.PatAlgos.tools.jetTools import *
## CA8 jets
switchJetCollection(process,
    cms.InputTag('ca8PFJetsCHS'),
    jetIdLabel='ca8',
    rParam = 0.8,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel = inputJetCorrLabelAK5,
    doType1MET   = False,
    genJetCollection = cms.InputTag("ca8GenJetsNoNu"),
    doJetID      = False,
)
## Filtered CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSFiltered'),
    'CA8','FilteredPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("ca8GenJetsNoNu")
)
## Filtered subjets of CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSFiltered','SubJets'),
    'CA8', 'FilteredSubjetsPFCHS',
    rParam = 0.8,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNuFiltered','SubJets')
)
## MassDrop-BDRS filtered CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSMDBDRSFiltered'),
    'CA8','MDBDRSFilteredPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("ca8GenJetsNoNu")
)
## MassDrop-BDRS filtered subjets of CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSMDBDRSFiltered','SubJets'),
    'CA8', 'MDBDRSFilteredSubjetsPFCHS',
    rParam = 0.8,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNuMDBDRSFiltered','SubJets')
)
## Kt-BDRS filtered CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSKtBDRSFiltered'),
    'CA8','KtBDRSFilteredPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("ca8GenJetsNoNu")
)
## Kt-BDRS filtered subjets of CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSKtBDRSFiltered','SubJets'),
    'CA8', 'KtBDRSFilteredSubjetsPFCHS',
    rParam = 0.8,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNuKtBDRSFiltered','SubJets')
)
## Pruned CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSPruned'),
    'CA8','PrunedPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("ca8GenJetsNoNu")
)
## Pruned subjets of CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSPruned','SubJets'),
    'CA8', 'PrunedSubjetsPFCHS',
    rParam = 0.8,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNuPruned','SubJets')
)
## Kt pruned CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSKtPruned'),
    'CA8','KtPrunedPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("ca8GenJetsNoNu")
)
## Kt subjets of CA8 jets
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsCHSKtPruned','SubJets'),
    'CA8', 'KtSubjetsPFCHS',
    rParam = 0.8,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNuKtPruned','SubJets')
)

#-------------------------------------
## N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

process.NjettinessCA8 = Njettiness.clone(
    src = cms.InputTag("ca8PFJetsCHS"),
    cone = cms.double(0.8)
)

process.patJets.userData.userFloats.src += ['NjettinessCA8:tau1','NjettinessCA8:tau2','NjettinessCA8:tau3']

#-------------------------------------
## Grooming ValueMaps
from RecoJets.JetProducers.ca8PFJetsCHS_groomingValueMaps_cfi import ca8PFJetsCHSPrunedLinks

process.ca8PFJetsCHSPrunedMass = ca8PFJetsCHSPrunedLinks.clone(
    src = cms.InputTag("ca8PFJetsCHS"),
    matched = cms.InputTag("ca8PFJetsCHSPruned"),
    distMax = cms.double(0.8),
    value = cms.string('mass')
)

process.ca8PFJetsCHSFilteredMass = ca8PFJetsCHSPrunedLinks.clone(
    src = cms.InputTag("ca8PFJetsCHS"),
    matched = cms.InputTag("ca8PFJetsCHSFiltered"),
    distMax = cms.double(0.8),
    value = cms.string('mass')
)

process.ca8PFJetsCHSTrimmedMass = ca8PFJetsCHSPrunedLinks.clone(
    src = cms.InputTag("ca8PFJetsCHS"),
    matched = cms.InputTag("ca8PFJetsCHSTrimmed"),
    distMax = cms.double(0.8),
    value = cms.string('mass')
)

process.patJets.userData.userFloats.src += ['ca8PFJetsCHSPrunedMass','ca8PFJetsCHSFilteredMass','ca8PFJetsCHSTrimmedMass']

#-------------------------------------
if options.useSVClustering:
    ## Enable clustering-based jet-SV association for IVF vertices and AK5 jets
    #process.inclusiveSecondaryVertexFinderTagInfosAODPFlow = process.inclusiveSecondaryVertexFinderTagInfosAODPFlow.clone(
        #useSVClustering = cms.bool(True),
        ##useSVMomentum   = cms.bool(True), # otherwise using SV flight direction
        #jetAlgorithm    = cms.string("AntiKt"),
        #rParam          = cms.double(0.5),
        #ghostRescaling  = cms.double(1e-18)
    #)
    ## Enable clustering-based jet-SV association for IVF vertices and subjets of CA8 jets
    process.inclusiveSecondaryVertexFinderTagInfosCA8FilteredSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCA8FilteredSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(0.8),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("ca8PFJetsCHS"),
        groomedFatJets   = cms.InputTag("ca8PFJetsCHSFiltered")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCA8MDBDRSFilteredSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCA8MDBDRSFilteredSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(0.8),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("ca8PFJetsCHS"),
        groomedFatJets   = cms.InputTag("ca8PFJetsCHSMDBDRSFiltered")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCA8KtBDRSFilteredSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCA8KtBDRSFilteredSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(0.8),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("ca8PFJetsCHS"),
        groomedFatJets   = cms.InputTag("ca8PFJetsCHSKtBDRSFiltered")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCA8PrunedSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCA8PrunedSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(0.8),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("ca8PFJetsCHS"),
        groomedFatJets   = cms.InputTag("ca8PFJetsCHSPruned")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCA8KtSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCA8KtSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(0.8),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("ca8PFJetsCHS"),
        groomedFatJets   = cms.InputTag("ca8PFJetsCHSKtPruned")
    )
#-------------------------------------
## New jet flavor still requires some cfg-level adjustments for subjets until it is better integrated into PAT
## Adjust the jet flavor for CA8 filtered subjets
process.patJetFlavourAssociationCA8FilteredSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("ca8PFJetsCHSFiltered"),
    subjets = cms.InputTag("ca8PFJetsCHSFiltered", "SubJets")
)
process.patJetsCA8FilteredSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCA8FilteredSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA8 MassDrop-BDRS filtered subjets
process.patJetFlavourAssociationCA8MDBDRSFilteredSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("ca8PFJetsCHSMDBDRSFiltered"),
    subjets = cms.InputTag("ca8PFJetsCHSMDBDRSFiltered", "SubJets")
)
process.patJetsCA8MDBDRSFilteredSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCA8MDBDRSFilteredSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA8 Kt-BDRS filtered subjets
process.patJetFlavourAssociationCA8KtBDRSFilteredSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("ca8PFJetsCHSKtBDRSFiltered"),
    subjets = cms.InputTag("ca8PFJetsCHSKtBDRSFiltered", "SubJets")
)
process.patJetsCA8KtBDRSFilteredSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCA8KtBDRSFilteredSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA8 pruned subjets
process.patJetFlavourAssociationCA8PrunedSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("ca8PFJetsCHSPruned"),
    subjets = cms.InputTag("ca8PFJetsCHSPruned", "SubJets")
)
process.patJetsCA8PrunedSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCA8PrunedSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA8 Kt subjets
process.patJetFlavourAssociationCA8KtSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("ca8PFJetsCHSKtPruned"),
    subjets = cms.InputTag("ca8PFJetsCHSKtPruned", "SubJets")
)
process.patJetsCA8KtSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCA8KtSubjetsPFCHS","SubJets")

#-------------------------------------
## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsCA8FilteredPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCA8FilteredPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCA8FilteredSubjetsPFCHS")
)

process.selectedPatJetsCA8MDBDRSFilteredPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCA8MDBDRSFilteredPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCA8MDBDRSFilteredSubjetsPFCHS")
)

process.selectedPatJetsCA8KtBDRSFilteredPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCA8KtBDRSFilteredPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCA8KtBDRSFilteredSubjetsPFCHS")
)

process.selectedPatJetsCA8PrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCA8PrunedPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCA8PrunedSubjetsPFCHS")
)

process.selectedPatJetsCA8KtPrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCA8KtPrunedPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCA8KtSubjetsPFCHS")
)

## Define BoostedJetMerger sequence
process.jetMergerSeq = cms.Sequence(
    process.selectedPatJetsCA8FilteredPFCHSPacked
    + process.selectedPatJetsCA8MDBDRSFilteredPFCHSPacked
    + process.selectedPatJetsCA8KtBDRSFilteredPFCHSPacked
    + process.selectedPatJetsCA8PrunedPFCHSPacked
    + process.selectedPatJetsCA8KtPrunedPFCHSPacked
)

#-------------------------------------
from PhysicsTools.PatAlgos.tools.coreTools import *
## Keep only jets in the default sequence
removeAllPATObjectsBut(process, ['Jets'])

## Remove MC matching when running over data
if options.runOnData:
    removeMCMatching( process, ['All'] )

#-------------------------------------
## GenParticles for GenJets
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu = genParticlesForJetsNoNu

#-------------------------------------
## Define instances of the RutgersJetAnalyzer
process.jetAnalyzerCA8FatJets_PrunedSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCA8PrunedPFCHSPacked'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(True),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    Bdiscriminator            = cms.string(options.csvDiscriminator),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCA8FatJets_FilteredSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCA8FilteredPFCHSPacked'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    Bdiscriminator            = cms.string(options.csvDiscriminator),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCA8MDBDRSFilteredPFCHSPacked'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    Bdiscriminator            = cms.string(options.csvDiscriminator),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCA8KtBDRSFilteredPFCHSPacked'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    Bdiscriminator            = cms.string(options.csvDiscriminator),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCA8FatJets_KtSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCA8KtPrunedPFCHSPacked'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    Bdiscriminator            = cms.string(options.csvDiscriminator),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)

## Define jet analyzer sequence
process.jetAnalyzerSequence = cms.Sequence(
    process.jetAnalyzerCA8FatJets_PrunedSubjets
    + process.jetAnalyzerCA8FatJets_FilteredSubjets
    + process.jetAnalyzerCA8FatJets_MDBDRSFilteredSubjets
    + process.jetAnalyzerCA8FatJets_KtBDRSFilteredSubjets
    + process.jetAnalyzerCA8FatJets_KtSubjets
)

#-------------------------------------
## If using explicit jet-track association
if options.useExplicitJTA:
    from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorExplicit
    for m in getattr(process,"patDefaultSequence").moduleNames():
        if m.startswith('jetTracksAssociatorAtVertex'):
            print 'Switching ' + m + ' to explicit jet-track association'
            setattr( process, m, ak5JetTracksAssociatorExplicit.clone(jets = getattr(getattr(process,m),'jets')) )

#-------------------------------------
## Remove tau stuff that really shouldn't be there (probably a bug in PAT)
process.patDefaultSequence.remove(process.kt6PFJetsForRhoComputationVoronoi)
for m in getattr(process,"patDefaultSequence").moduleNames():
    if m.startswith('hpsPFTau'):
        getattr(process,"patDefaultSequence").remove(getattr(process,m))

#-------------------------------------
## Define jet sequences
process.genJetSeq = cms.Sequence(
    process.ca8GenJetsNoNu
    + process.ca8GenJetsNoNuFiltered
    + process.ca8GenJetsNoNuMDBDRSFiltered
    + process.ca8GenJetsNoNuKtBDRSFiltered
    + process.ca8GenJetsNoNuPruned
    + process.ca8GenJetsNoNuKtPruned
)
process.jetSeq = cms.Sequence(
    (
    process.ca8PFJetsCHS
    + process.ca8PFJetsCHSFiltered
    + process.ca8PFJetsCHSMDBDRSFiltered
    + process.ca8PFJetsCHSKtBDRSFiltered
    + process.ca8PFJetsCHSPruned
    + process.ca8PFJetsCHSKtPruned
    + process.ca8PFJetsCHSTrimmed
    )
    * (
    process.NjettinessCA8
    + process.ca8PFJetsCHSFilteredMass
    + process.ca8PFJetsCHSPrunedMass
    + process.ca8PFJetsCHSTrimmedMass
    )
)

if not options.runOnData:
    process.jetSeq = cms.Sequence( process.genParticlesForJetsNoNu * process.genJetSeq + process.jetSeq )

#-------------------------------------
## Adapt primary vertex collection
from PhysicsTools.PatAlgos.tools.pfTools import *
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='jetSeq')
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='patDefaultSequence')

#-------------------------------------
## Add TagInfos to PAT jets
for m in ['patJets', 'patJetsCA8FilteredSubjetsPFCHS', 'patJetsCA8MDBDRSFilteredSubjetsPFCHS', 'patJetsCA8KtBDRSFilteredSubjetsPFCHS',
          'patJetsCA8PrunedSubjetsPFCHS', 'patJetsCA8KtSubjetsPFCHS']:
    if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )

#-------------------------------------
## Adapt fat jet b tagging
if options.doBTagging:
    # Set the cone size for the jet-track association to the jet radius
    process.jetTracksAssociatorAtVertex.coneSize = cms.double(options.jetRadius) # default is 0.5
    process.secondaryVertexTagInfosAOD.trackSelection.jetDeltaRMax = cms.double(options.jetRadius)   # default is 0.3
    process.secondaryVertexTagInfosAOD.vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.5
    # Set the jet-SV dR to the jet radius
    process.inclusiveSecondaryVertexFinderTagInfosAOD.extSVDeltaRToJet = cms.double(options.jetRadius) # default is 0.3
    # Set the JP track dR cut to the jet radius
    process.jetProbabilityCA8 = process.jetProbability.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.3
    process.jetProbabilityBJetTagsAOD.jetTagComputer = cms.string('jetProbabilityCA8')
    # Set the JBP track dR cut to the jet radius
    process.jetBProbabilityCA8 = process.jetBProbability.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.5
    process.jetBProbabilityBJetTagsAOD.jetTagComputer = cms.string('jetBProbabilityCA8')
    # Set the CSV track dR cut to the jet radius
    process.combinedSecondaryVertexCA8 = process.combinedSecondaryVertex.clone()
    process.combinedSecondaryVertexCA8.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexCA8.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexBJetTagsAOD.jetTagComputer = cms.string('combinedSecondaryVertexCA8')
    process.combinedInclusiveSecondaryVertexBJetTagsAOD.jetTagComputer = cms.string('combinedSecondaryVertexCA8')

#-------------------------------------
## Various additional options
for m in getattr(process,'jetAnalyzerSequence').moduleNames():
    if not options.doBosonMatching:
        setattr( getattr(process,m), 'DoBosonMatching', cms.bool(False) )
    if not options.applyBosonIsolation:
        setattr( getattr(process,m), 'ApplyBosonIsolation', cms.bool(False) )
    if options.useEventWeight:
        setattr( getattr(process,m), 'UseEventWeight', cms.bool(True) )

#-------------------------------------
## Path definition
process.p = cms.Path(
    ( process.jetSeq * process.patDefaultSequence )
    * process.jetMergerSeq
    * process.jetAnalyzerSequence
)

## Delete output module
del process.out

## Schedule definition
process.schedule = cms.Schedule(process.p)
