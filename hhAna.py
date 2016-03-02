import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

process = cms.Process('HHbbggAna')

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        #"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/mdonega/flashgg/RunIISpring15-50ns/Spring15BetaV2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150716_155016/0000/myMicroAODOutputFile_81.root",
        #"microAOD-GJet.root"
        #"microAOD-ggH125.root"
        "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15-50ns-Spring15BetaV4-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150813_120650/0000/myMicroAODOutputFile_1.root"

), ## mandatory
    maxEvents   = cms.int32(-1),                            ## optional
    outputEvery = cms.uint32(40000),                        ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('a_Histograms.root'),  ## mandatory
)


import flashgg.bbggTools.parameters as param
import flashgg.Taggers.flashggTags_cff as flashggTags

from flashgg.MetaData.JobConfig import customize
customize.crossSections.append("$CMSSW_BASE/src/APZ/fgg-ana/data/cross_sections.json")
customize.setDefault("maxEvents",200)
customize.setDefault("targetLumi",10e+3)
customize.parse()

print "Process ID = ", bcolors.OKBLUE + customize.processId + bcolors.ENDC

sample = 'None'
if 'QCD' in customize.processId:
    sample = 'QCD'
elif 'GJet' in customize.processId:
    sample = 'GJet'
elif 'DYJets' in customize.processId:
    sample = 'DYJets'
elif 'DiPhoton' in customize.processId:
    sample = 'DiPhoton'
elif 'DoubleEG' in customize.processId:
    sample = 'DoubleEG'
elif 'Graviton' in customize.processId:
    sample = 'Graviton'
elif 'Radion' in customize.processId:
    sample = 'Radion'
else: sample='Nonono'


process.HHbbggAnalyzer = cms.PSet(
    ## input specific for this analyzer

    # 1 - NCU; 2 - Rafael
    cutFlow = cms.untracked.uint32(2),
    useDiPhotons = cms.untracked.bool(True),
    diPhotonTag  = cms.InputTag('flashggDiPhotons'),

    # 1 - Cut Based XX WP; 2 Cut based from PAT object; 3 - MVA ID from Egamma;
    # 4 - Hgg MVA ID; 5 - High Pt ID to be implementd (but not sure if it's available)
    # (Only used if useDiPhotons==False...)
    phoIDtype = cms.untracked.uint32(3),

    #phoISOcutEB=param._phoISOlooseEB,
    #phoISOcutEE=param._phoISOlooseEE,
    #phoIDcutEB =param._phoIDlooseEB,
    #phoIDcutEE =param._phoIDlooseEE,

    phoISOcutEB=param._phoISOmediumEB,
    phoISOcutEE=param._phoISOmediumEE,
    phoIDcutEB =param._phoIDmediumEB,
    phoIDcutEE =param._phoIDmediumEE,

    lep  = cms.untracked.string('el'), # "mu" or "el"
    jetTag      = cms.InputTag('flashggFinalJets'),
    muonTag     = cms.InputTag('flashggSelectedMuons'),
    electronTag = cms.InputTag('flashggSelectedElectrons'),
    photonTag   = cms.InputTag('flashggPhotons'),
    #genTag      = cms.InputTag('prunedGenParticles'),
    genTag      = cms.InputTag('flashggPrunedGenParticles'),
    #inputTagJets= flashggTags.UnpackedJetCollectionVInputTag,
    lumiWeight = cms.double(1.),

    runSample = cms.untracked.string(sample),

)

customize(process)
