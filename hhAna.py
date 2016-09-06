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
    fileNames   = cms.vstring (
        "root://eoscms//eos/cms//store/group/phys_higgs/resonant_HH/RunII/MicroAOD/HHbbgg_Signal76X_RegVars/1_4_0/GluGluToBulkGravitonToHHTo2B2G_M-250_narrow_13TeV-madgraph/HHbbgg_Signal76X_RegVars-1_4_0-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160407_152523/0000/myMicroAODOutputFile_1.root"
        #"microAOD-GJet.root"
        #"microAOD-ggH125.root"

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
#customize = JobConfig(metaDataSrc=, crossSections=["$CMSSW_BASE/src//MetaData/data/cross_sections.json", "/path/to/my/additional/cross_sections.json" ])
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
elif 'HH_' in customize.processId:
    sample = str(customize.processId)
else: sample='Noname'


process.HHbbggAnalyzer = cms.PSet(
    ## input specific for this analyzer

    # 1 - Default; 2 - Alternative
    cutFlow = cms.untracked.uint32(1),
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

    # In the future replace this with param._myTriggers from bbggtools:
    myTriggers=cms.untracked.vstring(
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v",
        "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v"
        ),

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
