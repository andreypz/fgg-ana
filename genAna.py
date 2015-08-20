import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils


import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.MetaData.samples_utils import SamplesManager


process = cms.Process('GenAna')

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        #"root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/mdonega/flashgg/RunIISpring15-50ns/Spring15BetaV2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15-50ns-Spring15BetaV2-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150716_155016/0000/myMicroAODOutputFile_81.root",
        #"microAOD-GJet.root"
        #"microAOD-ggH125.root"
        "root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15-50ns-Spring15BetaV4-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150813_120650/0000/myMicroAODOutputFile_1.root"

), ## mandatory
    maxEvents   = cms.int32(-1),                            ## optional
    outputEvery = cms.uint32(10000),                            ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('a_Histograms.root'),  ## mandatory
)

process.GenAnalyzer = cms.PSet(
    ## input specific for this analyzer
    lep  = cms.untracked.string('el'), # "mu" or "el"
    myMuons = cms.InputTag('flashggSelectedMuons'),
    myElectrons = cms.InputTag('flashggSelectedElectrons'),
    myPhotons = cms.InputTag('flashggPhotons'),
    myGens = cms.InputTag('flashggPrunedGenParticles'),
    lumiWeight = cms.double(1.),
)

from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents",200)
customize.setDefault("targetLumi",10e+3)
customize.crossSections.append("$CMSSW_BASE/src/APZ/fgg-ana/cross_sections.json")

customize(process)

