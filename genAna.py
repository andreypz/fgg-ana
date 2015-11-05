import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils


import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.MetaData.samples_utils import SamplesManager


process = cms.Process('GenAna')

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
    #"microAOD-ggH125.root"
    "root://eoscms//eos/cms//store/group/phys_higgs/cmshgg/rateixei/flashgg/RunIISpring15-25ns/Spring15BetaV5/GluGluToBulkGravitonToHHTo2B2G_M-260_narrow_13TeV-madgraph/RunIISpring15-25ns-Spring15BetaV5-v0-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/151013_094152/0000/myMicroAODOutputFile_1.root"
    #"root://cmsxrootd.fnal.gov//store/mc/RunIISpring15DR74/GluGluToBulkGravitonToHHTo2B2G_M-260_narrow_13TeV-madgraph/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/5850C23B-6A12-E511-A4BD-0025905C96A6.root"

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
    muonTag    = cms.InputTag('flashggSelectedMuons'),
    electronTag= cms.InputTag('flashggSelectedElectrons'),
    photonTag  = cms.InputTag('flashggPhotons'),
    jetTag      = cms.InputTag('flashggFinalJets'),
    genTag     = cms.InputTag('flashggPrunedGenParticles'),
    lumiWeight = cms.double(1.),
)

from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents",200)
customize.setDefault("targetLumi",10e+3)
customize.crossSections.append("$CMSSW_BASE/src/APZ/fgg-ana/cross_sections.json")

customize(process)
