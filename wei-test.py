#!/usr/bin/env python
import ROOT

fName = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15-50ns-Spring15BetaV4-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150813_120650/0000/myMicroAODOutputFile_1.root"

pp = "root://eoscms//eos/cms"
fin = ROOT.TFile.Open("%s/%s" % (pp, fName))


totWei1 = ROOT.TH1D("totWei1","totWei1", 1, 0,2)
totWei2 = ROOT.TH1D("totWei2","totWei2", 1, 0,2)

lumi = fin.Get("LuminosityBlocks")
ne = lumi.Draw("1>>totEvents(1,0,2)","edmMergeableCounter_eventCount__FLASHggMicroAOD.obj.value","goff")
ne = lumi.Draw("1>>totWei1","edmMergeableDouble_weightsCount_totalWeight_FLASHggMicroAOD.obj.value","goff")

totEvents = ROOT.gDirectory.Get("totEvents").Integral()
totWeights = ROOT.gDirectory.Get("totWei1").Integral()

print  totEvents, totWeights

events = fin.Get("Events")
totEvents = events.GetEntries()
#ne = events.Draw("1>>totWeights(1,0,2)","GenEventInfoProduct_generator__GEN.obj.weights_[0]","goff")
ne = events.Draw("1>>totWei2","GenEventInfoProduct_generator__SIM.obj.weights_[0]","goff")
totWeights = ROOT.gDirectory.Get("totWei2").Integral()

print  totEvents, totWeights
