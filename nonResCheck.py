#! /usr/bin/env python
import argparse
from ROOT import *
gROOT.SetBatch()
TH1.SetDefaultSumw2(kTRUE)

parser =  argparse.ArgumentParser(description='Ploting my plots', usage="./nonResCheck.py")
#parser.add_argument("-f",  dest="fName", type=str, default='fName1.root', help="Filename with flat tree")

opt = parser.parse_args()


chain = TChain('fsDir/GenTree')
figDir = 'FigNodes'
dirName = 'NonResOut2/'
f = {}
h_mHH = {}

nodes = ['2','3','4','5','6','7','8','9','10','11','12','13','SM','box']
#nodes = ['box','2']

c1 = TCanvas("c1","small canvas",600,600)
c2 = TCanvas("c2","big canvas",  600,700);
c1.cd()
for n in nodes:
    fName = dirName+'/output_GluGluToHHTo2B2G_node_'+n+'_13TeV-madgraph.root'

    if n!='SM' and n!='box':
        chain.Add(fName)

    f[n] = TFile(fName,'OPEN')
    t=f[n].Get('fsDir/GenTree')
    hname = 'mHH_node_'+n
    t.Draw('mHH>>'+hname+'(400, 0,1800)','')
    h_mHH[n] = gDirectory.Get(hname).Clone()

    c1.SaveAs(figDir+'/node_'+n+'.png')


## Here we will read the parameter grid
#  and find out which benchmark node correspond to what point in the grid

transFile = '../generateHH/temperature/list_all_translation_1507.txt'
grid = {}
wNumToNode = {}

import csv
with open(transFile, 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for i,row in enumerate(spamreader):
        # print i, row, row[0]
        grid[str(i)] = {'kl':float(row[0]), 'kt':float(row[1]), 'c2':float(row[2]), 'cg':float(row[3]), 'c2g':float(row[4])}
        k = str(i)
        #print i,grid[k]

        # The numbers below are taken from the Table 1 in our AN-16-030:

        #if grid[k]['kl']==1.0 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
        #    wNumToNode['SM'] = i
        #if grid[k]['kl']==7.5 and grid[k]['kt']==2.5 and grid[k]['c2']==-0.5 and grid[k]['cg']==0 and grid[k]['c2g']==0:
        #    wNumToNode['2'] = i
        if grid[k]['kl']==0.0001 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
            wNumToNode['box'] = i

print wNumToNode


## Here is the moment of truth: will the re-weighted histogramms match the one of each node before reweighting?

gROOT.ProcessLine(".L ~/tdrstyle.C")
setTDRStyle()
#gROOT.ForceStyle()

gStyle.SetOptStat(0)
c2.cd()
pad1 = TPad("pad1","pad1",0,0.3,1,1);
pad2 = TPad("pad2","pad2",0,0,1,0.3);
pad1.SetBottomMargin(0)
pad1.Draw()
pad2.SetBottomMargin(0.25)
pad2.SetTopMargin(0)
pad2.Draw()

for n,wN in wNumToNode.iteritems():
    print n, wN
    if n not in nodes: continue

    name1 = 'mHH_wN_'+n
    pad1.cd()
    chain.Draw('mHH>>'+name1+'(400,0,1800)','NRWeights['+str(wN)+']','')
    # gDirectory.Print()

    h = gDirectory.Get(name1)

    h1_n = h.DrawNormalized('hist')
    h2_n = h_mHH[n].DrawNormalized('same',1.0)


    h1_n.SetLineColor(kRed+2)
    h1_n.SetMaximum(0.03)
    h1_n.SetTitle('; m(HH), GeV')

    c2.cd()
    pad2.cd()
    r = h1_n.Clone()
    r.Divide(h2_n)
    r.GetYaxis().SetTitle("Ratio")
    r.SetMaximum(1.4)
    r.SetMinimum(0.6)
    r.GetYaxis().SetNdivisions(206)
    r.GetYaxis().SetTitleOffset(0.4)
    r.SetTitleSize(0.1,"XYZ")
    r.SetLabelSize(0.1,"XY")
    r.SetLineColor(kBlack)
    r.Draw("e1p")
    gPad.RedrawAxis()

    c2.SaveAs(figDir+'/weighted_to_node_'+n+'.png')



c1.cd()

HowMany = 150
chain.Draw('mHH>>mHH_0(400,0,1800)','','')
N0 = gDirectory.Get('mHH_0').Integral()
#N0 = chain.GetEntries()

print 'Number of entries in the chain =', N0
hNEvt = TH1F('hNEvt','Normalization; (N_{0}-N)/N_{0};Number of Points', 50,-0.02,0.02)
for w in xrange(0,HowMany):
    hname = 'h_w'+str(w)
    chain.Draw('mHH>>'+hname+'(400,0,1800)','NRWeights['+str(w)+']','')

    hInt = gDirectory.Get(hname).Integral()
    print 'weight index=',w,'Integral = ', hInt

    hNEvt.Fill((N0-hInt)/N0)

hNEvt.UseCurrentStyle()
gStyle.SetOptStat(1111)
hNEvt.Draw('hist')
c1.SaveAs(figDir+'/Normalization.png')

'''

import numpy as np
sumOfWeights = np.zeros(1507)

for e in chain:
    for n in xrange(0,5):
        sumOfWeights[n]+=e.NRWeights[n]


'''
