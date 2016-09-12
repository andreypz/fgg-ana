#! /usr/bin/env python
import argparse
from ROOT import *
gROOT.SetBatch()

parser =  argparse.ArgumentParser(description='Ploting my plots', usage="./nonResCheck.py")
#parser.add_argument("-f",  dest="fName", type=str, default='fName1.root', help="Filename with flat tree")

opt = parser.parse_args()


chain = TChain('fsDir/GenTree')
figDir = 'FigNodes'
dirName = 'NonResOut/'
f = {}

nodes = ['2','3','4','5','6','7','8','9','10','11','12','13','SM','box']

c1 = TCanvas("c4","small canvas",600,600)
for n in nodes:
    fName = dirName+'/output_GluGluToHHTo2B2G_node_'+n+'_13TeV-madgraph.root'

    if n!='SM' and n!='box':
        chain.Add(fName)

    f[n] = TFile(fName,'OPEN')
    t=f[n].Get('fsDir/GenTree')
    t.Draw('mHH>>mHH(400, 0,2000)','')
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

        if grid[k]['kl']==1.0 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
            wNumToNode['SM'] = i
        if grid[k]['kl']==7.5 and grid[k]['kt']==1.0 and grid[k]['c2']==-1.0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
            wNumToNode['2'] = i
        if grid[k]['kl']==1.0 and grid[k]['kt']==1.0 and grid[k]['c2']==0.5 and grid[k]['cg']==-0.8 and grid[k]['c2g']==0.6:
            wNumToNode['3'] = i
        if grid[k]['kl']==1.0 and grid[k]['kt']==1.0 and grid[k]['c2']==-1.5 and grid[k]['cg']==0 and grid[k]['c2g']==-0.8:
            wNumToNode['4'] = i
        if grid[k]['kl']==-3.5 and grid[k]['kt']==1.5 and grid[k]['c2']==-3.0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
            wNumToNode['5'] = i
        if grid[k]['kl']==1.0 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==0.8 and grid[k]['c2g']==-1.0:
            wNumToNode['6'] = i
        if grid[k]['kl']==2.4 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==0.2 and grid[k]['c2g']==-0.2:
            wNumToNode['7'] = i
        if grid[k]['kl']==5.0 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==0.2 and grid[k]['c2g']==-0.2:
            wNumToNode['8'] = i
        if grid[k]['kl']==15.0 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==-1.0 and grid[k]['c2g']==1.0:
            wNumToNode['9'] = i
        if grid[k]['kl']==1.0 and grid[k]['kt']==1.0 and grid[k]['c2']==1.0 and grid[k]['cg']==-0.6 and grid[k]['c2g']==0.6:
            wNumToNode['10'] = i
        if grid[k]['kl']==10.0 and grid[k]['kt']==1.5 and grid[k]['c2']==-1.0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
            wNumToNode['11'] = i
        if grid[k]['kl']==2.4 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==1.0 and grid[k]['c2g']==-1.0:
            wNumToNode['12'] = i
        if grid[k]['kl']==15.0 and grid[k]['kt']==1.0 and grid[k]['c2']==1.0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
            wNumToNode['13'] = i
        if grid[k]['kl']==0.0001 and grid[k]['kt']==1.0 and grid[k]['c2']==0 and grid[k]['cg']==0 and grid[k]['c2g']==0:
            wNumToNode['box'] = i

print wNumToNode


## Here is the moment of truth: will the re-weighted histogramms match the one of each node before reweighting?

for n in ['box','SM']:
#for n in nodes:
    chain.Draw('mHH>>mHH(400, 0,2000)','NRWeights['+str(wNumToNode[n])+']')

    c1.SaveAs(figDir+'/weighted_to_node_'+n+'.png')
