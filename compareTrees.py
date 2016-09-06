#! /usr/bin/env python
import argparse
from ROOT import *
gROOT.SetBatch()

parser =  argparse.ArgumentParser(description='Ploting my plots', usage="./compare.py --fAnd fName1 --fRaf fName2")
parser.add_argument("--fAnd",  dest="fAnd", type=str, default='fName1.root', help="Filename with flat tree")
parser.add_argument("--fRaf",  dest="fRaf", type=str, default='fName2.root', help="Filename with flat tree")

opt = parser.parse_args()



fAnd = TFile(opt.fAnd,'OPEN')
fRaf = TFile(opt.fRaf,'OPEN')


tAnd = fAnd.Get('fsDir/bbggSelectionTree')
tRaf = fRaf.Get('fsDir/bbggSelectionTree')

evAnd = set()
for e in tAnd:
    evAnd.add((str(e.run), str(e.event)))

evRaf = set()
for e in tRaf:
    evRaf.add((str(e.run), str(e.event)))

print 'Total events in file %s: %i' %(opt.fAnd, len(evAnd))
print 'Total events in file %s: %i' %(opt.fRaf, len(evRaf))
print 'Intersection of the two: %i' %(len(evAnd.intersection(evRaf)))
print 'Events is %s but not in %s: %i' %(opt.fAnd, opt.fRaf, len(evAnd.difference(evRaf)))
print '\t And here is a first few of those:'
for i,a in enumerate(evAnd.difference(evRaf)):
    print i, a
    if i>5:
        break

print 'Events is %s but not in %s: %i' %(opt.fRaf, opt.fAnd, len(evRaf.difference(evAnd)))
print '\t And here is a first few of those:'
for i,a in enumerate(evRaf.difference(evAnd)):
    print i, a
    if i>5:
        break
