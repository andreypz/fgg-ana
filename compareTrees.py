#! /usr/bin/env python
import sys,argparse
from ROOT import *
gROOT.SetBatch()

parser =  argparse.ArgumentParser(description='Ploting my plots', usage="./compare.py --f1 fName1 --f2 fName2")
parser.add_argument("--f1",  dest="f1", type=str, required=True, help="Filename with flat tree")
parser.add_argument("--f2",  dest="f2", type=str, required=True, help="Filename with flat tree")

opt = parser.parse_args()


f1 = TFile(opt.f1,'OPEN')
f2 = TFile(opt.f2,'OPEN')
type1 = 0
type2 = 0

if f1.Get("fsDir/bbggSelectionTree"):
    t1 = f1.Get('fsDir/bbggSelectionTree')
    type1=1
elif f1.Get("bbggSelectionTree"):
    t1 = f1.Get('bbggSelectionTree')
    type1=2
else:
    print 'File', opt.f1, 'does not contain expected trees with data to compare.'
    sys.exit(0)
    
if f2.Get("fsDir/bbggSelectionTree"):
    t2 = f2.Get('fsDir/bbggSelectionTree')
    type2 = 1
elif f2.Get("bbggSelectionTree"):
    t2 = f2.Get('bbggSelectionTree')
    type2=2
else:
    print 'File', opt.f2, 'does not contain expected trees with data to compare.'
    sys.exit(0)
    
ev1 = set()
ev1More = []
for e in t1:
    mgg = (e.leadingPhoton+e.subleadingPhoton).M()
    mjj = (e.leadingJet+e.subleadingJet).M()
    if type1==2 and not e.isSignal: continue
    ev1.add((int(e.run), int(e.event)))
    #ev1More.append((int(e.run), int(e.event), mgg, mjj, e.leadingPhoton.Pt(), e.subleadingPhoton.Pt()))
    ev1More.append({'run':int(e.run), 'evt':int(e.event), 'mgg':mgg, 'mjj':mjj, 'ptg1':e.leadingPhoton.Pt(), 'ptg2':e.subleadingPhoton.Pt(),
                    'j1pt': e.leadingJet.Pt(), 'j2pt': e.subleadingJet.Pt(),
                    'j1-bDis': e.leadingJet_bDis, 'j2-bDis': e.subleadingJet_bDis,})

ev2 = set()
ev2More = []
for e in t2:
    mgg = (e.leadingPhoton+e.subleadingPhoton).M()
    mjj = (e.leadingJet+e.subleadingJet).M()
    if type2==2 and not e.isSignal: continue
    ev2.add((int(e.run), int(e.event)))
    ev2More.append({'run':int(e.run), 'evt':int(e.event), 'mgg':mgg, 'mjj':mjj, 'ptg1':e.leadingPhoton.Pt(), 'ptg2':e.subleadingPhoton.Pt()})

print 'Total events in file %s: %i' %(opt.f1, len(ev1))
print 'Total events in file %s: %i' %(opt.f2, len(ev2))
print 'Intersection of the two: %i' %(len(ev1.intersection(ev2)))

print 'Events in %s but not in %s: %i' %(opt.f1, opt.f2, len(ev1.difference(ev2)))
print '\t And here is a first few of those:'
for i,a in enumerate(ev1.difference(ev2)):
    print i, a
    for e in ev1More:
        if e['run']==a[0] and e['evt']==a[1]: print e        
    if i>5:
        break

print 'Events in %s but not in %s: %i' %(opt.f2, opt.f1, len(ev2.difference(ev1)))
print '\t And here is a first few of those:'
for i,a in enumerate(ev2.difference(ev1)):
    print i, a
    for e in ev2More:
        if e['run']==a[0] and e['evt']==a[1]: print e
    if i>5:
        break
