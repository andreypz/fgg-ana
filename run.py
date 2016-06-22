#!/usr/bin/env python

#from ROOT import *
import sys,os

import argparse
parser =  argparse.ArgumentParser(description='Ploting my plots', usage="./run.py -d DIR [--sig or --few or --bkg]")
parser.add_argument("dir", help="The name of the directory where all root files with plots are stored.")

parser.add_argument("--sig", default='res', dest="sig", choices=['all', 'res', 'nres'], help="Run over all of the signal samples")
parser.add_argument("--few", action="store_true", default=False,  dest="few", help="Run over a few signal samples (for sync etc.)")
parser.add_argument("--bkg", action="store_true", default=False,  dest="bkg", help="Run over background samples")
# Make those above a chice options
parser.add_argument("--lsf", action="store_true", default=False,  dest="lsf", help="submitting to lsf. By default it's running in local mode.")
opt = parser.parse_args()
print ' \t Options specified:\n', opt

jobFile = "jobs-test.json"

if opt.sig=='all':
  jobFile = "jobs-sig-hh-all.json"
if opt.sig=='res':
  jobFile = "jobs-sig-hh-res.json"
if opt.sig=='nres':
  jobFile = "jobs-sig-hh-nonres.json"
if opt.few:
  jobFile = "jobs-hh-sig-few.json"
if opt.bkg:
  jobFile = "jobs-bkg-hh.json"

myArgs = ''
if opt.bkg:
  myArgs+=' -H -D -n 1'

suffix = '--no-use-tarball'
if opt.lsf:
  suffix = '-q 1nh'
  myArgs+=' -m 2'
else:
  myArgs+=' -m 0'

os.system('fggRunJobs.py --load data/'+jobFile+' -d '+opt.dir+' '+myArgs+' -x hhRun hhAna.py '+ suffix)
