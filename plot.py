#! /usr/bin/env python
from optparse import OptionParser
import sys,os,datetime
from array import *
from ROOT import *
sys.path.append("utils")
import utils as u
import makeHTML as ht
import numpy as np
gROOT.SetBatch()

parser = OptionParser(usage="usage: %prog plots-version-path [options -c, -e, -p, -m]")
parser.add_option("-c","--cut",   dest="cut", type="int", default=5, help="Plots after a certain cut")
parser.add_option("--mass", dest="mass", type="int", default=125,    help="Signal sample (mass)")
parser.add_option("-m","--merge", dest="merge",action="store_true", default=False, help="Do merging?")
parser.add_option("-v","--verbose",dest="verbose",action="store_true", default=False, help="Verbose level")

parser.add_option("-p", "--period",dest="period", default="2012",  help="Year period; 2011 or 2012")
parser.add_option("--bkg",  dest="bkg",  action="store_true", default=False, help="Make plots from bkg sample")
parser.add_option("--qcd",  dest="qcd",  action="store_true", default=False, help="Include QCD samples")
parser.add_option("--sig",  dest="sig",  action="store_true", default=False, help="Signal MC")
parser.add_option("--data", dest="data", action="store_true", default=False, help="Data only")
parser.add_option("--vbf",  dest="vbf",  action="store_true", default=False, help="Use signal samples: ggH, vbf, vH")

parser.add_option("--spec",dest="spec",action="store_true", default=False, help="Make some special plots at the end.")
parser.add_option("--evt",dest="evt",action="store_true", default=False, help="Show raw events (not scaled to lumi) for MC samles")
parser.add_option("-e","--extra",dest="extra",action="store_true", default=False, help="Make all extra plots")
parser.add_option("--fit",dest="fit",action="store_true", default=False, help="Do the various fits")
parser.add_option("--apz",dest="apz",action="store_true", default=False, help="Discover new particle (requires apzTree)")
parser.add_option("--zjp",dest="zjp",action="store_true", default=False, help="Study J/Psi region and more (requires apzTree)")
parser.add_option("--hjp",dest="hjp",action="store_true", default=False, help="Study J/Psi region and more (requires apzTree)")
parser.add_option("-s", '--sel', dest="sel", type="string", default='mugamma',
                  help="Selection to be used. Options are: '4mu','2e2mu', 'zee','mugamma', 'elgamma'")

(opt, args) = parser.parse_args()

mass = opt.mass
sel  = opt.sel

comments = ["Bla",
            "Bla bla"]

HHresMass = ['300','650']
#HHresMass = ['250','300','340','350','400','450','500','600','700','800','900']

doLumiScale = 1
if opt.evt: doLumiScale=0

if __name__ == "__main__":
  timer = TStopwatch()
  timer.Start()

  if len(args) < 1:
    parser.print_usage()
    exit(1)

  ver    = sys.argv[1].rstrip('/')
  if 'vv/' in ver: ver = ver[3:].rstrip('/')
  cut=str(opt.cut)
  doMerge = opt.merge
  period  = opt.period
  doBkg   = opt.bkg

  u.setVerboseLevel(opt.verbose)
  #gROOT.ProcessLine(".L ../tdrstyle.C")
  #setTDRStyle()
  TH1.SetDefaultSumw2(kTRUE)

  import socket
  hostname = socket.gethostname()

  if '/tthome' in os.getcwd():
    pathBase = "/tthome/andrey/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath    = "/tthome/andrey/batch_output/zgamma/8TeV/"+ver
  elif 'cmslpc' in hostname:
    pathBase = "/uscms_data/d2/andreypz/html/zgamma/dalitz/"+ver+"_cut"+cut
    hPath    = "/eos/uscms/store/user/andreypz/batch_output/zgamma/8TeV/"+ver
  elif 'pcncu' in hostname or 'lxplus' in hostname:
    pathBase = "/afs/cern.ch/user/a/andrey/work/html/"+ver+"_cut"+cut
    hPath    = os.getcwd()+"/"+ver


  u.createDir(pathBase)

  dataFile  = None
  bkgFiles  = []
  bkgNames  = []

  u.setSelection(sel)
  subsel=sel


  qcdSamples = None
  if doBkg:
    #bkgFiles.append(TFile(hPath+"/m_ZG_"+subsel+"_"+period+".root","OPEN"))
    #bkgNames.append('ZG')
    #bkgFiles.append(TFile(hPath+"/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root","OPEN"))
    #bkgNames.append('DYJets50')

    bkgFiles.append(TFile(hPath+"/output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root","OPEN"))
    bkgNames.append('DiPhotons')

    bkgFiles.append(TFile(hPath+"/output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root"))
    bkgNames.append('GJets20')
    
    bkgFiles.append(TFile(hPath+"/output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root","OPEN"))
    bkgNames.append('GJets40')

    bkgFiles.append(TFile(hPath+"/output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root"))
    bkgNames.append('QCD_EM_Pt_40')
    
    bkgFiles.append(TFile(hPath+"/output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root"))
    bkgNames.append('QCD_EM_Pt_30to40')
    
    
    yields_bkg = []
    for b in xrange(len(bkgNames)):
      yields_bkg.append(u.getYields(bkgFiles[b],bkgNames[b],True))

    bkgZip = zip(bkgNames, bkgFiles)


    print '\t\t my BKG yields'
    print yields_bkg

  else: bkgZip = None

  sigGrav  = {}
  sigRadi  = {}
  sigGG  = {}
  sigVBF = {}
  sigVH  = {}
  for m in HHresMass:
    sigGrav[m]  = TFile(hPath+"/output_GluGluToBulkGravitonToHHTo2B2G_M-"+m+"_narrow_13TeV-madgraph.root", "OPEN")
    sigRadi[m]  = TFile(hPath+"/output_GluGluToRadionToHHTo2B2G_M-"+m+"_narrow_13TeV-madgraph.root", "OPEN")
    #sigGG[m]  = TFile(hPath+"/hhhh_ggH-mad"+m+"_1.root", "OPEN")
    #sigVBF[m] = TFile(hPath+"/hhhh_vbfH-mad"+m+"_1.root", "OPEN")
    #sigVH[m]  = TFile(hPath+"/hhhh_vH-mad"+ m+"_1.root", "OPEN")

  try:
    sigFileGG   = sigGG['125']
    sigFileVBF  = sigVBF['125']
    sigFileVH   = sigVH['125']
  except KeyError:
    print 'Key error in single H?  Who cares...'
    sigFileGG = sigFileVBF = sigFileVH = None

  if opt.zjp: sigFileZjp  = TFile(hPath+"/hhhh_ZtoJPsiGamma_1.root",     "OPEN")
  if opt.zjp: sigFileHjp  = TFile(hPath+"/hhhh_HiggsToJPsiGamma_1.root", "OPEN")
  #sigFileHjp  = TFile(hPath+"/hhhh_HiggsToJPsiGamma_RECO.root", "OPEN")

  #if opt.hzg: ggHZGFile   = TFile(hPath+"/hhhh_ggHZG-"+str(mass)+"_1.root", "OPEN")

  if opt.hjp: sigFile = sigFileHjp
  elif opt.zjp: sigFile = sigFileZjp
  else:         sigFile = sigFileGG

  if opt.data:
    if sel=='mu':
      dataFile = TFile(hPath+"/output_DoubleMuon.root","OPEN")
    elif sel=='el':
      dataFile = TFile(hPath+"/output_DoubleEG.root","OPEN")
    elif sel=='hhbbgg':
      dataFile = TFile(hPath+"/output_DoubleEG.root","OPEN")


  if opt.data:
    yields_data = u.getYields(dataFile)
    print '\t\t my data yields'
    print yields_data

  if opt.sig:
    if opt.zjp: yields_zjp  = u.getYields(sigFileZjp, 'ZtoJPsiGamma',  doLumiScale)
    if opt.hjp: yields_hjp  = u.getYields(sigFileHjp, 'HtoJPsiGamma',  doLumiScale)
    if sel in ['mugamma','elgamma']:
      yields_ggH  = u.getYields(sigFileGG,  'ggH-125',  doLumiScale)
      yields_vbf  = u.getYields(sigFileVBF, 'vbfH-125', doLumiScale)
      yields_vh   = u.getYields(sigFileVH,  'vH-125',   doLumiScale)
      yields_sig  = [sum(x) for x in zip(yields_ggH,yields_vbf,yields_vh)]

      if opt.verbose: print 'ggH yi', yields_ggH

    if sel=='hhbbgg':
      yields_grav = {}
      yields_radi = {}
      for m in HHresMass:
        yields_grav[m]   = u.getYields(sigGrav[m],'grav',doLumiScale)
        yields_radi[m]   = u.getYields(sigRadi[m],'radi',doLumiScale)

    # print 'Total sig yi', yields_sig

    # subdir = 'DY'
    # path = pathBase+"/"+subdir
    # if doBkg:
    #  path = pathBase+"/bkg_"+subdir
    #  u.createDir(path)
  #path = pathBase+"/"+subdir


  sigName = '#splitline{XX x Signal}{m_{H}=125 GeV}'
  if opt.zjp: sigName= '#splitline{XX x Signal}{Z #rightarrow J/Psi #gamma}'
  if opt.hjp: sigName= 'H #rightarrow J/Psi #gamma'


  sigZip = zip(['Grav 300','Radi 300','Grav 650', 'Radi 650'],
               [sigGrav['300'], sigRadi['300'], sigGrav['650'], sigRadi['650']])
  #sigZip = zip(['Grav 250','Grav 340','Grav 600'],
  #             [sigGrav['250'], sigGrav['340'], sigGrav['600']])


  if opt.bkg:
    print

    # u.drawAllInFile(dataFile, "Data", bkgZip, None, sigName, "DY", pathBase+'/DY', cut, "lumi")
    u.drawAllInFile(dataFile, "Data", bkgZip, None, sigName, "Main", pathBase+'/Main', cut, "lumi")

  else:
    print
    if opt.data and opt.sig:
      u.drawAllInFile(dataFile, "Data", None, sigZip, sigName, "Main", pathBase+'/DY', cut)
    if opt.data:
      u.drawAllInFile(dataFile, "Data", None, None, '', "Main", pathBase+'/Main-Data', cut)
    if opt.sig:
      u.drawAllInFile(None, None, None, sigZip, sigName,"Main", pathBase+'/Main-Sig', cut, 'norm')


  if opt.extra:
    # For the plots where cut number is set:
    if sel=='hhbbgg': dirList = ['Angles-bb','Angles-gg']
    else: dirList = ['Angles']

    for n in dirList:
      if opt.data:
        u.drawAllInFile(dataFile, "Data", bkgZip, None,"signal", n, pathBase+"/"+n, cut, "norm")
      else:
        u.drawAllInFile(None, None, None, sigZip, sigName, n, pathBase+"/"+n, cut, "norm")

    # For the cases without cut number:
    if sel=='hhbbgg': dirList = ['N','DiPho-Lead','DiPho-Sub', 'GEN', 'bJet-Lead','bJet-Sub']
    else: dirList = ['Photon','Ele','N', 'GEN']
    for n in dirList:
      if opt.data:
        u.drawAllInFile(dataFile, "Data", bkgZip, None,"signal", n, pathBase+"/"+n, None, "norm")




  if opt.spec:
    print 'Doing special plotting'
    u.createDir(pathBase+'/Spec')
    c1 = TCanvas("c4","small canvas",600,600);
    if sel=='hhbbgg':
      # DO the acc vs mass
      accept_grav = []
      accept_radi = []
      for m in HHresMass:
        print m,'Grav yields:\n', yields_grav[m]
        print m,'Radi yields:\n', yields_radi[m]
        accept_grav.append(yields_grav[m][8]/yields_grav[m][3])
        accept_radi.append(yields_radi[m][8]/yields_radi[m][3])

      massArray = np.array(HHresMass, dtype=float)
      print '\t DBG: mass array', massArray
      nPoints = len(massArray)
      accept_array = np.array(accept_grav)
      grAccGrav = TGraph(nPoints, massArray, accept_array)
      print '\t DBG: grav graph:'
      grAccGrav.Print('all')

      accept_array = np.array(accept_radi)
      grAccRadi = TGraph(nPoints, massArray, accept_array)
      print '\t DBG: Radi graph:'
      grAccRadi.Print('all')

      grAccGrav.Draw('AP')
      grAccRadi.Draw('same P')
      grAccGrav.SetMarkerStyle(20)
      grAccGrav.SetMarkerColor(kBlue)
      grAccRadi.SetMarkerStyle(22)
      grAccRadi.SetMarkerColor(kRed)

      grAccGrav.SetMinimum(0)
      grAccGrav.SetMaximum(1)
      leg = TLegend(0.55,0.75,0.8,0.89)
      leg.SetFillColor(kWhite)
      leg.SetBorderSize(0)
      leg.AddEntry(grAccGrav,'Graviton', 'p')
      leg.AddEntry(grAccRadi,'Radion', 'p')
      leg.Draw()
      grAccGrav.SetTitle(';m_{X} (GeV); Acceptance')
      c1.SaveAs(pathBase+'/Spec/acceptance_vs_mX.png')


    print 'End of special plotting'

  plot_types =[]
  dirlist = os.listdir(pathBase)
  for d in dirlist:
    if os.path.isdir(pathBase+"/"+d):
      plot_types.append(d)


  if doBkg:
    if sel=='hhbbgg':
      names = ['DiPhoton', 'GJets20','GJets40',"QCD_EM_Pt_40"]
      yields_all  = [yields_bkg[0],yields_bkg[1],yields_bkg[2],yields_bkg[3]]
      if opt.sig:
        names.append('Grav m = 300')
        yields_all.append(yields_grav['300'])
    else:
      names = ['DYJet50', 'Signal @125']
      yields_all = [yields_bkg[0],yields_sig]
  else:
    if opt.zjp or opt.hjp:
      names = ['H to J/Psi &gamma;','Z to J/Psi &gamma;', 'ggH-125']
      yields_all = [yields_hjp, yields_zjp, yields_ggH]
    #elif not opt.apz:
    elif opt.sig:
      if sel in ['mugamma','elgamma']:
        names = ['Sig: total','ggH','vbfH','VH']
        yields_all  = [yields_sig, yields_ggH,yields_vbf, yields_vh]
      elif sel=='hhbbgg':
        names = ['Grav M=%s'%m for m in HHresMass]
        yields_all  = [yields_grav[m] for m in HHresMass]
    else:
      yields_all = None

  if opt.data:
    names[:0] = ['Data']
    yields_all[:0] = [yields_data]

  table_all  = u.yieldsTable(yields_all, names)


  if opt.hjp: precision='%.4f'
  elif sel=='hhbbgg': precision='%.0f'
  else:       precision='%.2f'

  u.makeTable(table_all,"all", "html", precision)
  u.makeTable(table_all,"all", "html", precision)
  u.makeTable(table_all,"all", "twiki", precision)
  u.makeTable(table_all,"all", "tex")

  os.system("cat yields_all.html > yields.html")

  if sel=='hhbbgg':
    defaultPage = 'Main'
    ht.makeHTML("X &rarr; HH &rarr; bb &gamma;&gamma; decay plots",pathBase, plot_types, comments, defaultPage)
  else:
    defaultPage = 'Main'
    ht.makeHTML("X &rarr; HH &rarr; bb &gamma;&gamma; decay plots",pathBase, plot_types, comments, defaultPage)

  print "\n\t\t finita la comedia \n"
