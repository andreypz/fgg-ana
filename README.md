# Analyzers to run over flash-gg microAOD ntuples
To start you need to set up flashgg first, see here: 	https://github.com/cms-analysis/flashgg

Then, do:
```
 cd $CMSSW_BASE/src
 mkdir APZ
 cd APZ
 git clone git@github.com:andreypz/fgg-ana.git
 scram b -j9
```

Now, you have it compiled. Enjoy, by running:
```
 fggRunJobs.py --load jobs_dy-el.json -d testDir -x zgammaRun zgAna.py  maxEvents=100
```
where jobs_dy-el.json is the file prepared following instructions here: https://github.com/cms-analysis/flashgg/tree/master/MetaData
