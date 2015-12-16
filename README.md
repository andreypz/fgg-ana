# Analyzers to run over flash-gg microAOD ntuples
To begin with this analyzers you must to set up the *flashgg*  framework first, see here: https://github.com/cms-analysis/flashgg
Then we need *bbgTools*, since some of the methods from there are used: https://github.com/ResonantHbbHgg/bbggTools

Finally, to get this code:
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



## How to create your own analyzer:
 * Write _src/MyAna.cc_ and _interface/MyAna.h_ files. It is a good idea to inherit it from DummyAnalyzer (see my other analyzers: GenAna, ZgammaAna), but this is not necessary.
 * Make your executable file at _bin/myAnaRun.cc_, and include it in _bin/Buildfile.xml_
 * Create your config file, _myana.py_, and run it with fgg comand like this:
```
 fggRunJobs.py --load jobs-myAna.json -d testDir -x myAnaRun myana.py  maxEvents=100
```

If you want to run a job without fgg scripts, first comment out the ```customize(process)``` string at the end of your myana.py file. Then run it simply with ```myAnaRun myana.py```.   

