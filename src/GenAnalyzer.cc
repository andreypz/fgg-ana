#include "../interface/GenAnalyzer.h"

//bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

GenAnalyzer::GenAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs):
  DummyAnalyzer::DummyAnalyzer(cfg, fs)
{
  cout<<"\t GGGGGGGGG \t Constructructor in "<<__PRETTY_FUNCTION__<<endl;
  FHM = new FggHistMakerBase(hists);
}

void GenAnalyzer::beginJob()
{
  //edm::LogInfo("My info");
  cout<<"\t GGGGGGGGG \t Begin job in: "<<__PRETTY_FUNCTION__<<endl;

  cout<<"lumiWei="<<lumiWeight_<<endl;
  cout<<"lep="<<lep_<<endl;
}


void GenAnalyzer::analyze(const edm::EventBase& event)
{
  Double_t ww = 1;
  isRealData = event.isRealData();
  if (!isRealData){
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    event.getByLabel(edm::InputTag("generator"), genEvtInfo );
    ww = genEvtInfo->weight();
    if (totEvents==0)  W0 = fabs(ww); // This is an Etalon!

    if (fabs(ww)!=W0) throw cms::Exception("ETALON","The weights are not the same. Can not deal with this... ww=")<<ww<<"  W0="<<W0;

    ww/=W0; //Divide by the etalon, get +/-1
  }

  if (ww<0) count_neg++;
  else count_pos++;

  totEvents++;
  totWeights+=ww;
  
  CountEvents(0, "Ntuple events",ww,fcuts);
  FillHistoCounts(0, ww);

  edm::Handle<vector<reco::GenParticle> > genParts;
  event.getByLabel( myGen_, genParts );

  if (!isRealData){
    for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {

      hists->fill1DHist(igen->pdgId(),"pdg_id","All PDGs", 40,-20,20, ww, "PDG");
    }
  }



  
}

void GenAnalyzer::endJob()
{
  cout<<"\t GGGGGGGGG \t END JOB in "<<__PRETTY_FUNCTION__<<endl;
  DummyAnalyzer::endJob();
}
