#include "../interface/GenAnalyzer.h"

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
    if (totEvents==0)  W0 = fabs(ww); // This is the Etalon!

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


  std::vector<TLorentzVector> gen_photons;
  TLorentzVector gan_gamma1, gen_gamma2;

  if (!isRealData){

    //Bool_t yes = 0;
    TLorentzVector tmp;
    
    /*
    for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {
      
      hists->fill1DHist(igen->pdgId(),"pdg_id","All PDGs", 40,-20,20, ww, "PDG");

      if (igen->pdgId()==25){
	//std::cout<<"\t\t 777 Found the HIGGS! 777"<<std::endl;
	//std::cout<<"\t Its *mother* is: "<<igen->mother()->pdgId()<<std::endl;
	//std::cout<<"\t Its *daughters* are: "<<igen->daughter(0)->pdgId()<<" and "<<igen->daughter(1)->pdgId()
	//	 <<"\n \t and their statuses are "<<igen->daughter(0)->status()<<" and "<<igen->daughter(1)->status()<<std::endl;


	//if (igen->daughter(0)->pdgId()==22 && igen->daughter(0)->status()!=1)
	//{
	//  std::cout<<"The photon from the Higgs is not status 1"<<std::endl;
	//yes = 1;
	//}
      }
    }
    */

    
    for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {
      
      //if (igen->pdgId()==22 && igen->status()==1 && igen->mother()->pdgId()==25){
      //if (igen->pdgId()==22 && igen->status()==1 && igen->isHardProcess()){
      //if (igen->pdgId()==22 && igen->status()==1 && igen->fromHardProcess()){
      if (igen->pdgId()==22 && igen->status()==1 && 
	  (igen->fromHardProcessFinalState() || igen->mother()->pdgId()==25)){
	tmp.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());
	gen_photons.push_back(tmp);
      }

      //if (yes && igen->pdgId()==22 && igen->status()==1){
      //std::cout<<totEvents<<"\t\t 777 Found a Photon! 777  status="<<igen->status()
      //	 <<"\n IsHard Process = "<<igen->isHardProcess()
      //	 <<" fromHardProcessFinalState = "<<igen->fromHardProcessFinalState()<<std::endl;
      //std::cout<<"\t Its Mother id "<< igen->mother()->pdgId()
      //	 <<"\t Its grandMother id "<< igen->mother()->mother()->pdgId()<<std::endl;
      // //<<"\n uniqueMother ID = "<<igen->uniqueMother().pdgId()<<std::endl;
      //}
    }
  }
  else throw cms::Exception("Dude! This is Data. You can't run gen-level analysis on Data...");

  

  if (gen_photons.size()<2) return;
  CountEvents(1, "Two gen photons from HardProcess",ww,fcuts);
  FillHistoCounts(1, ww);
      
  sort(gen_photons.begin(), gen_photons.end(), P4SortCondition);



  
}

void GenAnalyzer::endJob()
{
  cout<<"\t GGGGGGGGG \t END JOB in "<<__PRETTY_FUNCTION__<<endl;
  DummyAnalyzer::endJob();
}
