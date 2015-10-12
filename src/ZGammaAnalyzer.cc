#include "../interface/ZGammaAnalyzer.h"

typedef std::pair<std::string,float> IdPair;
//typedef std::pair<std::string,Bool_t> IdPair;

ZGammaAnalyzer::ZGammaAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs):
  DummyAnalyzer::DummyAnalyzer(cfg, fs)
{
  cout<<"\t ZZZZZZZZ \t Constructructor in "<<__PRETTY_FUNCTION__<<endl;
  FHM = new FggHistMakerZgamma(hists);
}

void ZGammaAnalyzer::beginJob()
{
  //edm::LogInfo("My info");
  cout<<"\t ZZZZZZZZ \tBegin job in: "<<__PRETTY_FUNCTION__<<endl;

  cout<<"lumiWei="<<lumiWeight_<<endl;
  cout<<"lep="<<lep_<<endl;
}


void ZGammaAnalyzer::analyze(const edm::EventBase& event)
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
  //else
  //throw cms::Exception("DYJets", "This has to be an MC sample!")<<ww<<"  W0="<<W0;
  //if (fabs(ww)!=W0) throw cms::Exception("ETALONS!","The weights are not the same. Can not deal with this... ww=")<<ww<<"  W0="<<W0;
  //cout<<" evt wei = "<<ww<<endl;

  if (ww<0) count_neg++;
  else count_pos++;

  totEvents++;
  totWeights+=ww;
  
  CountEvents(0, "Ntuple events",ww,fcuts);
  FillHistoCounts(0, ww);

  FHM->Reset(1, 1);


  edm::Handle<vector<reco::GenParticle> > genParts;
  event.getByLabel( myGen_, genParts );

  if (!isRealData){
    for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {

      hists->fill1DHist(igen->pdgId(),"pdg_id","All PDGs", 40,-20,20, ww, "PDG");
    }
  }


  vector<TLorentzVector> myLeptons, myPhotons; 
  //vector<TLorentzVector> myEle, myMuo; 

  // Handle to the muon collection
  edm::Handle<std::vector<flashgg::Muon> > muons;
  event.getByLabel(muons_, muons);

  edm::Handle<std::vector<flashgg::Electron> > electrons;
  event.getByLabel(electrons_, electrons);


  
  // ***
  // LEPTONS (Mu or El)
  // ***

  if (lep_=="mu"){
    
    for(std::vector<flashgg::Muon>::const_iterator it=muons->begin(); it!=muons->end(); ++it){
      TLorentzVector tmp = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());
      if (it->pt()>20 && it->isLooseMuon()){
	myLeptons.push_back(tmp);
	FHM->MakeMuonPlots(*it);
	//cout<<"Muon "<<mu1->pt()<<endl;
      }
    }
  }
  
  else if (lep_=="el"){
    for(std::vector<flashgg::Electron>::const_iterator it=electrons->begin(); it!=electrons->end(); ++it){
      
      hists->fill1DHist(it->pt(),"ele_Pt",";Electrons p_{T}, GeV", 40,-20,20, ww, "ELE");
      
      //const std::vector<IdPair> IDs = it->electronIDs();
      //for (size_t t=0; t<IDs.size(); t++)
      //cout<<t<<"  ID name:"<<IDs[t].first<<"  ID pass?="<<IDs[t].second<<endl;
  
      TLorentzVector tmp = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());
      if (it->pt()>20 && it->electronID("cutBasedElectronID-CSA14-50ns-V1-standalone-loose")==1){
	//if (it->pt()>20 && it->electronID("eidRobustLoose")>=6){
	//Try HEEP id
	myLeptons.push_back(tmp);
	FHM->MakeElectronPlots(*it);
      }
    }

    //for (size_t s=0; s<myEle.size();s++){
    //cout<<s<<" pt="<<myEle[s].Pt()<<endl;
    //}

  }
  

  sort(myLeptons.begin(), myLeptons.end(), P4SortCondition);

  
  if (myLeptons.size()<2) return;
  CountEvents(1, "Two leptons with pT > 20 GeV",ww,fcuts);
  FillHistoCounts(1, ww);


  TLorentzVector l1 = myLeptons[0];
  TLorentzVector l2 = myLeptons[1];

  FHM->SetLeptons(l1,l2);
  
  Float_t Mll = (l1+l2).M();


  if (l1.Pt() < 35 || l2.Pt() < 35 )   return;
  CountEvents(2, "pT>35",ww,fcuts);
  FillHistoCounts(2, ww);
  FHM->MakeMainHistos(2, ww);

  hists->fill1DHist(Mll,"mll_cut2",";M(ll), GeV", 40, 60,120, ww, "DY");


  // ***
  // Z peak
  // ***

  if (Mll < 60 || Mll > 120) return;

  CountEvents(3, "60 < M(ll) < 120",ww,fcuts);
  FillHistoCounts(3, ww);
  FHM->MakeMainHistos(3, ww);

  hists->fill1DHist(Mll,"mll_cut3",";M(ll), GeV", 40, 60,120, ww, "DY");
  hists->fill1DHist(l1.DeltaR(l2),"dR_l1l2_cut3",";dR(ll)", 40, 0,6, ww, "DY");


  // ***
  // PHOTON
  // ***

  edm::Handle<std::vector<flashgg::Photon> > photons;
  event.getByLabel(photons_, photons);

  for(std::vector<flashgg::Photon>::const_iterator it=photons->begin(); it!=photons->end(); ++it){
    TLorentzVector tmp = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());

    Bool_t isMatchedToLepton = false;
    for(std::vector<TLorentzVector>::const_iterator l=myLeptons.begin(); l!=myLeptons.end(); ++l)
      if (tmp.DeltaR(*l) < 0.4) { isMatchedToLepton = true; break;};
    if (isMatchedToLepton) continue;

    if (it->pt()>15 && it->photonID("PhotonCutBasedIDLoose")){
      myPhotons.push_back(tmp);
      FHM->MakePhotonPlots(*it);
    }
    //const std::vector<IdPair> IDs = it->photonIDs();
    //for (size_t t=0; t<IDs.size(); t++)
    //cout<<t<<"  ID name:"<<IDs[t].first<<"  ID pass?="<<IDs[t].second<<endl;

  }

  sort(myPhotons.begin(), myPhotons.end(), P4SortCondition);

  //cout<<"Size of Leptons collection =  "<<myLeptons.size()<<endl;
  //cout<<"Size of Photons collection =  "<<myPhotons.size()<<endl;
 
  //cout<<"Size of Muo collection =  "<<myMuo.size()<<endl;
  //cout<<"Size of Ele collection = "<<myEle.size()<<endl;


  if (myPhotons.size()<1) return;
  TLorentzVector gamma = myPhotons[0];

  FHM->SetGamma(gamma);
  Float_t Mllg = (l1+l2+gamma).M();
  
  CountEvents(4, "Photon with pT > 15",ww,fcuts);
  FillHistoCounts(4, ww);
  FHM->MakeMainHistos(4, ww);

  hists->fill1DHist(Mll, "mll_cut4", ";M(ll), GeV", 40, 60,120, ww, "DY");
  hists->fill1DHist(Mllg,"mllg_cut4",";M(llg),GeV", 40,  60,400, ww, "DY");
  hists->fill1DHist(l1.DeltaR(l2),"dR_l1l2_cut4",";dR(ll)", 40, 0,6, ww, "DY");
  hists->fill1DHist(l1.DeltaR(gamma),"dR_l1g_cut4",";dR(l1g)", 40, 0,6, ww, "DY");
  hists->fill1DHist(l2.DeltaR(gamma),"dR_l2g_cut4",";dR(l2g)", 40, 0,6, ww, "DY");
  

  CountEvents(5, "",ww,fcuts);
  FillHistoCounts(5, ww);

}

void ZGammaAnalyzer::endJob()
{
  cout<<"\t ZZZZZZZZ \t END JOB in "<<__PRETTY_FUNCTION__<<endl;
  DummyAnalyzer::endJob();
}
