#include "../interface/HHbbggAnalyzer.h"

HHbbggAnalyzer::HHbbggAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs):
  DummyAnalyzer::DummyAnalyzer(cfg, fs),
  inputTagJets_(cfg.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
  phoIDcutEB_(cfg.getUntrackedParameter<std::vector<double > >("phoIDcutEB") ),
  phoIDcutEE_(cfg.getUntrackedParameter<std::vector<double > >("phoIDcutEE") )
{
  cout<<"\t HHHHHbbbbggggg \t Constructructor in "<<__PRETTY_FUNCTION__<<endl;
  FHM = new FggHistMakerHHbbgg(hists);
  tools_ = bbggTools();

}

void HHbbggAnalyzer::beginJob()
{
  cout<<"\t HHHHHHHHHH \t Begin job in: "<<__PRETTY_FUNCTION__<<endl;

  cout<<"lumiWei="<<lumiWeight_<<endl;
  cout<<"lep="<<lep_<<endl;
}


void HHbbggAnalyzer::analyze(const edm::EventBase& event)
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

  FHM->Reset(1, 1);


  /*
  edm::Handle<vector<reco::GenParticle> > genParts;
  event.getByLabel( myGen_, genParts );

  std::vector<TLorentzVector> gen_photons, gen_jets;
  TLorentzVector gen_gamma1, gen_gamma2;
  TLorentzVector gen_bjet1, gen_bjet2;
  TLorentzVector gen_bQ1, gen_bQ2;


  for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {

    if (abs(igen->pdgId())==5 && igen->mother()->pdgId()==25) {
      //std::cout<<totEvents<<"\t\t BBB Found b-quark from Higgs! BBB  Its status = "<<igen->status()<<std::endl;

      tmp.SetPxPyPzE(igen->px(), igen->py(), igen->pz(), igen->energy());
      if (igen->pdgId()==5)
	gen_bQ1 = tmp;
      else if (igen->pdgId()==-5)
	gen_bQ2 = tmp;
    }


    if (igen->pdgId()==22 && igen->isPromptFinalState()) {
      //if (igen->pdgId()==22 && igen->fromHardProcessFinalState()) {
      //if (igen->pdgId()==22 && igen->status()==1 &&
      //(igen->fromHardProcessFinalState() || igen->mother()->pdgId()==25 || igen->mother()->mother()->pdgId()==25)){
      tmp.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());
      gen_photons.push_back(tmp);

    }

  }

  */


  vector<TLorentzVector> myLeptons, myPhotons;

  edm::Handle<std::vector<flashgg::Photon> > photons;
  event.getByLabel(photons_, photons);

  for(std::vector<flashgg::Photon>::const_iterator it=photons->begin(); it!=photons->end(); ++it){
    TLorentzVector tmp = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());

    if (it->pt()>15 && it->photonID("PhotonCutBasedIDLoose")){
      myPhotons.push_back(tmp);
      FHM->MakePhotonPlots(*it);
    }
    //const std::vector<IdPair> IDs = it->photonIDs();
    //for (size_t t=0; t<IDs.size(); t++)
    //cout<<t<<"  ID name:"<<IDs[t].first<<"  ID pass?="<<IDs[t].second<<endl;
  }


  sort(myPhotons.begin(), myPhotons.end(), P4SortCondition);

  if (myPhotons.size()<2) return;
  TLorentzVector gamma1 = myPhotons[0];
  TLorentzVector gamma2 = myPhotons[1];

  CountEvents(1, "Two Photons reconstructed",ww,fcuts);
  FillHistoCounts(1, ww);


  //edm::Handle<reco::GenJetCollection> genJets;
  //event.getByLabel(edm::InputTag("slimmedGenJets"), genJets);

  vector<TLorentzVector> myJets, myJets24, bJets;
  TLorentzVector bjet1, bjet2;
  
  edm::Handle<std::vector<std::vector<flashgg::Jet>>> jetsCol;
  event.getByLabel(jets_, jetsCol);

  if (jetsCol->size()<=1) return;
 
  UInt_t nJets=0;
  for(UInt_t j = 0 ; j < jetsCol->at( 0 ).size() ; j++ ) {
    pat::Jet jet = jetsCol->at(0)[j];
        
    //std::cout<<jetsCol->at(0)[j].pt()<<"  from pat = "<<jet.pt()<<std::endl;    
    if (jet.pt() > 30){
      TLorentzVector tmp = TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy());
      myJets.push_back(tmp);
      
      if (fabs(jet.eta()) < 2.4)
	myJets24.push_back(tmp);
	
      nJets++;
    }
  }
  
    
  if (myJets.size()<2) return;
  CountEvents(2, "At least two Jets found in event (no b-tag)",ww,fcuts);
  FillHistoCounts(2, ww);



  if (myJets24.size()<2) return;
  CountEvents(3, "At least two Jets in |eta|<2.4",ww,fcuts);
  FillHistoCounts(3, ww);


  sort(myJets24.begin(), myJets24.end(), P4SortCondition);


  FHM->SetGamma1(gamma1);
  FHM->SetGamma2(gamma2);

  FHM->SetBJet1(myJets24[0]);
  FHM->SetBJet2(myJets24[1]);

  FHM->MakeMainHistos(3, ww);

  /*

  if (gen_gamma1.Pt() < 30 || gen_gamma2.Pt() < 30) return;
  if (fabs(gen_gamma1.Eta()) > 2.5 || fabs(gen_gamma2.Eta()) > 2.5) return;

  CountEvents(3, "Photons pT> 30 GeV and |eta|<2.5",ww,fcuts);
  FillHistoCounts(3, ww);
  FHM->MakeMainHistos(3, ww);

  Float_t Mgg = (gen_gamma1 + gen_gamma2).M();
  if ( Mgg < 100 || Mgg > 180) return;

  CountEvents(4, "Photon 100 < m(gg) < 180 GeV",ww,fcuts);
  FillHistoCounts(4, ww);
  FHM->MakeMainHistos(4, ww);

  if ( gen_gamma1.Pt() < Mgg/3 || gen_gamma2.Pt() < Mgg/4) return;

  CountEvents(5, "pT(g1)/m(gg) > 1/3  and pT(g2)/m(gg) > 1/4",ww,fcuts);
  FillHistoCounts(5, ww);
  FHM->MakeMainHistos(5, ww);

  if (gen_bjet1.Pt() < 25 || gen_bjet2.Pt() < 25) return;
  if (fabs(gen_bjet1.Eta()) > 2.5 || fabs(gen_bjet2.Eta()) > 2.5) return;

  CountEvents(6, "The 2 Jets pT > 25 GeV and |eta|<2.5",ww,fcuts);
  FillHistoCounts(6, ww);
  FHM->MakeMainHistos(6, ww);

  Float_t Mbjbj = (gen_bjet1 + gen_bjet2).M();

  if ( Mbjbj < 60 || Mbjbj > 180) return;

  CountEvents(7, "60 < m(jj) < 180 GeV",ww,fcuts);
  FillHistoCounts(7, ww);
  FHM->MakeMainHistos(7, ww);

  */

}

void HHbbggAnalyzer::endJob()
{
  cout<<"\t HHHHHbbbbbggggg \t END JOB in "<<__PRETTY_FUNCTION__<<endl;
  DummyAnalyzer::endJob(2);
}
