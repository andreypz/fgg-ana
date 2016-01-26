#include "../interface/HHbbggAnalyzer.h"

//typedef std::pair<std::string,bool> IdPair;
Bool_t doGen=0;

HHbbggAnalyzer::HHbbggAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs):
  DummyAnalyzer::DummyAnalyzer(cfg, fs),
  //inputTagJets_(cfg.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
  phoIDcutEB_(cfg.getUntrackedParameter<std::vector<double > >("phoIDcutEB") ),
  phoIDcutEE_(cfg.getUntrackedParameter<std::vector<double > >("phoIDcutEE") ),
  phoISOcutEB_(cfg.getUntrackedParameter<std::vector<double > >("phoISOcutEB") ),
  phoISOcutEE_(cfg.getUntrackedParameter<std::vector<double > >("phoISOcutEE") ),
  cutFlow_(cfg.getUntrackedParameter<UInt_t>("cutFlow") ),
  useDiPhotons_(cfg.getUntrackedParameter<Bool_t>("useDiPhotons") ),
  diPhotons_(cfg.getParameter<edm::InputTag>("diPhotonTag")),
  phoIDtype_(cfg.getUntrackedParameter<UInt_t>("phoIDtype") )
{

  cout<<red<<"\t HHHHbbbbgggg \t Constructructor in "<<__PRETTY_FUNCTION__<<def<<endl;
  FHM = new FggHistMakerHHbbgg(hists);
  tools = new bbggTools();

  bTag = "pfCombinedInclusiveSecondaryVertexV2BJetTags";

  rhoFixedGrid_ = edm::InputTag( "fixedGridRhoAll" ) ;

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


  edm::Handle<double> rhoHandle;
  event.getByLabel( rhoFixedGrid_, rhoHandle );
  const double rhoFixedGrd = *( rhoHandle.product() );

  FHM->Reset(rhoFixedGrd, 1, eventNumber);
  tools->setRho(rhoFixedGrd);

  if (!isRealData && doGen) {
    //&& sdample = signal...

    edm::Handle<vector<reco::GenParticle> > genParts;
    event.getByLabel( myGen_, genParts );

    std::vector<TLorentzVector> gen_photons, gen_jets;
    TLorentzVector gen_gamma1, gen_gamma2;
    //TLorentzVector gen_bjet1, gen_bjet2;
    TLorentzVector gen_bQ1, gen_bQ2;
    TLorentzVector tmp;

    for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {

      if (abs(igen->pdgId())==5 && igen->mother()->pdgId()==25) {

	//std::cout<<totEvents<<"\t"<<igen->pdgId()<<"\t\t BBB Found b-quark from Higgs! BBB  Its status = "<<igen->status()
	//<<" Its charge = "<<igen->charge()<<std::endl;
	tmp.SetPxPyPzE(igen->px(), igen->py(), igen->pz(), igen->energy());
	if (igen->pdgId()==5)
	  gen_bQ1 = tmp;
	if (igen->pdgId()==-5)
	  gen_bQ2 = tmp;
      }

      if (igen->pdgId()==22 && igen->isPromptFinalState()) {
	tmp.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());
	gen_photons.push_back(tmp);

      }

    }

    /*
    sort(gen_photons.begin(), gen_photons.end(), P4SortCondition);

    gen_gamma1 = gen_photons[0];
    gen_gamma2 = gen_photons[1];

    gen_bjet1 = gen_bQ1;
    gen_bjet2 = gen_bQ2;

    FHM->SetGamma1(gen_gamma1);
    FHM->SetGamma2(gen_gamma2);

    FHM->MakeMainHistos(0, ww, "GEN");
    */
  }

  vector<TLorentzVector> myLeptons, myPhotons;


  edm::Handle<std::vector<flashgg::Photon> > photons;
  event.getByLabel(photons_, photons);

  // Verteces
  edm::Handle<reco::VertexCollection> primaryVtcs;
  event.getByLabel(edm::InputTag("offlineSlimmedPrimaryVertices"), primaryVtcs);
  if (primaryVtcs->size()<1) return;
  const edm::Ptr<reco::Vertex> CandVtx(primaryVtcs, 0);

  edm::Handle<std::vector<flashgg::DiPhotonCandidate> > diPhotons;
  event.getByLabel(diPhotons_, diPhotons);

  // This is to store the selected di-photons:
  std::vector<flashgg::DiPhotonCandidate> myDiPho;

  if (useDiPhotons_) {
    for(std::vector<flashgg::DiPhotonCandidate>::const_iterator it=diPhotons->begin(); it!=diPhotons->end(); ++it){

      const flashgg::Photon *p1, *p2;
      p1 = it->leadingPhoton();
      p2 = it->subLeadingPhoton();

      if (p2->pt() < 15) continue;

      FHM->MakePhotonPlots(*p1, "DiPho-Lead");
      FHM->MakePhotonPlots(*p2, "DiPho-Sub");

      Bool_t passID = false;
      Bool_t passISO = false;
      switch ( phoIDtype_ ) {
      case 1 :
        // Cut based ID from bbggTools
	// Leading photon:
        passID = ( (p1->isEB() && tools->isPhoID(&(*p1), phoIDcutEB_)) ||
		   (p1->isEE() && tools->isPhoID(&(*p1), phoIDcutEE_))
		   );

	passISO = ((p1->isEB() && tools->isPhoISO(&(*it), 0, phoISOcutEB_)) ||
		   (p1->isEE() && tools->isPhoISO(&(*it), 0, phoISOcutEE_))
		   );
	// SubLeading photon:
        passID = passID && (( p2->isEB() && tools->isPhoID(&(*p2), phoIDcutEB_)) ||
			    ( p2->isEE() && tools->isPhoID(&(*p2), phoIDcutEE_))
			    );
	passISO = passISO && ((p2->isEB() && tools->isPhoISO(&(*it), 1, phoISOcutEB_)) ||
			      (p2->isEE() && tools->isPhoISO(&(*it), 1, phoISOcutEE_))
			      );

	passID = passID && passISO;
	break;

      case 2 :
	// Cut based ID stored in AOD/PAT/fgg
        passID = (p1->photonID("cutBasedPhotonID-Spring15-50ns-V1-standalone-medium") &&
		  p2->photonID("cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"));
        break;

      case 3 :
	// EGamma MVA ID stored in the PAT object

	// Cut on the value
	//passID = (( it->isEB() && it->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values") > 0.374 ) ||
	//        ( it->isEE() && it->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values") > 0.336 ));
	// Use the boolean
	passID = (p1->photonID("mvaPhoID-Spring15-25ns-nonTrig-V2-wp90") &&
		  p2->photonID("mvaPhoID-Spring15-25ns-nonTrig-V2-wp90"));
	break;

      case 4 :
        // Hgg MVA ID:
        //https://github.com/cms-analysis/flashgg/blob/f1ad5f8c6d9b02481656c8beec8fc5e10ffa35b7/DataFormats/interface/Photon.h#L136
        //// if lazy flag is true only compare key (needed since fwlite does not fill provenance info)
        passID = ((p1->phoIdMvaDWrtVtx(CandVtx, true) > 0.2) &&
		  (p2->phoIdMvaDWrtVtx(CandVtx, true) > 0.2));

        break;

      case 5 :
	// HEEP ID: TBD
        break;

      default : break;

      }

      passID = passID && p1->passElectronVeto() && p2->passElectronVeto();

      passID = passID && tools->HggHLTpreselection(&(*it));
      
      if (!passID) continue;

      myDiPho.push_back(*it);

      /* there is probably double-counting going on here. need to check
	 TLorentzVector tmp = TLorentzVector(p1->px(), p1->py(), p1->pz(), p1->energy());
	 myPhotons.push_back(tmp);
	 tmp = TLorentzVector(p2->px(), p2->py(), p2->pz(), p2->energy());
	 myPhotons.push_back(tmp);
      */
    }
  }
  else { // Use regular photons
    for(std::vector<flashgg::Photon>::const_iterator it=photons->begin(); it!=photons->end(); ++it){
      TLorentzVector tmp = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());

      Bool_t passID = false;

      if (it->pt() < 15) continue;

      FHM->MakePhotonPlots(*it);

      //std::vector<std::string> labels = it->userFloatNames();
      //for (size_t l = 0; l<labels.size(); l++)
      //cout<<labels[l]<<endl;
      //return;

      switch ( phoIDtype_ ) {
      case 1 :
        // Cut based ID from bbggTools
	// (Notice: there is no ISO here. It only works on diPhotons with bbggTools)
        passID = (( it->isEB() && tools->isPhoID(&(*it), phoIDcutEB_) ) ||
		  ( it->isEE() && tools->isPhoID(&(*it), phoIDcutEE_) ));
        break;

      case 2 :

        // Cut based ID stored in AOD/PAT/fgg
        passID = it->photonID("cutBasedPhotonID-Spring15-50ns-V1-standalone-medium");
        //passID = it->photonID("cutBasedPhotonID-Spring15-25ns-V1-standalone-medium");
        break;

      case 3 :
        // EGamma MVA ID stored in the PAT object

        // Cut on the value
        //passID = (( it->isEB() && it->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values") > 0.374 ) ||
        //	( it->isEE() && it->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values") > 0.336 ));
        // Use the boolean
        passID = it->photonID("mvaPhoID-Spring15-25ns-nonTrig-V2-wp90");
        break;

      case 4 :
        // Hgg MVA ID:
        //https://github.com/cms-analysis/flashgg/blob/f1ad5f8c6d9b02481656c8beec8fc5e10ffa35b7/DataFormats/interface/Photon.h#L136
        //// if lazy flag is true only compare key (needed since fwlite does not fill provenance info)
        passID = (it->phoIdMvaDWrtVtx(CandVtx, true) > 0.2);
	break;

      default : break;
      }

      //const std::vector<IdPair> IDs = it->photonIDs();
      //for (size_t t=0; t<IDs.size(); t++)
      //cout<<t<<"  ID name:"<<IDs[t].first<<"  ID pass?="<<IDs[t].second<<endl;

      passID = passID && it->passElectronVeto();
      if (!passID) continue;

      myPhotons.push_back(tmp);
    }
  }


  TLorentzVector gamma1, gamma2;

  if (useDiPhotons_){
    if (myDiPho.size()<1) return;
    const flashgg::Photon *p1, *p2;
    p1 = myDiPho[0].leadingPhoton();
    p2 = myDiPho[0].subLeadingPhoton();
    gamma1 = TLorentzVector(p1->px(), p1->py(), p1->pz(), p1->energy());
    gamma2 = TLorentzVector(p2->px(), p2->py(), p2->pz(), p2->energy());
  }
  else{

    sort(myPhotons.begin(), myPhotons.end(), P4SortCondition);
    if (myPhotons.size()<2) return;
    gamma1 = myPhotons[0];
    gamma2 = myPhotons[1];
  }

  CountEvents(1, "Two Photons reconstructed",ww,fcuts);
  FillHistoCounts(1, ww);


  //edm::Handle<reco::GenJetCollection> genJets;
  //event.getByLabel(edm::InputTag("slimmedGenJets"), genJets);

  vector<TLorentzVector> myJets, myJets25;
  vector<flashgg::Jet> bJets;
  TLorentzVector bjet1, bjet2;


  edm::Handle<std::vector<std::vector<flashgg::Jet>>> jetsCol;
  event.getByLabel(jets_, jetsCol);

  if (jetsCol->size()<1) return;

  UInt_t nJets=0;
  for(UInt_t j = 0 ; j < jetsCol->at( 0 ).size() ; j++ ) {
    flashgg::Jet jet = jetsCol->at(0)[j];

    //edm::Ptr<flashgg::Jet> jet(jetsCol, 0); Does not work

    //if (!tools->isJetID( &jet )) continue;
    if (!jet.passesJetID(flashgg::JetIDLevel::Loose)) continue;

    //std::cout<<j<<" Pass jet ID. Pt="<<jet.pt()<<std::endl;
    //std::cout<<j<<" Does not Pass jet ID. Pt="<<jet.pt()<<std::endl;


    // does not work in flashgg:
    //if (!jet.passesPuJetId( need Vtx here)) continue;

    //std::cout<<jetsCol->at(0)[j].pt()<<"  from pat = "<<jet.pt()<<std::endl;

    if (jet.pt() > 25){

      TLorentzVector tmp = TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy());
      myJets.push_back(tmp);

      if (fabs(jet.eta()) > 2.5) continue;
      if (tmp.DeltaR(gamma1) < 0.4 || tmp.DeltaR(gamma2) < 0.4) continue;

      myJets25.push_back(tmp);

      //
      // Order the Jets by the b-discriminator value:
      //
      //Only kep the b-tagged jets:
      if (jet.bDiscriminator(bTag) < 0) continue;
      // Begin with putting the first jet in the array
      if (bJets.size()==0){
	bJets.push_back(jet);
	continue;
      }
      // Now loop over all and order them by b-tag discriminator
      for (std::vector<flashgg::Jet>::const_iterator b = bJets.begin(); b!= bJets.end(); b++) {
	auto nx = std::next(b);
	if ( jet.bDiscriminator(bTag) > b->bDiscriminator(bTag)) {
	  bJets.insert(b, jet);
	  break;
	}
	else if (nx == bJets.end()){
	  bJets.push_back(jet);
	  break;
	}
      }
      // END of b-jet ordering
    }

    nJets++;

  }


  /*
  // Checks:
  for (size_t b = 0; b < myJets25.size(); b++) {
  cout<<b<<" before b-tag sorting  pt="<<myJets25[b].Pt()<<endl;
  }
  //sort(bJets.begin(), bJets.end(), bTagSort);
  for (size_t b = 0; b < bJets.size(); b++) {
  cout<<b<<" \t after sorting bdisc = "<< bJets[b].bDiscriminator(bTag) <<"  pt="<<bJets[b].pt()<<endl;
  }
  */


  FHM->SetGamma1(gamma1);
  FHM->SetGamma2(gamma2);
  Float_t Mgg = (gamma1 + gamma2).M();

  if (bJets.size()>=1){
    bjet1.SetPxPyPzE(bJets[0].px(), bJets[0].py(), bJets[0].pz(), bJets[0].energy());
    FHM->SetBJet1(bjet1);
    if (bJets.size()>=2) {
      bjet2.SetPxPyPzE(bJets[1].px(), bJets[1].py(), bJets[1].pz(), bJets[1].energy());
      FHM->SetBJet2(bjet2);
    }
  }

  Float_t Mbjbj = (bjet1 + bjet2).M();

  switch ( cutFlow_ ) {
  case 1 :
    // NCU cutflow

    if (myJets.size()<2) return;
    CountEvents(2, "At least two Jets w/ pT>25 (no b-tag)",ww,fcuts);
    FillHistoCounts(2, ww);
    FHM->MakeNPlots(2, myPhotons.size(), myJets.size(), bJets.size(), ww);


    if (myJets25.size()<2) return;
    CountEvents(3, "At least two Jets w/ pT>25 and |eta|<2.5",ww,fcuts);
    FillHistoCounts(3, ww);
    //FHM->MakeMainHistos(3, ww);

    sort(myJets25.begin(), myJets25.end(), P4SortCondition);


    if (bJets.size()<2) return;

    CountEvents(4, "Two jets with bDissc > 0",ww,fcuts);
    FillHistoCounts(4, ww);
    FHM->MakeMainHistos(4, ww);

    CountEvents(5, "Reserved",ww,fcuts);
    FillHistoCounts(5, ww);
    //FHM->MakeMainHistos(5, ww);


    FHM->MakeJetPlots(bJets[0], "bJet-Lead");
    FHM->MakeJetPlots(bJets[1], "bJet-Sub");


    if (gamma1.DeltaR(gamma2) < 0.4) return;
    if (gamma1.Pt() < 30 || gamma2.Pt() < 30) return;
    if (fabs(gamma1.Eta()) > 2.5 || fabs(gamma2.Eta()) > 2.5) return;

    CountEvents(6, "Two Photons pT>30GeV, |eta|<2.5 and dR>0.4",ww,fcuts);
    FillHistoCounts(6, ww);
    FHM->MakeMainHistos(6, ww);


    if ( gamma1.Pt() < Mgg/3 || gamma2.Pt() < Mgg/4) return;
    CountEvents(7, "pT(g1)/m(gg) > 1/3  and pT(g2)/m(gg) > 1/4",ww,fcuts);
    FillHistoCounts(7, ww);
    FHM->MakeMainHistos(7, ww);


    if ( Mgg < 100 || Mgg > 180) return;
    CountEvents(8, "100 < m(gg) < 180 GeV",ww,fcuts);
    FillHistoCounts(8, ww);
    FHM->MakeMainHistos(8, ww);


    if (bjet1.Pt() < 25 || bjet2.Pt() < 25) return;
    if (fabs(bjet1.Eta()) > 2.5 || fabs(bjet2.Eta()) > 2.5) return;

    CountEvents(9, "The 2 Jets pT > 25 GeV and |eta|<2.5",ww,fcuts);
    FillHistoCounts(9, ww);
    FHM->MakeMainHistos(9, ww);

    FHM->MakeNPlots(9, myPhotons.size(), myJets.size(), bJets.size(), ww);


    if ( Mbjbj < 60 || Mbjbj > 180) return;

    CountEvents(10, "60 < m(jj) < 180 GeV",ww,fcuts);
    FillHistoCounts(10, ww);
    FHM->MakeMainHistos(10, ww);


    break;


  case 2 :
    // Rafaels' cutflow

    if (gamma1.DeltaR(gamma2) < 0.4) return;
    if (gamma1.Pt() < 30 || gamma2.Pt() < 30) return;
    if (fabs(gamma1.Eta()) > 2.5 || fabs(gamma2.Eta()) > 2.5) return;

    CountEvents(2, "Two Photons pT>30GeV, |eta|<2.5 and dR>0.4",ww,fcuts);
    FillHistoCounts(2, ww);
    FHM->MakeMainHistos(2, ww);


    if ( gamma1.Pt() < Mgg/3 || gamma2.Pt() < Mgg/4) return;
    CountEvents(3, "pT(g1)/m(gg) > 1/3  and pT(g2)/m(gg) > 1/4",ww,fcuts);
    FillHistoCounts(3, ww);
    FHM->MakeMainHistos(3, ww);


    if ( Mgg < 100 || Mgg > 180) return;
    CountEvents(4, "100 < m(gg) < 180 GeV",ww,fcuts);
    FillHistoCounts(4, ww);
    FHM->MakeMainHistos(4, ww);


    if (myJets25.size()<2) return;
    CountEvents(5, "At least two Jets w/ pT>25 and |eta|<2.5",ww,fcuts);
    FillHistoCounts(5, ww);
    FHM->MakeNPlots(5, myPhotons.size(), myJets.size(), bJets.size(), ww);

    sort(myJets25.begin(), myJets25.end(), P4SortCondition);

    if (bJets.size()<2) return;

    CountEvents(6, "Two jets with bDissc > 0",ww,fcuts);
    FillHistoCounts(6, ww);
    FHM->MakeMainHistos(6, ww);


    FHM->MakeJetPlots(bJets[0], "bJet-Lead");
    FHM->MakeJetPlots(bJets[1], "bJet-Sub");

    if (bjet1.Pt() < 25 || bjet2.Pt() < 25) return;
    if (fabs(bjet1.Eta()) > 2.5 || fabs(bjet2.Eta()) > 2.5) return;

    CountEvents(7, "The 2 Jets pT > 25 GeV and |eta|<2.5",ww,fcuts);
    FillHistoCounts(7, ww);
    FHM->MakeMainHistos(7, ww);

    FHM->MakeNPlots(7, myPhotons.size(), myJets.size(), bJets.size(), ww);


    if ( Mbjbj < 60 || Mbjbj > 180) return;

    CountEvents(8, "60 < m(jj) < 180 GeV",ww,fcuts);
    FillHistoCounts(8, ww);
    FHM->MakeMainHistos(8, ww);


    break;

  default :
    return;
  }

}

void HHbbggAnalyzer::endJob()
{
  cout<<"\t HHHHHbbbbbggggg \t END JOB in "<<__PRETTY_FUNCTION__<<endl;
  DummyAnalyzer::endJob(2);
}
