#include "../interface/GenAnalyzer.h"

Int_t ncnt = 0;

GenAnalyzer::GenAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs):
  DummyAnalyzer::DummyAnalyzer(cfg, fs)
{
  cout<<"\t GGGGGGGGG \t Constructructor in "<<__PRETTY_FUNCTION__<<endl;
  FHM = new FggHistMakerHHbbgg(hists);
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
  eventNumber  = event.id().event();
  runNumber    = event.id().run();

  if (isRealData) throw cms::Exception("Dude! This is Data. You can't run gen-level analysis on Data...");

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

  edm::Handle<vector<reco::GenParticle> > genParts;
  event.getByLabel( myGen_, genParts );


  std::vector<TLorentzVector> gen_photons, gen_jets;
  TLorentzVector gen_gamma1, gen_gamma2;
  TLorentzVector gen_bjet1, gen_bjet2;
  TLorentzVector gen_bQ1, gen_bQ2;


  //Bool_t yes = 0;
  TLorentzVector tmp;

  // --
  // GEN Particles.
  // Select b-qurks and photons from the Higgses:
  
  for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {

    //if (igen->pdgId()==22 && igen->status()==1 && igen->mother()->pdgId()==25){
    //if (igen->pdgId()==22 && igen->status()==1 && igen->isHardProcess()){
    //if (igen->pdgId()==22 && igen->status()==1 && igen->fromHardProcess()){
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

    //if (yes && igen->pdgId()==22 && igen->status()==1){
    //std::cout<<totEvents<<"\t\t 777 Found a Photon! 777  status="<<igen->status()
    //	 <<"\n IsHard Process = "<<igen->isHardProcess()
    //	 <<" fromHardProcessFinalState = "<<igen->fromHardProcessFinalState()<<std::endl;
    //std::cout<<"\t Its Mother id "<< igen->mother()->pdgId()
    //	 <<"\t Its grandMother id "<< igen->mother()->mother()->pdgId()<<std::endl;
    // //<<"\n uniqueMother ID = "<<igen->uniqueMother().pdgId()<<std::endl;
    //}
  }

  /*
  if (gen_photons.size()<2){
    ncnt++;
    for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {
      
      hists->fill1DHist(igen->pdgId(),"pdg_id","All PDGs", 40,-20,20, ww, "PDG");
      
      if (igen->pdgId()==25){
	std::cout<<ncnt<<"\t\t 777 Found the HIGGS! 777"<<std::endl;
	//std::cout<<"\t Its *mother* is: "<<igen->mother()->pdgId()<<std::endl;
	//std::cout<<"\t Its *daughters* are: "<<igen->daughter(0)->pdgId()<<" and "<<igen->daughter(1)->pdgId()
	//	 <<"\n \t and their statuses are "<<igen->daughter(0)->status()<<" and "<<igen->daughter(1)->status()<<std::endl;


	for (UInt_t d=0; d<2; d++){
	  if (igen->daughter(d)->pdgId()==22){
	    std::cout<<d<<" daugheter status = "<<igen->daughter(d)->status()<<std::endl;
	    for (UInt_t dd=0; dd < igen->daughter(d)->numberOfDaughters(); dd++){
	      std::cout<<"  grand-daugheter "<<dd<<"  PDG = "<<igen->daughter(d)->daughter(dd)->pdgId()
		       <<"  status = "<<igen->daughter(d)->daughter(dd)->status()
		       <<" isPromptFS="<<igen->isPromptFinalState()<<" fromHardFS="<<igen->fromHardProcessFinalState()<<std::endl;
	      for (UInt_t ddd=0; ddd < igen->daughter(d)->daughter(dd)->numberOfDaughters(); ddd++){
		std::cout<<"\t grand-daugheter "<<ddd<<"  PDG"<<igen->daughter(d)->daughter(dd)->daughter(ddd)->pdgId()
			 <<" status = "<<igen->daughter(d)->daughter(dd)->daughter(ddd)->status()<<std::endl;
	      }
	      
	    }
	  }
	}
      }
    }
 }
  */
 


  
  if (gen_photons.size()<2) return;
  CountEvents(1, "Two gen photons from HardProcess",ww,fcuts);
  FillHistoCounts(1, ww);

  sort(gen_photons.begin(), gen_photons.end(), P4SortCondition);

  gen_gamma1 = gen_photons[0];
  gen_gamma2 = gen_photons[1];

  gen_bjet1 = gen_bQ1;
  gen_bjet2 = gen_bQ2;
  
  FHM->SetGamma1(gen_gamma1);
  FHM->SetGamma2(gen_gamma2);

  // For now set hte b-jets to be the b-quarks:
  FHM->SetBJet1(gen_bjet1);
  FHM->SetBJet2(gen_bjet2);


  FHM->MakeMainHistos(1, ww);

  hists->fill1DHist(gen_bQ1.Pt(), Form("04_gen_bQ1_pt_%s_cut%i", "", 1),";p_{T}^{b}", 100,0,300,  ww, "GEN");
  hists->fill1DHist(gen_bQ2.Pt(), Form("04_gen_bQ2_pt_%s_cut%i", "", 1),";p_{T}^{#bar{b}}",100,0,300,  ww, "GEN");
  
  hists->fill1DHist(gen_bQ1.Eta(), Form("04_gen_bQ1_eta_%s_cut%i", "", 1),";#eta^{b}", 50,-5,5,  ww, "GEN");
  hists->fill1DHist(gen_bQ2.Eta(), Form("04_gen_bQ2_eta_%s_cut%i", "", 1),";#eta^{#bar{b}}",50,-5,5,  ww, "GEN");


  // --
  // GEN JETS: match them to b-quarks
  // --
  edm::Handle<reco::GenJetCollection> genJets;
  event.getByLabel(edm::InputTag("slimmedGenJets"), genJets);

  UInt_t jetInd = 0;
  Float_t dRb1 = 9999, dRb2 = 9999;

  for( reco::GenJetCollection::const_iterator jet = genJets->begin(); jet != genJets->end(); ++jet ) {

    if (jet->pt() < 25) continue;
    //std::cout <<"Gen Jet: "<<jetInd<<", pt="<<jet->pt()<<", eta="<<jet->eta()<<", phi="<<jet->phi() <<std::endl;
    tmp.SetPxPyPzE(jet->px(), jet->py(), jet->pz(), jet->energy());
    gen_jets.push_back(tmp);


    //UInt_t ndau = jet->numberOfDaughters();
    //std::cout<<"Daughters: "<<jet->numberOfDaughters()<<std::endl;

    //if (fabs(jet->eta())<2){ 
      //std::cout<<"\t XXX Event number = "<<eventNumber<<std::endl;
      //for (UInt_t n = 0; n<ndau; n++){ 
	//CandidatePtr daughterPtr( size_type i ) const 
	
	//const reco::Candidate *jetDau = jet->daughter(n);      
	//std::cout<<n<<"  PDG="<<jet->daughter(n)->pdgId()<<"  status="<<jet->daughter(n)->status()<<std::endl;
	//std::cout<<n<<"  PDG="<<jetDau->pdgId()<<"  status="<<jetDau->status()<<std::endl;
    //}
    //}

    // Does not work from Gen-jets:
    //std::cout<<"Jet Hadron Flavour "<<jet->hadronFlavour()
    //<<"Parton Flavour "<<jet->partonFlavour()<<std::endl;
    
    if (tmp.DeltaR(gen_bQ1) < dRb1) {
      gen_bjet1 = tmp;
      dRb1 = tmp.DeltaR(gen_bQ1);
    }
    if (tmp.DeltaR(gen_bQ2) < dRb2){
      gen_bjet2 = tmp;
      dRb2 = tmp.DeltaR(gen_bQ2);
    }

    //std::vector <const reco::GenParticle*> mcparts = jet->getGenConstituents();
    //for (size_t i = 0; i < mcparts.size(); i++) {
    //const reco::GenParticle* mcpart = mcparts[i];
    //if (mcpart) {
    //std::cout << "      #" << i << "  PDG code:" << mcpart->pdgId()
    //	    << ", p/pt/eta/phi: " << mcpart->p() << '/' << mcpart->pt() << '/' << mcpart->eta() << '/' << mcpart->phi() << std::endl;
    //}
    //else {
    //std::cout << "      #" << i << "  No information about constituent" << std::endl;
    //}

    //}

    jetInd++;

  }

  hists->fill1DHist(dRb1, "dR_jetb1",";#DeltaR(jet, b-quark)", 50, 0,1.2, ww, "GEN");
  hists->fill1DHist(dRb2, "dR_jetb2",";#DeltaR(jet, #bar{b}-quark)", 50, 0,1.2, ww, "GEN");

  Float_t Mbb = (gen_bQ1+gen_bQ2).M();
  Float_t dRbb = gen_bQ1.DeltaR(gen_bQ2);

  hists->fill1DHist(Mbb, "Mbb",";m(bb), GeV", 100, 124.5, 125.5, ww, "GEN");
  hists->fill1DHist(dRbb, "dR_bb",";#DeltaR(b, #bar{b})", 50, 0,5, ww, "GEN");


  hists->fill2DHist(gen_bQ1.Eta()-gen_bQ2.Eta(), gen_bQ1.DeltaPhi(gen_bQ2), 
		    "dEta_dPhi",";#DeltaEta(b, #bar{b});#DeltaPhi(b, #bar{b})", 50, -5,5, 50, -3.15, 3.15, ww, "GEN");

  
  if (gen_jets.size()<2) return;
  CountEvents(2, "Two gen Jets in event",ww,fcuts);
  FillHistoCounts(2, ww);

  sort(gen_jets.begin(), gen_jets.end(), P4SortCondition);

  //gen_bjet1 = gen_jets[0];
  //gen_bjet2 = gen_jets[1];

  // Here swithch to the actual jets.
  // The jets aer those that match to the given quark 1=b, 2=bbar 
  FHM->SetBJet1(gen_bjet1);
  FHM->SetBJet2(gen_bjet2);

  FHM->MakeMainHistos(2, ww);


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

}

void GenAnalyzer::endJob()
{
  cout<<"\t GGGGGGGGG \t END JOB in "<<__PRETTY_FUNCTION__<<endl;
  DummyAnalyzer::endJob(2);
}
