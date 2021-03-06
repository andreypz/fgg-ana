#include "../interface/HHbbggAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"


//typedef std::pair<std::string,bool> IdPair;
const UInt_t hisEVTS[] = {33842,40479};
Int_t evSize = sizeof(hisEVTS)/sizeof(int);

HHbbggAnalyzer::HHbbggAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs):
  DummyAnalyzer::DummyAnalyzer(cfg, fs),
  //vertexes_( cfg.getParameter<edm::InputTag> ( "vertexes" ) ),
  myTriggers_(cfg.getUntrackedParameter<std::vector<std::string> >("myTriggers")),
  //inputTagJets_(cfg.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
  phoIDcutEB_(cfg.getUntrackedParameter<std::vector<double > >("phoIDcutEB") ),
  phoIDcutEE_(cfg.getUntrackedParameter<std::vector<double > >("phoIDcutEE") ),
  phoISOcutEB_(cfg.getUntrackedParameter<std::vector<double > >("phoISOcutEB") ),
  phoISOcutEE_(cfg.getUntrackedParameter<std::vector<double > >("phoISOcutEE") ),
  phoISO_nhCorEB_(cfg.getUntrackedParameter<std::vector<double > >("phoISO_nhCorEB") ),
  phoISO_nhCorEE_(cfg.getUntrackedParameter<std::vector<double > >("phoISO_nhCorEE") ),
  phoISO_phCorEB_(cfg.getUntrackedParameter<std::vector<double > >("phoISO_phCorEB") ),
  phoISO_phCorEE_(cfg.getUntrackedParameter<std::vector<double > >("phoISO_phCorEE") ),
  cutFlow_(cfg.getUntrackedParameter<UInt_t>("cutFlow") ),
  useDiPhotons_(cfg.getUntrackedParameter<Bool_t>("useDiPhotons") ),
  diPhotons_(cfg.getParameter<edm::InputTag>("diPhotonTag") ),
  phoIDtype_(cfg.getUntrackedParameter<UInt_t>("phoIDtype") ),
  doBJetRegression_(cfg.getUntrackedParameter<Bool_t>("doBJetRegression") ),
  bRegFile_(cfg.getUntrackedParameter<edm::FileInPath>("bRegFile") ),
  doNonResWeights_(cfg.getUntrackedParameter<Bool_t>("doNonResWeights") )
{

  cout<<red<<"\t HHHHbbbbgggg \t Constructructor in "<<__PRETTY_FUNCTION__<<def<<endl;

  FHM = new FggHistMakerHHbbgg(hists);
  tools  = new bbggTools();
  angles = new Angles();
  jetReg = new bbggJetRegression();
  NRW    = new NonResWeights();

  jetReg->SetupRegression("BDTG method", bRegFile_.fullPath().data());

  bTagName = "pfCombinedInclusiveSecondaryVertexV2BJetTags";

  rhoFixedGrid_ = edm::InputTag( "fixedGridRhoAll" ) ;

  globVar_ = new flashgg::GlobalVariablesDumper(cfg);

  flatTree = fs.make<TTree>("bbggSelectionTree","Flat tree for HH->bbgg analysis");
  std::map<std::string, std::string> replacements; // It's not doing nothing..
  globVar_->bookTreeVariables(flatTree, replacements); // This fills event and run numbers
  //flatTree->Branch("run", &o_run, "o_run/i");
  //flatTree->Branch("evt", &o_evt, "o_evt/l");

  // These are actually for Limit trees (can't we have them in one single tree)?
  //flatTree->Branch("cut_based_ct", &o_category, "o_category/B"); //0: 2btag, 1: 1btag
  //flatTree->Branch("evWeight", &o_weight, "o_weight/D");
  flatTree->Branch("mjj",  &o_bbMass, "o_bbMass/D");
  flatTree->Branch("mgg",  &o_ggMass, "o_ggMass/D");
  flatTree->Branch("mtot", &o_bbggMass, "o_bbggMass/D"); //

  flatTree->Branch("gen_mHH",      &gen_mHH,      "gen_mHH/D");
  flatTree->Branch("gen_cosTheta", &gen_cosTheta, "gen_cosTheta/D");
  //if (doNonResWeights_)
  //flatTree->Branch("NRWeights", NRWeights, "NRWeights[1519]/F");
  // Here is all vars from bbggTools code:
  flatTree->Branch("genWeights", &genWeights);
  flatTree->Branch("genTotalWeight", &genTotalWeight, "genTotalWeight/D");
  flatTree->Branch("leadingPhoton", &leadingPhoton);
  flatTree->Branch("leadingPhotonID", &leadingPhotonID);
  flatTree->Branch("leadingPhotonISO", &leadingPhotonISO);
  flatTree->Branch("leadingPhotonEVeto", &leadingPhotonEVeto, "leadingPhotonEVeto/I");
  flatTree->Branch("leadingPhotonIDMVA", &leadingPhotonIDMVA, "leadingPhotonIDMVA/f");
  flatTree->Branch("subleadingPhoton", &subleadingPhoton);
  flatTree->Branch("subleadingPhotonID", &subleadingPhotonID);
  flatTree->Branch("subleadingPhotonISO", &subleadingPhotonISO);
  flatTree->Branch("subleadingPhotonEVeto", &subleadingPhotonEVeto, "subleadingPhotonEVeto/I");
  flatTree->Branch("subleadingPhotonIDMVA", &subleadingPhotonIDMVA, "subleadingPhotonIDMVA/f");
  flatTree->Branch("diphotonCandidate", &diphotonCandidate);
  flatTree->Branch("nPromptInDiPhoton", &nPromptInDiPhoton, "nPromptInDiPhoton/I");
  flatTree->Branch("leadingJet", &leadingJet);
  flatTree->Branch("leadingJet_KF", &leadingJet_KF);
  flatTree->Branch("leadingJet_Reg", &leadingJet_Reg);
  flatTree->Branch("leadingJet_RegKF", &leadingJet_RegKF);
  flatTree->Branch("leadingJet_bDis", &leadingJet_bDis, "leadingJet_bDis/F");
  flatTree->Branch("leadingJet_flavour", &leadingJet_flavour, "leadingJet_flavour/I");
  flatTree->Branch("leadingJet_hadFlavour", &leadingJet_hadFlavour, "leadingJet_hadFlavour/I");
  flatTree->Branch("subleadingJet", &subleadingJet);
  flatTree->Branch("subleadingJet_KF", &subleadingJet_KF);
  flatTree->Branch("subleadingJet_Reg", &subleadingJet_Reg);
  flatTree->Branch("subleadingJet_RegKF", &subleadingJet_RegKF);
  flatTree->Branch("subleadingJet_bDis", &subleadingJet_bDis, "subleadingJet_bDis/F");
  flatTree->Branch("subleadingJet_flavour", &subleadingJet_flavour, "subleadingJet_flavour/I");
  flatTree->Branch("subleadingJet_hadFlavour", &subleadingJet_hadFlavour, "subleadingJet_hadFlavour/I");
  flatTree->Branch("dijetCandidate", &dijetCandidate);
  flatTree->Branch("dijetCandidate_KF", &dijetCandidate_KF);
  flatTree->Branch("dijetCandidate_Reg", &dijetCandidate_Reg);
  flatTree->Branch("dijetCandidate_RegKF", &dijetCandidate_RegKF);
  flatTree->Branch("diHiggsCandidate", &diHiggsCandidate);
  flatTree->Branch("diHiggsCandidate_KF", &diHiggsCandidate_KF);
  flatTree->Branch("diHiggsCandidate_Reg", &diHiggsCandidate_Reg);
  flatTree->Branch("diHiggsCandidate_RegKF", &diHiggsCandidate_RegKF);
  flatTree->Branch("isSignal", &isSignal, "isSignal/I");
  flatTree->Branch("isPhotonCR", &isPhotonCR, "isPhotonCR/I");
  flatTree->Branch("CosThetaStar", &CosThetaStar, "CosThetaStar/F");
  flatTree->Branch("TriggerResults", &myTriggerResults);
  flatTree->Branch("DiJetDiPho_DR_1", &DiJetDiPho_DR_1, "DiJetDiPho_DR_1/F");
  flatTree->Branch("DiJetDiPho_DR_2", &DiJetDiPho_DR_2, "DiJetDiPho_DR_2/F");
  flatTree->Branch("PhoJetMinDr", &PhoJetMinDr, "PhoJetMinDr/F");



  nodesOfHH = false;
  if (runSample_=="HH_SM" || runSample_.find("HH_NonRes")!=std::string::npos) {
    nodesOfHH = true;

    genTree = fs.make<TTree>("GenTree","A tree for signal reweighting study");
    genTree->Branch("run", &o_run, "run/i");
    genTree->Branch("evt", &o_evt, "evt/l");

    genTree->Branch("mHH",  &gen_mHH,  "mHH/D");
    genTree->Branch("ptH1", &gen_ptH1, "ptH1/D");
    genTree->Branch("ptH2", &gen_ptH2, "ptH2/D");
    genTree->Branch("cosTheta",  &gen_cosTheta, "cosTheta/D");
    genTree->Branch("cosTheta2", &gen_cosTheta2, "cosTheta2/D");

    if (doNonResWeights_) {
      genTree->Branch("NRWeights", NRWeights, "NRWeights[1519]/F");

      /*
      //std::string::size_type sz;   // alias of size_t
      if (runSample_=="HH_SM") nodeFileNum=0;
      else if (runSample_=="HH_NonRes_Box") nodeFileNum=1;
      //else if (runSample_.find("HH_NonRes_")!=std::string::npos)
      //cout<<runSample_.substr(10,1)<<"  "<<runSample_.substr(10,2)<<endl;
      else if (runSample_.size()==11) nodeFileNum = stoi(runSample_.substr(10,1)); //One digit file number: 2-9
      else if (runSample_.size()==12) nodeFileNum = stoi(runSample_.substr(10,2)); //Two digits: 10,11,12,13
      else nodeFileNum = 99;
      cout<<"\t sample = "<<runSample_<<"    nodeFileNum="<<nodeFileNum<<endl;

      //flatTree->Branch("file", &nodeFileNum, "file/i"); //
      //genTree->Branch("file",  &nodeFileNum, "file/i"); //
      */

      std::string fileNameWei1 = edm::FileInPath("APZ/fgg-ana/data/weights_v1_1507_points.root").fullPath();
      //NRwFile = new TFile(fileNameWei1.c_str(), "OPEN");

      // This file includes 12 benchmark v3 weights.
      std::string fileNameWei2 = edm::FileInPath("APZ/fgg-ana/data/weights_v3_bench12_points.root").fullPath();
      //NRwFile2 = new TFile(fileNameWei2.c_str(), "OPEN");

      NRW->LoadHists(fileNameWei1, fileNameWei2);
      
      /*
      if (NRwFile->IsZombie() || NRwFile2->IsZombie() ){
	cout<<" Input file does not exist!"<<endl;
	exit(1);
      }
      NRwFile->Print();
      NRwFile2->Print();

      TList *histList = NRwFile->GetListOfKeys();
      for (UInt_t n=0; n<1507; n++){
	if (histList->Contains(Form("point_%i_weights",n)))
	  NR_Wei_Hists[n] = (TH2F*)NRwFile->Get(Form("point_%i_weights",n));
	else
	  cout<<"This one does not existe pas: "<<n<<endl;
      }
      

      TList *histList2 = NRwFile2->GetListOfKeys();
      for (UInt_t n=0; n<12; n++){
	if (histList2->Contains(Form("point_%i_weights",n)))
	  NR_Wei_Hists[1507+n] = (TH2F*)NRwFile2->Get(Form("point_%i_weights",n));
	else
	  cout<<"This one does not existe pas: "<<1507+n<<endl;
      }
      */
    }
  }
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

  hists_["DummyHisto"]->Fill(123);
  globVar_->fill(event);

  if (!isRealData){
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    event.getByLabel(edm::InputTag("generator"), genEvtInfo );
    ww = genEvtInfo->weight();

    //if (totEvents==0)  W0 = fabs(ww); // This is the Etalon!
    //if (fabs(ww)!=W0) throw cms::Exception("ETALON","The weights are not the same. Can not deal with this... ww=")<<ww<<"  W0="<<W0;

    //ww/=W0; //Divide by the etalon, get +/-1
  }

  genTotalWeight = ww;

  if (ww<0) count_neg++;
  else count_pos++;

  totEvents++;
  totWeights+=ww;

  CountEvents(0, "Ntuple events",ww,fcuts);
  FillHistoCounts(0, ww);

  eventNumber = event.id().event();
  runNumber = event.id().run();

  o_run = runNumber;
  o_evt = eventNumber;

  edm::Handle<double> rhoHandle;
  event.getByLabel( rhoFixedGrid_, rhoHandle );
  const double rhoFixedGrd = *( rhoHandle.product() );
  //const double rhoFixedGrd = globVar_->valueOf(globVar_->indexOf("rho"));


  FHM->Reset(rhoFixedGrd, 1, eventNumber);
  tools->setRho(rhoFixedGrd);

  isSignal   = 1;
  isPhotonCR = 0;
  CosThetaStar = -999;
  nvtx = globVar_->valueOf("nvtx");
  //if (nvtx!=nvtx2) cout<<"\t nvtx = "<< nvtx<<" nvtx2 = "<<nvtx2<<endl;;

  diphotonCandidate.SetPxPyPzE(0,0,0,0);
  leadingPhoton.SetPxPyPzE(0,0,0,0);
  subleadingPhoton.SetPxPyPzE(0,0,0,0);
  leadingJet.SetPxPyPzE(0,0,0,0);
  leadingJet_bDis = 0;
  leadingJet_flavour = 0;
  leadingJet_hadFlavour = 0;
  subleadingJet.SetPxPyPzE(0,0,0,0);
  subleadingJet_bDis = 0;
  subleadingJet_flavour = 0;
  subleadingJet_hadFlavour = 0;
  dijetCandidate.SetPxPyPzE(0,0,0,0);
  diHiggsCandidate.SetPxPyPzE(0,0,0,0);

  gen_mHH = 0;
  gen_cosTheta = -99;

  std::vector<TLorentzVector> gen_photons, gen_jets;
  //TLorentzVector gen_gamma1, gen_gamma2;

  if (!isRealData) {
    //&& sample = signal...

    edm::Handle<vector<reco::GenParticle> > genParts;
    event.getByLabel( myGen_, genParts );

    //TLorentzVector gen_bjet1, gen_bjet2;
    TLorentzVector gen_bQ1, gen_bQ2;
    TLorentzVector tmp;

    for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {

      if (abs(igen->pdgId())==5 && igen->mother() && igen->mother()->pdgId()==25) {

	//std::cout<<totEvents<<"\t"<<igen->pdgId()<<"\t\t BBB Found b-quark from Higgs! BBB  Its status = "<<igen->status()
	//<<" Its charge = "<<igen->charge()<<std::endl;
	tmp.SetPxPyPzE(igen->px(), igen->py(), igen->pz(), igen->energy());
	if (igen->pdgId()==5)
	  gen_bQ1 = tmp;
	if (igen->pdgId()==-5)
	  gen_bQ2 = tmp;
      }

      //if (igen->pdgId()==22)
      //std::cout<<totEvents<<"\t"<<igen->pdgId()<<"\t\t Found a photon. Its status = "<<igen->status()
      //<<" isPrompt = "<<igen->isPrompt()
      //<<" isPrompt FS = "<<igen->isPromptFinalState()
      //<<" from HP FS = "<<igen->fromHardProcessFinalState()<<std::endl;
      if (igen->pdgId()==22) {
	//if (igen->pdgId()==22 && igen->isPromptFinalState()) {
	tmp.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());
	gen_photons.push_back(tmp);

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

    if (nodesOfHH) {

      // ----
      // Save info for HH nodes samples, for re-weighting etc.
      // ---

      TLorentzVector H1, H2;
      UInt_t nH = 0;

      for( vector<reco::GenParticle>::const_iterator igen = genParts->begin(); igen != genParts->end(); ++igen ) {

	if (igen->pdgId()==25 && igen->isHardProcess()){

	  //std::cout<<eventNumber<<"  Higgs it is!   Pt="<<igen->pt()<<"  M="<<igen->mass()
	  //<<"\n \t Is Hard ="<<igen->isHardProcess()<<" daug-ter = "<<igen->daughter(0)->pdgId()<<endl;
	  if (nH==0)
	    H1.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());
	  if (nH==1)
	    H2.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());

	  nH++;
	  if (nH==2) break;

	  /*
	    if (igen->daughter(0)->pdgId()==22)
	    H1.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());

	    if (abs(igen->daughter(0)->pdgId())==5)
	    H2.SetXYZM(igen->px(), igen->py(), igen->pz(), igen->mass());
	  */
	}
      }

      if (nH==2){

	gen_ptH1 = H1.Pt();
	gen_ptH2 = H2.Pt();
	gen_mHH  = (H1+H2).M();
	gen_cosTheta = angles->getCosThetaStar_CS(H1,H2,6500);

	// Another version of costheta star:
	TLorentzVector P1boost = H1; // take one higgs
	TLorentzVector P12 = H1 + H2; // this is the total vectorial momenta of the system
	P1boost.Boost(-P12.BoostVector());
	gen_cosTheta2 = P1boost.CosTheta(); // this is the costTheta

	hists->fill1DHist(fabs(gen_cosTheta2)-fabs(gen_cosTheta), "diffCosThetaStarrr",
			  ";|#cos{#theta*}^{1}| - |#cos{#theta*}^{2}|", 100,-0.01,0.01, 1, "GEN");
	hists->fill1DHist(H1.M(), "gen_mH1", "m(H_{1})", 100,120,130, 1, "GEN");
	hists->fill1DHist(H2.M(), "gen_mH2", "m(H_{2})", 100,120,130, 1, "GEN");
	hists->fill1DHist(P12.Pt(), "gen_HH_PT", "p_{T}(HH)", 100,0,400, 1, "GEN");

	if (doNonResWeights_ && nodesOfHH){
	  for (UInt_t n=0; n<1519; n++){

	    NRWeights[n] = NRW->GetWeight(n, gen_mHH, gen_cosTheta); 

	    /*
	    if (n==324 || n==910 || n==985 || n==990){
	      // The points above do not exist in the input file provided by Alexandra (and wont ever be added)

	      //cout<<"This one was not existing in the input file: "<<n<<endl;
	      NRWeights[n]=0;
	    }
	    else {
	      UInt_t binNum = NR_Wei_Hists[n]->FindBin(gen_mHH, fabs(gen_cosTheta));
	      NRWeights[n]  = NR_Wei_Hists[n]->GetBinContent(binNum);
	      // Just print out for one n:
	      //if (n==100 || n==1510) 
	      //	cout<<n<<" **  mHH = "<<gen_mHH<<"   cosT*="<<fabs(gen_cosTheta)
	      //    <<"  bin="<<binNum<<" wei="<<NRWeights[n]<<endl;
	    }

	    */
	  }
	}
	genTree->Fill();
      } // if nH==2
      else {
	throw cms::Exception("Non Res Nodes")<<"Number of Higges found is not equal 2; nH = "<<nH;
      }

    }  //End of NodeOfHH
  }  // End of isRealData

  // ---------
  // Verteces
  // ----------
  edm::Handle<reco::VertexCollection> primaryVtcs;
  event.getByLabel(edm::InputTag("offlineSlimmedPrimaryVertices"), primaryVtcs);
  if (primaryVtcs->size()<1) return;
  const edm::Ptr<reco::Vertex> CandVtx(primaryVtcs, 0);

  nvtx2=0;
  for(reco::VertexCollection::const_iterator iVtx = primaryVtcs->begin(); iVtx != primaryVtcs->end(); ++iVtx){
    reco::Vertex myVtx = reco::Vertex(*iVtx);
    if(!myVtx.isValid() || myVtx.isFake()) continue;
    if (fabs(myVtx.z()) > 24) continue;
    if (myVtx.ndof()<4) continue;
    if (myVtx.position().rho()>2) continue;
    ++nvtx2;
  }


  // ---------
  // TRIGGER
  //---------

  if(myTriggers_.size() > 0){
    edm::Handle<edm::TriggerResults> trigResults;
    event.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), trigResults);
    const edm::TriggerNames &names = event.triggerNames(*trigResults);

    myTriggerResults = tools->TriggerSelection(myTriggers_, names, trigResults);
  }


  CountEvents(1, "HLT di-photon Triggers",ww,fcuts);
  FillHistoCounts(1, ww);

  for(Int_t ev=0; ev<evSize;ev++){
    if (eventNumber==hisEVTS[ev])
      {cout<<"First  ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	break;}}

  // --------
  // Di-Photon selection
  // ---------
  vector<TLorentzVector> myPhotons;

  edm::Handle<std::vector<flashgg::DiPhotonCandidate> > diPhotons;
  event.getByLabel(diPhotons_, diPhotons);

  //std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> SeldiPhotons = tools->DiPhoton76XPreselection(diPhotons, myTriggerResults);

  // This is to store the selected di-photons:
  std::vector<flashgg::DiPhotonCandidate> myDiPho;

  Bool_t kinema = false;
  if (useDiPhotons_) {
    for(std::vector<flashgg::DiPhotonCandidate>::const_iterator it=diPhotons->begin(); it!=diPhotons->end(); ++it){

      if (!tools->passHgg76XPreselection(&(*it), myTriggerResults)) continue;

      const flashgg::Photon *p1, *p2;
      p1 = it->leadingPhoton();
      p2 = it->subLeadingPhoton();
      /*
      for(Int_t ev=0; ev<evSize;ev++){
	if (eventNumber==hisEVTS[ev])
	  {cout<<"\t Loop over di-photons ---> "<<eventNumber<<" <--- "<<endl;
	    cout<<"mgg = "<<it->mass()
		<<" pt1 ="<<p1->pt()<<" eta="<<p1->eta()<<" mva="<<p1->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values")
		<<" pt2 ="<<p2->pt()<<" eta="<<p2->eta()<<" mva="<<p2->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values")<<endl;
	      
	    break;}}
      */
      
      Float_t tmp_mgg = it->mass();
      if (p1->pt() > 0.333*tmp_mgg && p2->pt() > tmp_mgg/4 &&
	  tmp_mgg > 100 && tmp_mgg < 180)
	//&&	  deltaR(*p1,*p2) > 0.4)
	kinema = true;
      else continue;
	//if (p2->pt() < 15) continue;

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

	passISO = ((p1->isEB() && tools->isPhoISO(&(*it), 0, phoISOcutEB_, phoISO_nhCorEB_, phoISO_phCorEB_)) ||
		   (p1->isEE() && tools->isPhoISO(&(*it), 0, phoISOcutEE_, phoISO_nhCorEE_, phoISO_phCorEE_))
		   );
	// SubLeading photon:
        passID = passID && (( p2->isEB() && tools->isPhoID(&(*p2), phoIDcutEB_)) ||
			    ( p2->isEE() && tools->isPhoID(&(*p2), phoIDcutEE_))
			    );
	passISO = passISO && ((p2->isEB() && tools->isPhoISO(&(*it), 1, phoISOcutEB_, phoISO_nhCorEB_, phoISO_phCorEB_)) ||
			      (p2->isEE() && tools->isPhoISO(&(*it), 1, phoISOcutEE_, phoISO_nhCorEE_, phoISO_phCorEE_))
			      );

	passID = passID && passISO;
	break;

      case 2 :
	// Cut based ID stored in AOD/PAT/fgg
	// Warning: it is for 50ns! Need to switch to 25ns, once available
        passID = (p1->photonID("cutBasedPhotonID-Spring15-50ns-V1-standalone-medium") &&
		  p2->photonID("cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"));
        break;

      case 3 :
	// EGamma MVA ID stored in the PAT object

	// Cut on the value
	passID = (( p1->isEB() && p1->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values") > 0.374 ) ||
		  ( p1->isEE() && p1->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values") > 0.336 ));

	passID = (passID &&
		  (( p2->isEB() && p2->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values") > 0.374 ) ||
		   ( p2->isEE() && p2->userFloat("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values") > 0.336 )));

	// Use the boolean
	//passID = (p1->photonID("mvaPhoID-Spring15-25ns-nonTrig-V2p1-wp90") &&
	//	  p2->photonID("mvaPhoID-Spring15-25ns-nonTrig-V2p1-wp90"));
	break;

      case 4 :
        // Hgg MVA ID:
        //https://github.com/cms-analysis/flashgg/blob/f1ad5f8c6d9b02481656c8beec8fc5e10ffa35b7/DataFormats/interface/Photon.h#L136
        //// if lazy flag is true only compare key (needed since fwlite does not fill provenance info)
        passID = ((p1->phoIdMvaDWrtVtx(CandVtx, true) > 0.2) &&
		  (p2->phoIdMvaDWrtVtx(CandVtx, true) > 0.2));

        break;

      case 5 :
	// High Pt ID: TBD
        break;

      default : break;

      }

      passID = passID && p1->passElectronVeto() && p2->passElectronVeto();

      // Z+Jets CR: reverted electron vetoes
      //passID = passID && !p1->passElectronVeto() && !p2->passElectronVeto();


      if (!passID) continue;

      myDiPho.push_back(*it);

    }
  }
  else { // Use regular photons

    edm::Handle<std::vector<flashgg::Photon> > photons;
    event.getByLabel(photons_, photons);

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

  // -----------
  // SET PHOTONS
  //------------

  sort(myDiPho.begin(), myDiPho.end(), fggSortCondition);

  TLorentzVector gamma1, gamma2;
  // Bools for removing double counting photons from the bkg MC samples
  Bool_t isPrompt1=false, isPrompt2=false;
    
  if (useDiPhotons_){
    if (myDiPho.size()<1) return;

    CountEvents(2, "Di-Photon exists",ww,fcuts);
    FillHistoCounts(2, ww);

    if (!kinema) return;
    CountEvents(3, "Kin selection of di-photon",ww,fcuts);
    FillHistoCounts(3, ww);

    const flashgg::Photon *p1, *p2;
    p1 = myDiPho[0].leadingPhoton();
    p2 = myDiPho[0].subLeadingPhoton();
    gamma1 = TLorentzVector(p1->px(), p1->py(), p1->pz(), p1->energy());
    gamma2 = TLorentzVector(p2->px(), p2->py(), p2->pz(), p2->energy());


    if (p1->genMatchType()==1) isPrompt1=1;
    if (p2->genMatchType()==1) isPrompt2=1;

    /* //checks and printouts
       if (runSample_=="DiPhoton" && (p1->genMatchType()!=1 || p2->genMatchType()!=1))
       {
       cout<<"\t lead pho match type = "<<p1->genMatchType()<<endl;
       cout<<"\t sub pho match type = "<<p2->genMatchType()<<endl;
       }

       if (runSample_=="GJets" && (p1->genMatchType()==1 && p2->genMatchType()==1))
       {
       cout<<"\t lead pho match type = "<<p1->genMatchType()<<endl;
       cout<<"\t sub pho match type = "<<p2->genMatchType()<<endl;
       }
    */
  }
  else{

    sort(myPhotons.begin(), myPhotons.end(), P4SortCondition);
    if (myPhotons.size()<2) return;
    gamma1 = myPhotons[0];
    gamma2 = myPhotons[1];


    CountEvents(2, "Two photons pass ID",ww,fcuts);
    FillHistoCounts(2, ww);

  }


  
  for(Int_t ev=0; ev<evSize;ev++){
    if (eventNumber==hisEVTS[ev])
      {cout<<"\t Second  ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	cout<<"DiPho information:"<<endl;
	cout<<"mgg ="<<myDiPho[0].mass()<<"  ptDiPho="<<myDiPho[0].pt()<<endl;
	cout<<"g1 Pt="<<gamma1.Pt()<<"  g2 Pt="<<gamma2.Pt()<<endl;
	break;}}
  
  // ---------
  //  JETS
  // --------

  //edm::Handle<reco::GenJetCollection> genJets;
  //event.getByLabel(edm::InputTag("slimmedGenJets"), genJets);

  vector<TLorentzVector> myJets, myJets25;
  vector<flashgg::Jet> bJets;  // These will be ordered by b-tag score
  flashgg::Jet bjet1, bjet2; // Those are the two highest scores of b-tag
  //TLorentzVector bjet1, bjet2;
  UInt_t ind1=99, ind2=99;     // Using these indexes to order two jets by Pt

  edm::Handle<std::vector<std::vector<flashgg::Jet>>> jetsCols;
  event.getByLabel(jets_, jetsCols);
  if (jetsCols->size()<1) return;

  UInt_t jetCollIndex = 0;
  if (useDiPhotons_)
    jetCollIndex = myDiPho[0].jetCollectionIndex();

  std::vector<flashgg::Jet> theJets = jetsCols->at(jetCollIndex);

  // -------------------
  // B-jet energy regression
  // ---------------
  if(doBJetRegression_) {
    // std::cout << "DOING REGRESSION! " << std::endl;

    // First of all we need the Missing Et as an input:
    edm::Handle<pat::METCollection> METs;
    event.getByLabel(edm::InputTag("slimmedMETs"), METs);

    if( METs->size() != 1 )
      { std::cout << "WARNING number of MET is not equal to 1" << std::endl; }

    //edm::Ptr<pat::MET> theMET = METs->ptrAt( 0 );
    const edm::Ptr<pat::MET> theMET(METs, 0);
    jetReg->RegressedJets(theJets, nvtx, theMET->p4() );

  }  // End of b-jet regression


  UInt_t nJets=0;
  for(UInt_t j = 0 ; j < theJets.size() ; j++ ) {
    nJets++;

    const flashgg::Jet *jet = &(theJets[j]);

    //flashgg::Jet jet = theJets[j];
    //const flashgg::Jet *jetRef = &(jet);

    if (!tools->isJetID( jet )) continue;
    //if (!jet->passesJetID(flashgg::JetIDLevel::Loose)) continue;

    //std::cout<<j<<" Pass jet ID. Pt="<<jet->pt()<<std::endl;
    //std::cout<<j<<" Does not Pass jet ID. Pt="<<jet->pt()<<std::endl;


    // does not work in flashgg:
    //if (!jet->passesPuJetId( need Vtx here)) continue;

    //std::cout<<jetsCol->at(0)[j].pt()<<"  from pat = "<<jet->pt()<<std::endl;


    TLorentzVector tmp = TLorentzVector(jet->px(), jet->py(), jet->pz(), jet->energy());

    /*
    for(Int_t ev=0; ev<evSize;ev++){
      if (eventNumber==hisEVTS[ev])
	{cout<<"\t Inside the Jets loop  ---> "<<eventNumber<<" <--- Found an event.   jetCollIndex="<<jetCollIndex<<endl;
	  cout<<j<<" Pt="<<jet->pt()<<"  eta="<<jet->eta()<<"  bTag="<<jet->bDiscriminator(bTagName)
	      <<" dR(g1) "<<tmp.DeltaR(gamma1)<< " dR(g2) "<<tmp.DeltaR(gamma2)<<endl;
	  break;}}
    */
    
    if (jet->pt() > 25){

      // TLorentzVector tmp = TLorentzVector(jet->px(), jet->py(), jet->pz(), jet->energy());
      myJets.push_back(tmp);

      if (fabs(jet->eta()) > 2.4) continue;
      if (tmp.DeltaR(gamma1) < 0.4 || tmp.DeltaR(gamma2) < 0.4) continue;

      //
      // Order the Jets by the b-discriminator value:
      //
      if (jet->bDiscriminator(bTagName) < 0) continue;

      myJets25.push_back(tmp);

      // Begin with putting the first jet in the array
      if (bJets.size()==0){
	bJets.push_back(*jet);
	continue;
      }

      // Now loop over all and order them by b-tag discriminator
      for (std::vector<flashgg::Jet>::const_iterator b = bJets.begin(); b!= bJets.end(); b++) {
	auto nx = std::next(b);
	if ( jet->bDiscriminator(bTagName) > b->bDiscriminator(bTagName)) {
	  bJets.insert(b, *jet);
 	  break;
	}
	else if (nx == bJets.end()){
	  bJets.push_back(*jet);
	  break;
	}
      }
      // END of b-jet ordering
    }
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

  // Set up gammas for Histogramming tool:
  FHM->SetGamma1(gamma1);
  FHM->SetGamma2(gamma2);
  Float_t Mgg = (gamma1 + gamma2).M();

  //if (totEvents>100) return;

  // ----------------- 
  // Two B-Jet selection
  // Getting first two b-jet indexes ordered by pT:
  // ---------- 

  UInt_t NbJets = bJets.size();
  Bool_t foundPair = false; 
  for (UInt_t j1=0; j1<NbJets; j1++){
    for (UInt_t j2=j1+1; j2<NbJets; j2++){
      Float_t mjj = (bJets[j1].p4() + bJets[j2].p4()).mass();
      if ( mjj > 80 && mjj < 200 ){
	if (bJets[j1].pt() > bJets[j2].pt()) {
	  ind1=j1; ind2=j2;}
	else {
	  ind1=j2; ind2=j1;}
	foundPair = true;
	break;
      } 
    }
    if (foundPair) break;
  }
  // -------------------- 

  Double_t Mbjbj = 0;
  Double_t Mtot  = 0;
  
  //if (isPrompt1 && isPrompt2 && Mtot && Mbjbj && ind1 && ind2) { /*Do nothing to make it compile*/}   


  // --------------------
  // ---- Final cut flow ----
  // ---------------------

  switch ( cutFlow_ ) {
  case 1 :
    // Default cutflow, for sync and final numbers


    /* //Double-counting removal on dR matching
       // Does not work since gen_photons are not prompt (isPromptFinalState method don't work)
      for (size_t g = 0;  g<gen_photons.size(); g++)
      {
      Float_t dR1 = gamma1.DeltaR(gen_photons[g]);
      Float_t dR2 = gamma2.DeltaR(gen_photons[g]);
      hists->fill1DHist(dR1, "GEN_dR1_gamma-prompt",
      ";#DeltaR(gen prompt #gamma, reco #gamma)", 100,0,6,  ww, "GEN");
      hists->fill1DHist(dR2, "GEN_dR2_gamma-prompt",
      ";#DeltaR(gen prompt #gamma, reco #gamma)", 100,0,6,  ww, "GEN");

      if (dR1 < 0.3) isPrompt1 = 1;
      if (dR2 < 0.3) isPrompt2 = 1;
      }
    */
    // Finished with double-counting removal


    //if (gamma1.DeltaR(gamma2) < 0.4) return;
    //if (gamma1.Pt() < 30 || gamma2.Pt() < 30) return;
    //if (fabs(gamma1.Eta()) > 2.5 || fabs(gamma2.Eta()) > 2.5) return;


    CountEvents(4, " Empty",ww,fcuts);
    FillHistoCounts(4, ww);


    if (myJets25.size()<2) return;
    if (NbJets<2) return;
    if (ind1==99 || ind2==99) return;

    bjet1 = bJets[ind1];
    bjet2 = bJets[ind2];
    
    FHM->SetBJet1(TLorentzVector(bjet1.px(), bjet1.py(), bjet1.pz(), bjet1.energy()));
    FHM->SetBJet2(TLorentzVector(bjet2.px(), bjet2.py(), bjet2.pz(), bjet2.energy()));
    //TLorentzVector tmp;
    //tmp.SetPxPyPzE(bjet1.px(), bjet1.py(), bjet1.pz(), bjet1.energy());
    //tmp.SetPxPyPzE(bjet2.px(), bjet2.py(), bjet2.pz(), bjet2.energy());
    //FHM->SetBJet2(tmp);
      
    CountEvents(5, "At least two Jets w/ pT>25 and |eta|<2.4",ww,fcuts);
    FillHistoCounts(5, ww);
    FHM->MakeMainHistos(5, ww);
    FHM->MakeNPlots(5, myPhotons.size(), myJets.size(), bJets.size(), ww);

    sort(myJets25.begin(), myJets25.end(), P4SortCondition);

    /*
    for(Int_t ev=0; ev<evSize;ev++){
      if (eventNumber==hisEVTS[ev])
      {cout<<"\t\t Third  ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	cout<<"\t DiPho information:"<<endl;
	cout<<"mgg ="<<myDiPho[0].mass()<<"  ptDiPho="<<myDiPho[0].pt()<<endl;
	cout<<"\t Jets information:"<<endl;
	cout<<"bJets.size ="<<bJets.size()<<"  myJets25.size = "<<myJets25.size()<<endl;
	cout<<"mjj ="<<Mbjbj<<"  pt1="<<bjet1.pt()<<"  pt2="<<bjet2.pt()<<endl;
	cout<<"btag1 = "<<bjet1.bDiscriminator(bTagName)<<" btag2="<<bjet2.bDiscriminator(bTagName)<<endl;
	//if (bJets.size()>2)
	//cout<<"btag3 = "<<bJets[2].bDiscriminator(bTagName)<<" pt3="<<bJets[2].pt()<<endl;

	break;}}
    */
    
    if (runSample_=="DiPhoton"   && !(isPrompt1 & isPrompt2)) return;
    else if (runSample_=="GJets" && !(isPrompt1 ^ isPrompt2)) return;
    else if (runSample_=="QCD"   &&  (isPrompt1 & isPrompt2)) return;

    CountEvents(6, "After double counting removal",ww,fcuts);
    FillHistoCounts(6, ww);
    FHM->MakeMainHistos(6, ww);


    CountEvents(7, "Empty",ww,fcuts);
    FillHistoCounts(7, ww);
    FHM->MakeMainHistos(7, ww);


    FHM->MakeJetPlots(bjet1, "bJet-Lead");
    FHM->MakeJetPlots(bjet2, "bJet-Sub");
    hists->fill1DHist(bjet1.bDiscriminator(bTagName) + bjet2.bDiscriminator(bTagName),
		      "bJets_sumOfBtags",";Sum of two bTag Discriminators", 100, 0, 2.3, ww, "bJets");

    if (bjet1.pt() < 25 || bjet2.pt() < 25) return;
    if (fabs(bjet1.eta()) > 2.4 || fabs(bjet2.eta()) > 2.4) return;
    
    CountEvents(8, "The 2 Jets pT > 25 GeV and |eta|<2.4",ww,fcuts);
    FillHistoCounts(8, ww);
    FHM->MakeMainHistos(8, ww);
    FHM->MakeNPlots(8, myPhotons.size(), myJets.size(), bJets.size(), ww);
    

    /*
    for(Int_t ev=0; ev<evSize;ev++){
      if (eventNumber==hisEVTS[ev])
	{cout<<"\t\t\t Fourth  ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	cout<<"Jets information:"<<endl;
	cout<<"mjj ="<<Mbjbj<<"  pt1="<<bjet1.pt()<<"  pt2="<<bjet2.pt()<<endl;
	break;}}
    */

    CountEvents(9, "... Reserved",ww,fcuts);
    FillHistoCounts(9, ww);
    FHM->MakeMainHistos(9, ww);

    /*
    for(Int_t ev=0; ev<evSize;ev++){
      if (eventNumber==hisEVTS[ev])
      {cout<<"\t\t\t\t Fifth  ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	cout<<"Jets information:"<<endl;
	cout<<"mjj ="<<Mbjbj<<"  pt1="<<bjet1.pt()<<"  pt2="<<bjet2.pt()<<endl;
	break;}}

    */

    // -----
    // Filling all the objects for flatTree
    // ----------
    o_weight = ww;
    o_bbMass = Mbjbj;
    o_ggMass = Mgg;
    o_bbggMass = Mtot;

    leadingPhoton     = LorentzToLorentz(gamma1);
    subleadingPhoton  = LorentzToLorentz(gamma2);
    diphotonCandidate = LorentzToLorentz(gamma1+gamma2);

    leadingJet = bJets[ind1].p4();
    leadingJet_bDis = bJets[ind1].bDiscriminator(bTagName);
    leadingJet_flavour = bJets[ind1].partonFlavour();
    leadingJet_hadFlavour = bJets[ind1].hadronFlavour();
    
    subleadingJet = bJets[ind2].p4();
    subleadingJet_bDis = bJets[ind2].bDiscriminator(bTagName);
    subleadingJet_flavour = bJets[ind2].partonFlavour();
    subleadingJet_hadFlavour = bJets[ind2].hadronFlavour();
  
    dijetCandidate   = leadingJet + subleadingJet;

    diHiggsCandidate = diphotonCandidate + dijetCandidate;
    
    Mbjbj = dijetCandidate.mass();
    Mtot  = diHiggsCandidate.mass();
    


    if ( Mbjbj < 80 || Mbjbj > 200) return;

    flatTree->Fill();
    // End of Filling.
    // All flatTree objects must be filled now

    CountEvents(10, "80 < m(jj) < 200 GeV",ww,fcuts);
    FillHistoCounts(10, ww);
    FHM->MakeMainHistos(10, ww);


    for(Int_t ev=0; ev<evSize;ev++){
      if (eventNumber==hisEVTS[ev])
	{cout<<"\t\t\t\t\t Sixth and Final  ---> "<<eventNumber<<" <--- Found an event after cut "<<endl;
	  cout<<"Jets information:"<<endl;
	  cout<<"mjj ="<<Mbjbj<<"  pt1="<<bjet1.pt()<<"  pt2="<<bjet2.pt()<<endl;
	  break;}}

    // Signal Region: High Purity
    if (bjet1.bDiscriminator(bTagName) > 0.8 && bjet2.bDiscriminator(bTagName) > 0.8){
      CountEvents(11, "SR: High Purity",ww,fcuts);
      FillHistoCounts(11, ww);
      FHM->MakeMainHistos(11, ww);

      o_category = 0;
    }


    // Signal Region: Medium Purity
    if (bjet1.bDiscriminator(bTagName) > 0.8 && bjet2.bDiscriminator(bTagName) < 0.8){
      CountEvents(12, "SR: Medium Purity",ww,fcuts);
      FillHistoCounts(12, ww);
      FHM->MakeMainHistos(12, ww);

      o_category = 1;
    }

    // Control Region: Light Jets
    if (bjet1.bDiscriminator(bTagName) < 0.8 && bjet2.bDiscriminator(bTagName) < 0.8){
      CountEvents(13, "CR: Light Jets",ww,fcuts);
      FillHistoCounts(13, ww);
      FHM->MakeMainHistos(13, ww);
    }
    
    break;

  case 2 :
    // Alternative  cutflow for checks and stuff

    if (myJets25.size()<2) return;
    CountEvents(3, "At least two Jets w/ pT>25 (no b-tag)",ww,fcuts);
    FillHistoCounts(3, ww);
    FHM->MakeNPlots(3, myPhotons.size(), myJets.size(), bJets.size(), ww);


    if (myJets25.size()<2) return;
    CountEvents(4, "At least two Jets w/ pT>25 and |eta|<2.4",ww,fcuts);
    FillHistoCounts(4, ww);
    //FHM->MakeMainHistos(4, ww);

    sort(myJets25.begin(), myJets25.end(), P4SortCondition);

    if (NbJets<2) return;
    if (ind1==99 || ind2==99) return;

    CountEvents(5, "Two jets with bDissc > 0",ww,fcuts);
    FillHistoCounts(5, ww);
    FHM->MakeMainHistos(5, ww);


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

   
    if (bjet1.pt() < 25 || bjet2.pt() < 25) return;
    if (fabs(bjet1.eta()) > 2.4 || fabs(bjet2.eta()) > 2.4) return;

    CountEvents(9, "The 2 Jets pT > 25 GeV and |eta|<2.4",ww,fcuts);
    FillHistoCounts(9, ww);
    FHM->MakeMainHistos(9, ww);

    FHM->MakeNPlots(9, myPhotons.size(), myJets.size(), bJets.size(), ww);

    
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


LorentzVector HHbbggAnalyzer::LorentzToLorentz(const TLorentzVector& v)
{
  LorentzVector LoL;
  LoL.SetPxPyPzE(v.Px(),v.Py(),v.Pz(),v.E());
  return LoL;
}

TLorentzVector HHbbggAnalyzer::LorentzToLorentz(const LorentzVector& v)
{
  TLorentzVector LoL;
  LoL.SetPxPyPzE(v.Px(),v.Py(),v.Pz(),v.E());
  return LoL;
}
