#include "../interface/FggHistMakerBase.h"

FggHistMakerBase::FggHistMakerBase(HistManager *h)
{
  hists = h;
}

FggHistMakerBase::~FggHistMakerBase(){}


void FggHistMakerBase::MakeMuonPlots(const flashgg::Muon& mu, string dir)
{

  //hists->fill1DHist(mu.PixelLayersWithMeasurement(),dir+"_mu_PixelLayersWithMeasurement",";PixelLayersWithMeasurement",  15,0,15, 1, dir);
  //hists->fill1DHist(mu.TrackLayersWithMeasurement(),dir+"_mu_TrackLayersWithMeasurement",";TrackerLayersWithMeasurement",30,0,30, 1, dir);
  //hists->fill1DHist(mu.NumberOfMatchedStations(),   dir+"_mu_NumberOfMatchedStations", ";NumberOfMatchedStations", 10, 0,10, 1, dir);
  //hists->fill1DHist(mu.numberOfValidHits(),     dir+"_mu_NumberOfValidHits",   ";NumberOfValidHits",   60, 0,60, 1, dir);
  //hists->fill1DHist(mu.NumberOfValidTrackerHits(),  dir+"_mu_NumberOfValidTrackerHits",";NumberOfValidTrackerHits",40, 0,40, 1, dir);
  //hists->fill1DHist(mu.NumberOfValidPixelHits(),    dir+"_mu_NumberOfValidPixelHits",  ";NumberOfValidPixelHits",  15, 0,15, 1, dir);
  //hists->fill1DHist(mu.NormalizedChi2_tracker(),    dir+"_mu_NormalizedChi2_tracker",  ";NormalizedChi2_tracker", 100, 0,4,  1, dir);
  //hists->fill1DHist(mu.normChi2(),            dir+"_mu_NormalizedChi2",          ";NormalizedChi2",         100, 0,4,  1, dir);
  //hists->fill1DHist(mu.Dxy(&_pv), dir+"_mu_dxy", ";dxy", 50, -0.02,0.02, 1, dir);
  //hists->fill1DHist(mu.Dz(&_pv),  dir+"_mu_dxz", ";dz",  50, -0.1,  0.1, 1, dir);
  //hists->fill1DHist(mu.PtError()/mu.Pt(), dir+"_mu_ptErrorOverPt", ";ptErrorOverPt", 50, 0,0.1, 1, dir);
  
  //cout<<"Muon plot making"<<endl;



  //Float_t muIso = ObjID->CalculateMuonIso(&mu);
  //hists->fill1DHist(muIso, dir+"_mu_iso", ";rel isolation, pu corr", 50, 0,0.7, 1, dir);
  //hists->fill2DHist(muIso, global_Mll, dir+"_mu_iso_vs_Mll", ";pfIsolationR04;M(l_{1},l_{2})", 50, 0,0.7, 50, 0,20, 1, dir);
  //Float_t iso2 = (mu.IsoMap("pfChargedHadronPt_R04") + mu.IsoMap("pfNeutralHadronEt_R04") + mu.IsoMap("pfPhotonEt_R04"))/mu.Pt();
  //hists->fill1DHist(iso2, dir+"_mu_iso_official, ";pfIsolationR04", 50, 0,0.7, 1, dir);

}

template <class EGammaObj>
void FggHistMakerBase::MakeEGammaCommonPlots(const EGammaObj& egm, TString n)
{
  const string dir = n.Data();

  
  hists->fill1DHist(egm.superCluster()->eta(),           dir+"-egm_SCEta",        ";SCEta",        100, -2.5, 2.5,  1,dir);
  hists->fill1DHist(egm.r9(),              dir+"-egm_R9",           ";R9",           100,    0.2, 1,  1,dir);
  /*
  hists->fill1DHist(egm.E1x3()/egm.E5x5(), dir+"-egm_E1x3OverE5x5", ";E1x3OverE5x5", 100,    0.1, 1,  1,dir);
  hists->fill1DHist(egm.E2x2()/egm.E5x5(), dir+"-egm_E2x2OverE5x5", ";E2x2OverE5x5", 100,    0.1, 1,  1,dir);
  hists->fill1DHist(egm.E1x5()/egm.E5x5(), dir+"-egm_E1x5OverE5x5", ";E1x5OverE5x5", 100,    0.1, 1,  1,dir);
  hists->fill1DHist(egm.E2x5()/egm.E5x5(), dir+"-egm_E2x5OverE5x5", ";E2x5OverE5x5", 100,    0.5, 1,  1,dir);
  //hists->fill1DHist(egm.E2x5Max()/egm.R9(),dir+"-egm_E2x5MaxOverR9",";E2x5MaxOverR9",100,    0, 1,  1,dir);

  hists->fill1DHist(egm.SigmaIEtaIPhi(),   dir+"-egm_SigmaIEtaIPhi",";SigmaIEtaIPhi", 100,-2e-4,2e-4,  1,dir);
  */
  hists->fill1DHist(egm.sigmaIetaIeta(),   dir+"-egm_SigmaIEtaIEta",";SigmaIEtaIEta", 100, 0, 0.016,   1,dir);
  //hists->fill1DHist(egm.sigmaIPhiIPhi(),   dir+"-egm_SigmaIPhiIPhi",";SigmaIPhiIPhi", 100, 0, 0.04,    1,dir);
  /*
  hists->fill1DHist(egm.SCEnergy(),        dir+"-egm_SCEnergy",     ";SCEnergy",      100, 10, 120,    1,dir);
  hists->fill1DHist(egm.EnergyRegression(),dir+"-egm_EnergyRegress",";Energy Regress",100, 10, 120,    1,dir);
  //  hists->fill1DHist(egm.GetNCrystals(),    dir+"-egm_nCrystals",    ";nCrystals",    160, 0, 160,     1,dir);

  if (!n.Contains("BaseSC")) {
    hists->fill1DHist(egm.HadOverEm(),     dir+"-egm_HadOverEm",    ";HadOverEm",    100,0.001, 0.05, 1,dir);
    hists->fill1DHist(egm.SCEtaWidth(),    dir+"-egm_SCEtaWidth",   ";SCEtaWidth",   100, 0,  0.03,   1,dir);
    hists->fill1DHist(egm.SCPhiWidth(),    dir+"-egm_SCPhiWidth",   ";SCPhiWidth",   100, 0,  0.12,   1,dir);
    hists->fill1DHist(egm.SCRawEnergy(),   dir+"-egm_SCRawEnergy",  ";SCRawEnergy",  100, 10, 120,    1,dir);
    if (fabs(egm.SCEta()) > 1.56)
      hists->fill1DHist(egm.PreShowerOverRaw(),dir+"-egm_PreShowerOverRaw",";PreShowerOverRaw",100,0,0.2, 1,dir);
  }

  if (n.Contains("Ele") && !n.Contains("BaseSC")){ //Not filled for photons
    hists->fill1DHist(egm.SCDeltaPhi(),   dir+"-egm_SCDeltaPhi",   ";SCDeltaPhiAtVtx", 100, -0.1, 0.1, 1, dir);
    hists->fill1DHist(egm.SCDeltaEta(),   dir+"-egm_SCDeltaEta",   ";SCDeltaEtaAtVtx", 100,-0.01,0.01, 1, dir);
  }
  */
  //if (n.Contains("BaseSC")){
  // egm.Dump();
  //}

}


void FggHistMakerBase::MakePhotonPlots(const flashgg::Photon& ph, string dir)
{
  //cout<<"DBG make photon plots  "<<dir<<endl;

  MakeEGammaCommonPlots(ph, dir+"-EGamma");

  hists->fill1DHist(ph.passElectronVeto(),      dir+"ph_trackVeto",     ";track veto",      3, 0, 3,   1,dir);
  //hists->fill1DHist(ph.ConversionVeto(), dir+"ph_ConversionVeto",";ConversionVeto",  3, 0, 3,   1,dir);
  hists->fill1DHist(ph.pfPhoIso04(),       dir+"ph_PfIsoPhoton",   ";PfIsoPhoton",   100, 0.01,5, 1,dir);
  hists->fill1DHist(ph.pfChgIsoWrtChosenVtx02(), dir+"ph_PfIsoCharged",  ";PfIsoCharged",  100, 0.01,2, 1,dir);
  //if (fabs(ph.SCEta())>1.6){ //These are only for EE
  if (fabs(ph.eta())>1.6){ //These are only for EE
    hists->fill1DHist(ph.esEffSigmaRR(),dir+"ph_ESEffSigmaRR_x",";ESEffSigmaRR_x",100, -5e-8,5e-8, 1,dir);
    //hists->fill1DHist(ph.ESEffSigmaRR()[1],dir+"ph_ESEffSigmaRR_y",";ESEffSigmaRR_y",100, -5e-8,5e-8, 1,dir);
  }
  //hists->fill1DHist(ph.IdMap("mvaScore"),    dir+"ph_mvaScore",      ";MVA score",     100, -1,1,   1,dir);
  //hists->fill1DHist(ph.CiCPF4chgpfIso02()[0],dir+"ph_CiCPF4chgpfIso02", ";ph_CiCPF4chgpfIso02", 100, 0,5, 1,dir);

  /*
  if (ph.R9() > 0.9){
    hists->fill1DHist(ph.IdMap("HadIso_R03")-0.005*ph.Et(), dir+"ph_HadIso_R03_R9high", ";ph_HadIso_R03_R9high", 100, 0, 60, 1,dir);
    hists->fill1DHist(ph.IdMap("TrkIso_R03")-0.002*ph.Et(), dir+"ph_TrkIso_R03_R9high", ";ph_HadIso_R03_R9high", 100, 0, 60, 1,dir);
  }
  else{
    hists->fill1DHist(ph.IdMap("HadIso_R03")-0.005*ph.Et(), dir+"ph_HadIso_R03_R9low", ";ph_HadIso_R03_R9low", 100, 0, 6, 1,dir);
    hists->fill1DHist(ph.IdMap("TrkIso_R03")-0.002*ph.Et(), dir+"ph_TrkIso_R03_R9low", ";ph_HadIso_R03_R9low", 100, 0, 6, 1,dir);
  }
  */
  //cout<<"Photon plot making"<<endl;

  //hists->fill1DHist(ph.(), dir+"ph_", ";", 100, 0, 100,1, dir);
}


void FggHistMakerBase::MakeElectronPlots(const flashgg::Electron& el, string dir)
{

  MakeEGammaCommonPlots(el, dir+"-EGamma");

  //hists->fill1DHist(el.IdMap("mvaScore"),dir+"el_mvaScore",";Dalitz MVA score",100, -1,1, 1,dir);
  //hists->fill1DHist(el.MvaID(), dir+"_el_mvaID",";HZZ mva ID", 100, -1,1, 1,dir);
  //hists->fill1DHist(el.FBrem(), dir+"_el_fbrem",";fbrem",      100, 0, 1, 1,dir);

  //hists->fill1DHist(el.InverseEnergyMomentumDiff(), dir+"_el_fabsEPDiff",";|1/E - 1/p|",100, 0, 0.06, 1, dir);
  //hists->fill1DHist(el.PtError()/el.Pt(),     dir+"_el_ptErrorOverPt",";ptErrorOverPt", 100,    0, 1, 1,dir);
  //hists->fill1DHist(el.NormalizedChi2(),      dir+"_el_gsfChi2", ";gsfChi2", 100, 0, 3, 1,dir);
  //hists->fill1DHist(el.NormalizedChi2Kf(),    dir+"_el_kfChi2",  ";kfChi2",  100, 0, 3, 1,dir);
  //hists->fill1DHist(el.DeltaEtaSeedCluster(), dir+"_el_deltaEtaSeedAtCalo", ";DeltaEtaSeedAtCalo", 100, -0.02, 0.02, 1,dir);
  //hists->fill1DHist(el.DeltaPhiSeedCluster(), dir+"_el_deltaPhiSeedAtCalo", ";DeltaPhiSeedAtCalo", 100, -0.06, 0.06, 1,dir);
  //hists->fill1DHist(1 - el.E1x5()/el.E5x5(),  dir+"_el_ome1x5oe5x5",";1 - e1x5/e5x5",100, 0, 1, 1,dir);
  //hists->fill1DHist(el.EoP(),    dir+"_el_EoP",    ";EoP",    100, 0, 6,   1,dir);
  //hists->fill1DHist(el.EoPout(), dir+"_el_EoPout", ";EoPout", 100, 0, 6,   1,dir);
  hists->fill1DHist(el.ip3d(),   dir+"_el_ip3d",   ";ip3d",   100, 0,0.025,1,dir);
  //hists->fill1DHist(el.IP3dSig(),dir+"_el_ip3dSig",";ip3dSig",100, 0, 3,   1,dir);
  //hists->fill1DHist(el.NumberOfValidHits(), dir+"_el_NumberOfValidHits", ";NumberOfValidHits", 25, 0,25, 1,dir);
  //hists->fill1DHist(el.TrackerLayersWithMeasurement(), dir+"_el_TrackerLayersWithMeasurement",
  //		    ";TrackerLayersWithMeasurement", 20, 0,20,1,dir);

  hists->fill1DHist(el.passConversionVeto(), dir+"_el_conv_passConv",";pass conv veto",  3, 0, 3, 1, dir);
  //hists->fill1DHist(el.ConversionMissHits(), dir+"_el_conv_MissHits",";conv miss hits",  5, 0, 5, 1, dir);

  /*
  hists->fill1DHist(el.ConversionDist(),   dir+"_el_convDist",  ";conv Dist",  100, 0,0.1,  1, dir);
  hists->fill1DHist(el.ConversionDcot(),   dir+"_el_convDcot",  ";conv Dcot",100,-0.6, 0.6, 1, dir);
  hists->fill1DHist(el.ConversionRadius(), dir+"_el_convRadius",";conv Radius", 100, 0, 20, 1, dir);

  hists->fill1DHist(el.GetTracks().size(), dir+"_el_nTracks", ";nTracks", 10, 0,10, 1,dir);
  */

  //cout<<"Electron plot making"<<endl;

  /*
  vector<flashgg::Electron::Track> trk = el.GetTracks();
  sort(trk.begin(), trk.end(), P4SortCondition);

  if (trk.size()>=2){
    Float_t mll = (trk[0]+trk[1]).M();
    Float_t dR  = trk[0].DeltaR(trk[1]);
    Float_t ptr = trk[1].Pt()/trk[0].Pt();
    Float_t EoverPt = el.SCRawEnergy()/(trk[0]+trk[1]).Pt();
    hists->fill1DHist(EoverPt, dir+"_el_EoverPt_all", ";E_{SC}^{raw}/p_{T}^{gsf1+gsf2}", 100,0,2.5, 1,dir);
    hists->fill1DHist(ptr, dir+"_el_ptRatio", ";p_{T}^{gsf2}/p_{T}^{gsf1}", 50,0,1, 1,dir);
    hists->fill1DHist(dR,  dir+"_el_dR_all",  ";#Delta R(gsf1, gsf2)",   100,   0,2, 1,dir);
    hists->fill1DHist(dR,  dir+"_el_dR_low",  ";#Delta R(gsf1, gsf2)",   100, 0,0.1, 1,dir);
    hists->fill1DHist(mll, dir+"_el_mll_full",";m_{ll} (GeV)",  100, 0,120, 1,dir);
    hists->fill1DHist(mll, dir+"_el_mll_low", ";m_{ll} (GeV)",  100,  0,20, 1,dir);
    hists->fill1DHist(mll, dir+"_el_mll_jpsi",";m_{ll} (GeV)",  100,  1, 5, 1,dir);
  }

  for (UInt_t i=0; i<trk.size() && i<3; i++){
    string subdir = dir+Form("-track-%i",i+1);
    hists->fill1DHist(trk[i].Pt(), subdir+Form("_el_01_track_%i",i+1)+"_Pt",
		      Form(";gsf track %i Pt (GeV)",i+1), 200, 0,100, 1,subdir);
    hists->fill1DHist(trk[i].PtError()/trk[i].Pt(), subdir+Form("_el_01_track_%i",i+1)+"_PtErrorOverPt",
		      Form(";Pt Error /Pt for track %i",i+1), 200, 0,1.5, 1,subdir);
    hists->fill1DHist(trk[i].NormalizedChi2(), subdir+Form("_el_track_%i",i+1)+"_01_NormalizedChi2",
		      Form(";NormalizedChi2 for track %i",i+1), 100, 0,3, 1,subdir);

    flashgg::Track::ConversionInfo inf = trk[i].GetConversionInfo();
    hists->fill1DHist(inf.isValid, subdir+Form("_el_conv_%i",i+1)+"_isValid",Form(";isValid trk %i", i+i),3, 0,3, 1,subdir);
    if (inf.isValid){
      hists->fill1DHist(inf.nHitsMax,subdir+Form("_el_conv_%i",i+1)+"_nHitsMax",Form(";nHitsMax trk %i",i+i),  20, 0,20, 1,subdir);
      hists->fill1DHist(inf.vtxProb, subdir+Form("_el_conv_%i",i+1)+"_vtxProb", Form(";vtxProb trk %i", i+i), 200,0,0.01,1,subdir);
      //hists->fill1DHist(inf.lxyPV,   subdir+Form("_el_conv_%i",i+1)+"_lxyPV",   Form(";lxyPV trk %i",   i+i), 200, 0,50, 1,subdir);
      hists->fill1DHist(inf.lxyBS,   subdir+Form("_el_conv_%i",i+1)+"_lxyBS",   Form(";lxyBS trk %i",   i+i), 200, 0,50, 1,subdir);
    }
  }

  vector<flashgg::EGamma> sc = el.BaseSC();
  sort(sc.begin(),  sc.end(),  SCESortCondition);
  hists->fill1DHist(sc.size(), dir+"_el_nBaseSC", ";nBaseSC", 10, 0,10, 1,dir);
  for (UInt_t i=0; i<sc.size() && i<2; i++)
    MakeEGammaCommonPlots(sc[i], dir+Form("-BaseSC-%i",i+1));

  */
}


/*
void FggHistMakerBase::MakePhotonEnergyCorrPlots(const flashgg::Photon& pho, Float_t corrPhoEnReco, Float_t corrPhoEnSC)
{
  const string dir = "Photon-EnergyCorr";
  hists->fill1DHist(pho.M(), "ph_recoMass",";M_{#gamma}", 100,-0.05,0.05, 1,dir);

  hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco",
		    ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

  hists->fill1DHist((pho.SCEnergy() - corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC",
		    ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,dir);

  hists->fill1DHist((pho.E() - pho.SCEnergy())/pho.E(),"ph_energyRecoVsSC",
		    ";(E_{#gamma}^{reco} - E_{#gamma}^{SC})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

  if (pho.Pt()>40 && fabs(pho.SCEta())<1.444){
    hists->fill1DHist((pho.E() - pho.SCEnergy())/pho.E(),"ph_energyRecoVsSC_PtEtaCut",
		      ";(E_{#gamma}^{reco} - E_{#gamma}^{SC})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

    hists->fill1DHist((pho.E() - corrPhoEnReco)/pho.E(),"ph_energyCorrectionReco_PtEtaCut",
		      ";(E_{#gamma}^{reco} - E_{#gamma}^{reco corr})/E_{#gamma}^{reco}", 100,-0.1,0.1, 1,dir);

    hists->fill1DHist((pho.SCEnergy() - corrPhoEnSC)/pho.SCEnergy(),"ph_energyCorrectionSC_PtEtaCut",
		      ";(E_{#gamma}^{SC} - E_{#gamma}^{SC corr})/E_{#gamma}^{SC}", 100,-0.1,0.1, 1,dir);

  }
}
*/

void FggHistMakerBase::MakeZeePlots(const flashgg::Photon& p1, const flashgg::Photon& p2)
{
  /*
  hists->fill2DHist(p1.R9(), p2.R9(), "zee_p1R9_p2R9",";#gamma_{1} R9; #gamma_{2} R9", 100, 0,1, 100,0,1,  1, "Zee");
  Float_t mZ = (p1+p2).M();
  hists->fill1DHist(mZ, "zee_M",";", 100, 70, 110, 1, "Zee");
  if (p1.R9() > 0.94 && p2.R9() > 0.94)
    hists->fill1DHist(mZ, "zee_M_highR9",";", 100, 70, 110, 1, "Zee");
  else
    hists->fill1DHist(mZ, "zee_M_lowR9",";", 100, 70, 110, 1, "Zee");

  //hists->fill1DHist(mZ, "zee_M",";", 200, 60, 130, 1, "Zee");
  //hists->fill1DHist(mZ, "zee_M",";", 200, 60, 130, 1, "Zee");
  */
}


void FggHistMakerBase::MakeNPlots(Int_t num, Int_t nmu, Int_t nele1, Int_t nele2, Int_t nphoHZG,
			   Int_t nphoTight, Int_t nphoMVA, Int_t nJets, Double_t w)
{
  hists->fill1DHist(nmu,       Form("size_mu_cut%i",  num), ";Number of muons",              5,0,5, w, "N");
  hists->fill1DHist(nele1,     Form("size_el_cut%i",  num), ";Number of electrons (HZZ)",    5,0,5, w, "N");
  hists->fill1DHist(nele2,     Form("size_el0_cut%i", num), ";Number of electrons (Loose)",  5,0,5, w, "N");
  hists->fill1DHist(nphoHZG,   Form("size_phHZG_cut%i",  num), ";Number of photons (HZG)",   5,0,5, w, "N");
  hists->fill1DHist(nphoTight, Form("size_phTight_cut%i",num), ";Number of photons (Tight)", 5,0,5, w, "N");
  hists->fill1DHist(nphoMVA,   Form("size_phMVA_cut%i",  num), ";Number of photons (MVA)",   5,0,5, w, "N");
  hists->fill1DHist(nJets,     Form("size_jets_cut%i",   num), ";Number of Jets",            5,0,5, w, "N");
}

float FggHistMakerBase::Zeppenfeld(const TLorentzVector& p, const TLorentzVector& pj1, const TLorentzVector& pj2)
{
  float zep = p.Eta()-(pj1.Eta()+pj2.Eta())/2.;
  return zep;
}
