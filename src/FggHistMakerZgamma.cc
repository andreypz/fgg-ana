#include "../interface/FggHistMakerZgamma.h"

FggHistMakerZgamma::FggHistMakerZgamma(HistManager *h):
  FggHistMakerBase::FggHistMakerBase(h),
  _isVtxSet(false),
  _isLepSet(false),
  _isGammaSet(false),
  _isGamma2Set(false),
  _isDaleSet(false),
  _isJet1Set(false),
  _isJet2Set(false),
  _isMetSet(false),
  _rhoFactor(-1),
  _nVtx(0)
{
  cout<<"\t HHHHHHH \t ZGammaHistmaker constractor"<<endl;
  angles = new Angles();
}

FggHistMakerZgamma::~FggHistMakerZgamma(){}

void FggHistMakerZgamma::Reset(float rho, UInt_t n){
  _isLepSet = 0; _isGammaSet=0; _isGamma2Set=0;
  _isDaleSet =0; _isJet1Set=0; _isJet2Set=0; _isMetSet=0;
  _rhoFactor = rho; _nVtx=n;}

void FggHistMakerZgamma::SetVtx(reco::Vertex v){
  _pv=v; _isVtxSet = 1;}
void FggHistMakerZgamma::SetLeptons(TLorentzVector l1, TLorentzVector l2){
  _lPt1 = l1;  _lPt2 = l2; _isLepSet = 1;}
void FggHistMakerZgamma::SetGamma(TLorentzVector g){
  _gamma = g; _isGammaSet=1;}
void FggHistMakerZgamma::SetGamma2(TLorentzVector g){
  _gamma2 = g; _isGamma2Set=1;}
void FggHistMakerZgamma::SetDalectron(TLorentzVector d){
  _dale = d; _isDaleSet=1;}
void FggHistMakerZgamma::SetJet1(TLorentzVector j){
  _jet1 = j; _isJet1Set=1;}
void FggHistMakerZgamma::SetJet2(TLorentzVector j){
  _jet2 = j; _isJet2Set=1;}
void FggHistMakerZgamma::SetMet(TVector2 m){
  _met = m; _isMetSet=1;}


void FggHistMakerZgamma::MakeMainHistos(Int_t num, Double_t weight, string dir)
{
  //cout<<num<<" Making plots "<<dir<<endl;
  char * d = new char [dir.length()+1];
  std::strcpy (d, dir.c_str());

  TLorentzVector diLep, tri, four;
  Float_t Mll(0), Mllg(0);

  if (_isLepSet){

    diLep = _lPt1 + _lPt2;
    Mll   = diLep.M();

    //hists->fill1DHist(Mll, Form("01_diLep_mass_0to20_%s_cut%i",    d, num),";m_{#mu#mu} (GeV)", 100, 0,20,  weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_0to20_b50_%s_cut%i",d, num),";m_{#mu#mu} (GeV)",  50, 0,20,  weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_0to20_b40_%s_cut%i",d, num),";m_{#mu#mu} (GeV)",  40, 0,20,  weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_0to50_%s_cut%i",    d, num),";m_{#mu#mu} (GeV)", 100, 0,50,  weight, dir);
    hists->fill1DHist(Mll, Form("01_diLep_mass_full_%s_cut%i", d, num),";m_{\\ell\\ell} (GeV)", 100, 0, 250, weight, dir);
    hists->fill1DHist(Mll, Form("01_diLep_mass_Z_%s_cut%i",    d, num),";m_{\\ell\\ell} (GeV)", 100, 60,120, weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_low1_%s_cut%i", d, num),";m_{ll} (GeV)",  50,   0,1.5,   weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_low2_%s_cut%i", d, num),";m_{ll} (GeV)", 100, 1.5,  8,   weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_low3_%s_cut%i", d, num),";m_{ll} (GeV)",  30,   8, 20,   weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_bump1_%s_cut%i",d, num),";m_{ll} (GeV)",  40,  15, 35,   weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_bump2_%s_cut%i",d, num),";m_{ll} (GeV)",  30,  35, 50,   weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_jpsi_%s_cut%i", d, num),";m_{ll} (GeV)",  50, 2.7,3.5,   weight, dir);
    //hists->fill1DHist(Mll, Form("01_diLep_mass_jpsi2_%s_cut%i",d, num),";m_{ll} (GeV)",  50, 2.5,3.7,   weight, dir);

    //hists->fill1DHist(TMath::Log10(Mll), Form("01_diLep_mass_log1_%s_cut%i", d, num),";log10(m_{ll})", 100,   -1,7, weight, dir);
    //hists->fill1DHist(TMath::Log10(Mll), Form("01_diLep_mass_log2_%s_cut%i", d, num),";log10(m_{ll})", 100,   -1,3, weight, dir);
    //hists->fill1DHist(TMath::Log10(Mll), Form("01_diLep_mass_log3_%s_cut%i", d, num),";log10(m_{ll})", 100,-0.7,1.3,weight, dir);

    hists->fill1DHist(diLep.Pt(),  Form("011_diLep_pt_%s_cut%i",   d, num),";di-Lepton p_{T}", 50, 0,120, weight, dir);
    hists->fill1DHist(diLep.Eta(), Form("011_diLep_eta_%s_cut%i",  d, num),";di-Lepton eta",   50, -3.5,3.5, weight, dir);
    hists->fill1DHist(diLep.Phi(), Form("011_diLep_phi_%s_cut%i",  d, num),";di-Lepton phi",   50, -TMath::Pi(),TMath::Pi(), weight, dir);
    //hists->fill1DHist(diLepCM.E(), Form("011_diLep_Ecom_%s_cut%i", d, num),";E(ll) in CoM",    50, 0,200,  weight, dir);

    hists->fill1DHist(_lPt1.Pt(),  Form("02_lPt1_pt_%s_cut%i", d, num), ";Leading lepton p_{T} (GeV)", 50, 0,100,    weight, dir);
    hists->fill1DHist(_lPt1.Eta(), Form("02_lPt1_eta_%s_cut%i",d, num), ";Leading lepton eta", 50, -3.5,3.5, weight, dir);
    hists->fill1DHist(_lPt1.Phi(), Form("02_lPt1_phi_%s_cut%i",d, num), ";Leading lepton phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);
    hists->fill1DHist(_lPt2.Pt(),  Form("02_lPt2_pt_%s_cut%i", d, num), ";Trailing lepton p_{T} (GeV)",50, 0,100,    weight, dir);
    hists->fill1DHist(_lPt2.Eta(), Form("02_lPt2_eta_%s_cut%i",d, num), ";Trailing lepton eta", 50, -3.5,3.5, weight, dir);
    hists->fill1DHist(_lPt2.Phi(), Form("02_lPt2_phi_%s_cut%i",d, num), ";Trailing lepton phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);

    hists->fill2DHist(_lPt1.Pt(), _lPt2.Pt(), Form("h2D_Pt2_vs_Pt1_%s_cut%i", d, num),
		      ";Leading lepton p_{T};Trailing lepton p_{T}",  50, 0,100, 50,0,100, weight, dir);



    //hists->fill1DHist(_lPt1.DeltaR(_lPt2), Form("04_ll_deltaR_low_%s_cut%i",    d, num),";#Delta R(l_{1}, l_{2})",25,0,0.05, weight,dir);
    //hists->fill1DHist(_lPt1.DeltaR(_lPt2), Form("04_ll_deltaR_%s_cut%i",        d, num),";#Delta R(l_{1}, l_{2})",50,  0,1,  weight,dir);
    hists->fill1DHist(_lPt1.DeltaR(_lPt2), Form("04_ll_deltaR_full_%s_cut%i",   d, num),";#Delta R(l_{1}, l_{2})",50,  0,4,  weight,dir);

  }

  if (_isLepSet && _isGammaSet){
    tri   = diLep + _gamma;

    Mllg  = tri.M();

    //hists->fill1DHist(Mllg, Form("00_tri_mass_%s_cut%i",         d, num),";m_{ll#gamma}",  100,40,240,  weight, dir);
    hists->fill1DHist(Mllg, Form("00_tri_mass_3GeVBin_%s_cut%i", d, num),";m_{ll#gamma}",  20,110,170,  weight, dir);

    //hists->fill1DHist(Mllg, Form("00_tri_mass125_%s_cut%i",      d, num),";m_{ll#gamma}",  30,110,140,  weight, dir);
    //hists->fill1DHist(Mllg, Form("00_tri_mass80_%s_cut%i",       d, num),";m_{ll#gamma}",  100,80,200,  weight, dir);
    hists->fill1DHist(Mllg, Form("00_tri_massZ_%s_cut%i",        d, num),";m_{ll#gamma}",  100,60,120,  weight, dir);
    hists->fill1DHist(Mllg, Form("00_tri_mass_longTail_%s_cut%i",d, num),";m_{ll#gamma}",  100, 60,360,  weight, dir);

    hists->fill1DHist(_gamma.Pt()/Mllg, Form("03_gamma_ptOverMllg_%s_cut%i",   d, num),
		      ";Photon p_{T}/m_{ll#gamma}",  100, 0,1, weight, dir);


    hists->fill1DHist(diLep.Pt()/Mllg,  Form("011_diLep_ptOverMllg_%s_cut%i",   d, num),
		      ";di-Lepton p_{T}/m_{ll#gamma}",100, 0,1, weight, dir);

    hists->fill1DHist(_lPt1.DeltaR(_gamma),Form("04_lPt1_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{1})",100,0,5, weight,dir);
    hists->fill1DHist(_lPt2.DeltaR(_gamma),Form("04_lPt2_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,l_{2})",100,0,5, weight,dir);
    hists->fill1DHist(diLep.DeltaR(_gamma),Form("04_ll_gamma_deltaR_%s_cut%i",d, num),";#Delta R(#gamma,ll)",100,0,5, weight,dir);



    /*
    TLorentzVector lPt1Gamma = _lPt1+_gamma;
    TLorentzVector lPt2Gamma = _lPt2+_gamma;
    hists->fill1DHist(lPt1Gamma.M(), Form("lPt1Gamma_%s_cut%i",d, num),";M(l1#gamma)",  100,60,160,  weight, dir);
    hists->fill1DHist(lPt2Gamma.M(), Form("lPt2Gamma_%s_cut%i",d, num),";M(l2#gamma)",  100, 0,150,  weight, dir);
    
    hists->fill2DHist(diLep.M2(), lPt2Gamma.M2(), Form("h2D_dalitzPlot_ll_vs_gl_%s_cut%i",  d, num),
		      ";m^{2}_{ll} (GeV^{2});m^{2}_{l_{2}#gamma} (GeV^{2})", 100,0,20000, 100, 0,10000,  weight, dir);
    
    hists->fill2DHist(lPt1Gamma.M2(), lPt2Gamma.M2(), Form("h2D_dalitzPlot_l1g_vs_l2g_%s_cut%i",  d, num),
		      ";m^{2}_{l_{1}#gamma} (GeV^{2});m^{2}_{l_{2}#gamma} (GeV^{2})", 100,0,20000, 100, 0,10000,  weight, dir);
    
    double rotation = atan(-1.);
    double sin_rotation = sin(rotation);
    double cos_rotation = cos(rotation);
    double x_prime = (lPt1Gamma.M2()*cos_rotation - lPt2Gamma.M2()*sin_rotation) + (1.-cos_rotation)*pow(125.,2);
    double y_prime =  lPt1Gamma.M2()*sin_rotation + lPt2Gamma.M2()*cos_rotation;
    
    hists->fill2DHist(x_prime, y_prime, Form("h2D_dalitzPlot_l1g_vs_l2g_rotation_%s_cut%i",  d, num),
		      ";M(l1#gamma)^{2};M(l2#gamma)^{2}", 100,0,22000, 100, -14000,5000,  weight, dir);
    
    */

    /*
    const UInt_t nBins1 = 60;
    hists->fill2DHist(Mllg, Mll, Form("h2D_tri_vs_diLep_mass0_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 0,20,  weight, dir);
    hists->fill2DHist(Mllg, Mll, Form("h2D_tri_vs_diLep_mass1_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 0,1.5, weight, dir);
    hists->fill2DHist(Mllg, Mll, Form("h2D_tri_vs_diLep_mass2_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 1.5,8, weight, dir);
    hists->fill2DHist(Mllg, Mll, Form("h2D_tri_vs_diLep_mass3_%s_cut%i",  d, num),
		      ";m_{ll#gamma} (GeV);m_{ll} (GeV)", nBins1,0,200, nBins1, 8,20,  weight, dir);
    */

    /*
    if (_isGamma2Set){
      four  = diLep + _gamma + _gamma2;
      hists->fill1DHist(four.M(), Form("four_massZ_%s_cut%i",d, num),";m_{ll#gamma#gamma}",  100,70,120,  weight, dir);
    }
    */


    TLorentzVector l1,l2;
    l1 = _lPt1; l2 = _lPt2;

    /*
    // l1 has to be positive and l2 is negative! charge!

    if (_lPt1.Charge()==1 && _lPt2.Charge()==-1)
      {l1 = _lPt1; l2 = _lPt2;}
    else if (_lPt1.Charge()==-1 && _lPt2.Charge()==1)
      {l1 = _lPt2; l2 = _lPt1;}
    else{
      return;
      //co1=co2=phi=co3=-5;
      //cout<<"Something is wrong: leptons have the same charge!"<<endl;
    }
    //cout<<eventNumber<<"RECO  Angles: c1= "<<co1<<"  c2="<<co2<<"   phi="<<phi<<"   coco="<<co3<<endl;

    */

    
    //ZGanlgles:
    double co1,co2,phi,co3;
    angles->GetZGAngles(l1, l2, _gamma,co1,co2,phi,co3);
    
    hists->fill1DHist(co1, Form("co1_cut%i", num), ";cos(#vec{l-},#vec{Z}) in CM", 100,-1,1, weight,"Angles");
    hists->fill1DHist(co2, Form("co2_cut%i", num), ";cos(#vec{l+},#vec{Z}) in CM", 100,-1,1, weight,"Angles");
    hists->fill1DHist(co3, Form("co3_cut%i", num), ";cos(#Theta)",                 100,-1,1, weight,"Angles");
    hists->fill1DHist(phi, Form("phi_cut%i", num), ";#phi(l+)",50, -TMath::Pi(), TMath::Pi(),weight,"Angles");
    
    
    
  }

  /*
  Float_t MdalG = 0;

  if(_isGammaSet)
    {
      if (_isDaleSet){
	MdalG = (_dale+_gamma).M();

	hists->fill1DHist(MdalG-Mllg, Form("01_mllg-diff_%s_cut%i", d, num),
			  ";m_{e'#gamma} - m_{ll#gamma} (GeV)", 60,-30,30, weight, dir+"-Dale");
	hists->fill1DHist(MdalG,   Form("01_mDalG_seeZ_%s_cut%i",d, num),";m_{e'#gamma} (GeV)",  50, 80,160, weight, dir+"-Dale");
	hists->fill1DHist(MdalG,   Form("01_mDalG_full_%s_cut%i",d, num),";m_{e'#gamma} (GeV)", 100, 50,200, weight, dir+"-Dale");
	hists->fill1DHist(MdalG,   Form("01_mDalG_fit_%s_cut%i", d, num),";m_{e'#gamma} (GeV)", 50, 110,170, weight, dir+"-Dale");
	hists->fill1DHist(MdalG,   Form("01_mDalG_m125_%s_cut%i",d, num),";m_{e'#gamma} (GeV)", 50, 110,140, weight, dir+"-Dale");
	hists->fill1DHist(_gamma.Pt()/MdalG, Form("02_gamma_ptOverMdalG_%s_cut%i",  d, num),
			  ";Photon p_{T}/m_{e'#gamma}",  100, 0,1, weight, dir+"-Dale");
	hists->fill1DHist(_dale.Pt()/MdalG,  Form("02_dale_ptOverMdalG_%s_cut%i",   d, num),
			  ";Dalitz Electron p_{T}/m_{e'#gamma}",  100, 0,1, weight, dir+"-Dale");

	hists->fill1DHist(_dale.DeltaR(_gamma), Form("03_deltaR_dale_gamma_%s_cut%i",d, num),
			  ";#Delta R(e', #gamma)",  50, 0,5, weight, dir+"-Dale");
	hists->fill1DHist(_dale.Pt(), Form("03_dale_pt_%s_cut%i",  d, num),
			  ";Dalitz Electron p_{T} (GeV)",  50, 0,120, weight, dir+"-Dale");
	hists->fill1DHist(_dale.Eta(),Form("03_dale_eta_%s_cut%i", d, num),
			  ";Dalitz Electron eta", 50, -3.5,3.5,       weight, dir+"-Dale");
	hists->fill1DHist(_dale.Phi(),Form("03_dale_phi_%s_cut%i", d, num),
			  ";Dalitz Electron phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir+"-Dale");

	hists->fill1DHist(_dale.Pt()/_gamma.Pt(),  Form("03_ptDaleOverGamma_%s_cut%i",   d, num),
			  ";p_{T}^{e'}/p_{T}^{#gamma}",  100, 0,3, weight, dir+"-Dale");

      }

      TLorentzVector gammaCM = _gamma;
      //TLorentzVector diLepCM = diLep;
      TVector3 b1  = -tri.BoostVector();
      gammaCM.Boost(b1);
      //diLepCM.Boost(b1);

      hists->fill1DHist(gammaCM.E(), Form("03_gamma_Ecom_%s_cut%i",d, num),";Photon Energy in CoM",50, 0,120, weight, dir);
      hists->fill1DHist(_gamma.E(),  Form("03_gamma_E_%s_cut%i",   d, num),";Photon Energy",       50, 0,200, weight, dir);
      hists->fill1DHist(_gamma.E(),  Form("03_gamma_E_hi_%s_cut%i",d, num),";Photon Energy",       50, 0,400, weight, dir);
      hists->fill1DHist(_gamma.Pt(), Form("03_gamma_pt_%s_cut%i",  d, num),";Photon p_{T} (GeV)",  50, 0,120, weight, dir);
      hists->fill1DHist(_gamma.Eta(),Form("03_gamma_eta_%s_cut%i", d, num),";Photon eta", 50, -3.5,3.5,       weight, dir);
      hists->fill1DHist(_gamma.Phi(),Form("03_gamma_phi_%s_cut%i", d, num),";Photon phi", 50, -TMath::Pi(),TMath::Pi(), weight, dir);

      hists->fill2DHist(diLep.Pt(),_gamma.Pt(),Form("h2D_diLep_vs_gamma_%s_cut%i",d, num),
			";q_{T}^{ll} (GeV); p_{T}^{#gamma} (GeV)", 50, 0,100, 50,0,100, weight, dir);
      hists->fill2DHist(diLep.Pt()/Mllg,_gamma.Pt()/Mllg,Form("h2D_diLep_vs_gamma_ptOverMllg_%s_cut%i",d, num),
			";q_{T}^{ll}/M_{ll#gamma}; p_{T}^{#gamma}/M_{ll#gamma}",  100, 0,1, 100,0,1, weight, dir);
      hists->fill2DHist(gammaCM.E(), _gamma.Pt(),Form("h2D_gamma_Ecom_vs_Pt_%s_cut%i",  d, num),
			";E_{#gamma} in CoM; p_{T}(#gamma)", 50, 0,100, 50,0,100, weight, dir);
      hists->fill2DHist(gammaCM.E(), Mllg,   Form("h2D_gamma_Ecom_vs_triM_%s_cut%i",d, num),
			";E_{#gamma} in CoM; m_{ll#gamma}",   50, 0,100, 50,0,200, weight, dir);
      hists->fill2DHist(diLep.Eta()-_gamma.Eta(), _gamma.Pt(), Form("h2D_deltaEta_vs_gammaPt_deltaEta_%s_cut%i",d, num),
			";#Delta#eta(ll, #gamma);p_{T} of #gamma", 100, -5,5, 100,0,130, weight, dir);

    }


  if (_isJet1Set){
    hists->fill1DHist(_jet1.Pt(),  Form("050_jet1_pt_%s_cut%i",  d, num), ";p_{T}^{jet1} (GeV)", 50,0,200,  weight,dir);
    hists->fill1DHist(_jet1.Eta(), Form("050_jet1_eta_%s_cut%i", d, num), ";#eta^{jet1}",        50, -5,5,  weight,dir);

    hists->fill1DHist(_jet1.DeltaR(diLep), Form("051_jet1_deltaR_diLep_%s_cut%i", d, num),
		      ";#Delta R(j_{1}, di-Lepton)", 50,0,6,  weight,dir);
    if(_isGammaSet)
      hists->fill1DHist(_jet1.DeltaR(_gamma), Form("051_jet1_deltaR_gamma_%s_cut%i", d, num),
			";#Delta R(j_{1}, #gamma)", 50,0,6,  weight,dir);
  }

  if (_isJet2Set){
    hists->fill1DHist(_jet2.Pt(),  Form("050_jet2_pt_%s_cut%i",  d, num), ";p_{T}^{jet2} (GeV)", 50,0,200,  weight,dir);
    hists->fill1DHist(_jet2.Eta(), Form("050_jet2_eta_%s_cut%i", d, num), ";#eta^{jet2}",        50, -5,5,  weight,dir);

    hists->fill1DHist(_jet2.DeltaR(diLep), Form("051_jet2_deltaR_diLep_%s_cut%i", d, num),
		      ";#Delta R(j_{2}, di-Lepton)", 50,0,6,  weight,dir);
    if(_isGammaSet)
      hists->fill1DHist(_jet2.DeltaR(_gamma), Form("051_jet2_deltaR_gamma_%s_cut%i", d, num),
			";#Delta R(j_{2}, #gamma)", 50,0,6,  weight,dir);
  }

  if (_isJet1Set && _isJet2Set){
    hists->fill1DHist(fabs(_jet1.Eta() -_jet2.Eta()), Form("052_deltaEta_j1_j2_%s_cut%i", d, num),
		      ";|#Delta #eta(j_{1}, j_{2})|", 50,0,7,  weight,dir);
    hists->fill1DHist(_jet1.DeltaR(_jet2), Form("052_deltaR_j1_j2_%s_cut%i", d, num),
		      ";#Delta R(j_{1}, j_{2})",   50,0,7,  weight,dir);
    hists->fill1DHist((_jet1+_jet2).M(), Form("052_mass_j1j2_%s_cut%i", d, num),
		      ";m(j_{1}, j_{2}) (GeV)", 100,0,800,  weight,dir);
    if (_isGammaSet){
      float zep = Zeppenfeld((_lPt1+_lPt2+_gamma),_jet1,_jet2);
      hists->fill1DHist(zep, Form("052_zep_%s_cut%i", d, num),
			";Zeppenfeld: #eta_{ll#gamma} - #frac{1}{2}(#eta_{j1} + #eta_{j2})", 50,-6,6,  weight,dir);
      hists->fill1DHist(fabs((_jet1+_jet2).DeltaPhi(_lPt1+_lPt2+_gamma)), Form("052_dPhi_diJet_higgs_%s_cut%i", d, num),
			";#Delta#phi(di-Jet, ll#gamma)", 50,0,3,  weight,dir);

    }
  }
  if (_isMetSet){
    hists->fill1DHist(_met.Mod(), Form("06_met_Et_%s_cut%i",  d, num),";E_{T}^{miss} (GeV)", 50,0,100,  weight,dir);
    hists->fill1DHist(_met.Phi()-TMath::Pi(), Form("06_met_phi_%s_cut%i", d, num),
		      ";#phi(E_{T}^{miss})", 50,-TMath::Pi(),TMath::Pi(),  weight,dir);

    if (_isJet1Set)
      hists->fill1DHist(fabs(_met.DeltaPhi(_jet1.Vect().XYvector())), Form("06_met_dPhiJet1_%s_cut%i", d, num),
			";#Delta#phi(E_{T}^{miss}, j_{1})", 50,0,TMath::Pi(),  weight,dir);
    if (_isJet2Set)
      hists->fill1DHist(fabs(_met.DeltaPhi(_jet2.Vect().XYvector())), Form("06_met_dPhiJet2_%s_cut%i", d, num),
			";#Delta#phi(E_{T}^{miss}, j_{2})", 50,0,TMath::Pi(),  weight,dir);

  }
  */

  hists->fill1DHist(weight, Form("weight_%s_cut%i", d, num), ";weight", 200, 0,3, 1, dir);
  hists->fill1DHist(_rhoFactor, Form("rhoFactor_%s_cut%i",     d, num), ";#rho-factor",     100, 0,35, weight, dir);
  hists->fill1DHist(_nVtx,      Form("vtx_nPV_weight_%s_cut%i",d, num), ";nPV, re-weighted", 40, 0,40, weight, dir);
  hists->fill1DHist(_nVtx,      Form("vtx_nPV_raw_%s_cut%i",   d, num), ";nPV, Raw",         40, 0,40,      1, dir);
  //hists->fill1DHist(nVtxTotal, Form("vtx_nPV_tot_%s_cut%i",   d, num), "vtx_nPV_tot",    40, 0,40, weight, dir);

  hists->fill1DHist(_pv.ndof(), Form("vtx_ndof1_%s_cut%i",  d, num), ";First vertex nDof", 50, 0,200, weight, dir);
  //hists->fill1DHist(nDofVtx2, Form("vtx_ndof2_%s_cut%i",  d, num), ";Second vertex ndof", 50, 0,200, weight, dir);


}
    //Phi* variable
  /*
    float phi_acop   = TMath::Pi() - fabs(l1.DeltaPhi(l2));
    float theta_star = TMath::ACos(TMath::TanH( (l2.Eta() - l1.Eta())/2));
    float phi_star   = sin(theta_star)*TMath::Tan(phi_acop/2);

    hists->fill1DHist(phi_acop, Form("star_phi_acop_cut%i", num), ";#phi_{acop}",  50, 0, TMath::Pi(), weight,"Angles");
    hists->fill1DHist(phi_star, Form("star_phi_star_cut%i", num), ";#phi*",        50, 0,         150, weight,"Angles");
    hists->fill1DHist(cos(theta_star), Form("star_cos_theta_cut%i", num), ";cos(#theta*)",  50, -1, 1, weight,"Angles");


    //Forward-Backward assymetry angles
    int sign = (diLep.Z() > 0) ? 1 : -1;
    float m = Mll;
    float cosCS = sign*2*(l2.Pz()*l1.E() - l1.Pz()*l2.E())/(m*sqrt(m + pow(diLep.Pt(),2)));

    hists->fill1DHist(cosCS, Form("FB_cosSC-ll_cut%i",  num), ";cosCS(l+,l-)", 100, -1, 1, weight,"Angles");
    hists->fill2DHist(cosCS, m, Form("FB_cosSC-ll_vs_Mll_cut%i",  num), ";cosCS(l+,l-);m_{ll}",  100, -1,1, 100,0,20, weight,"Angles");

    if (_isGammaSet){
      sign = (tri.Z() > 0) ? 1 : -1;
      m = Mllg;
      cosCS = sign*2*(diLep.Pz()*_gamma.E() - _gamma.Pz()*diLep.E())/(m*sqrt(m + pow(tri.Pt(),2)));

      hists->fill1DHist(cosCS, Form("FB_cosSC-diG_cut%i", num), ";cosCS(diLep,#gamma)",  100, -1, 1, weight,"Angles");
      hists->fill2DHist(cosCS, m, Form("FB_cosSC-diG_vs_Mllg_cut%i", num),
			";cosCS(diLep,#gamma);m_{ll#gamma}",  100, -1,1, 100,50,140, weight,"Angles");

      float cosJ = _gamma.Dot(l1-l2)/_gamma.Dot(l1+l2);
      hists->fill1DHist(cosJ, Form("FB_cosJ_cut%i", num), ";cosJ",  100, -1, 1, weight,"Angles");
      hists->fill2DHist(cosJ, Mllg, Form("FB_cosJ_vs_Mllg_cut%i", num), ";cosJ;m_{ll#gamma}",  100, -1,1, 100,50,140, weight,"Angles");
    }
}
  */
