#include "../interface/FggHistMakerHHbbgg.h"

FggHistMakerHHbbgg::FggHistMakerHHbbgg(HistManager *h):
  FggHistMakerBase::FggHistMakerBase(h),
  _isVtxSet(false),
  _isGamma1Set(false),
  _isGamma2Set(false),
  _isBJet1Set(false),
  _isBJet2Set(false),
  _isMetSet(false),
  _rhoFactor(-1),
  _nVtx(0)
{
  cout<<"\t HHHHHHH \t HHbbggHistmaker constractor"<<endl;
  //angles = new ZGAngles();
}

FggHistMakerHHbbgg::~FggHistMakerHHbbgg(){}

void FggHistMakerHHbbgg::Reset(float rho, UInt_t n){
  _isGamma1Set=0; _isGamma2Set=0;
  _isBJet1Set=0; _isBJet2Set=0; _isMetSet=0;
  _rhoFactor = rho; _nVtx=n;}

void FggHistMakerHHbbgg::SetVtx(reco::Vertex v){
  _pv=v; _isVtxSet = 1;}
void FggHistMakerHHbbgg::SetGamma1(TLorentzVector g){
  _gamma1 = g; _isGamma1Set=1;}
void FggHistMakerHHbbgg::SetGamma2(TLorentzVector g){
  _gamma2 = g; _isGamma2Set=1;}
void FggHistMakerHHbbgg::SetBJet1(TLorentzVector j){
  _bjet1 = j; _isBJet1Set=1;}
void FggHistMakerHHbbgg::SetBJet2(TLorentzVector j){
  _bjet2 = j; _isBJet2Set=1;}
void FggHistMakerHHbbgg::SetMet(TVector2 m){
  _met = m; _isMetSet=1;}


void FggHistMakerHHbbgg::MakeMainHistos(Int_t num, Double_t weight, string dir)
{
  //cout<<num<<" Making plots "<<dir<<endl;
  char * d = new char [dir.length()+1];
  std::strcpy (d, dir.c_str());


  if (_isGamma1Set && _isGamma2Set){
    TLorentzVector Hgg = _gamma1 + _gamma2;
    Float_t Mgg  = Hgg.M();
    hists->fill1DHist(Mgg, Form("01_Mgg_%s_cut%i", d, num),";m_{#gamma#gamma}", 60,110,170,  weight, dir);

    hists->fill1DHist(Hgg.Pt(), Form("01_pT_gg_%s_cut%i", d, num),";p_{T}^{#gamma#gamma}", 100,0,600,  weight, dir);

  }

  if (_isBJet1Set && _isBJet2Set){

    TLorentzVector Hbb = _bjet1 + _bjet2;
    Float_t Mbjbj  = Hbb.M();

    hists->fill1DHist(Mbjbj, Form("02_Mbjbj6_%s_cut%i", d, num),";m(jj), GeV", 50, 50,350, weight, dir);
    hists->fill1DHist(Hbb.Pt(), Form("02_pT_bjbj_%s_cut%i", d, num),";p_{T}^{jj}", 100,0,600,  weight, dir);

    hists->fill1DHist(fabs(_bjet1.Eta() -_bjet2.Eta()), Form("03_deltaEta_j1_j2_%s_cut%i", d, num),
		      ";|#Delta #eta(j_{1}, j_{2})|", 50,0,4,  weight,dir);
    hists->fill1DHist(_bjet1.DeltaR(_bjet2), Form("03_deltaR_j1_j2_%s_cut%i", d, num),
		      ";#Delta R(j_{1}, j_{2})",   50,0,5,  weight,dir);

    hists->fill2DHist(_bjet1.Eta()-_bjet2.Eta(), _bjet1.DeltaPhi(_bjet2), Form("03_dEta_dPhi_%s_cut%i", d, num),
		      ";#DeltaEta(j1,j2);#DeltaPhi(j1,j2)", 50, -4,4, 50, -3.15, 3.15, weight, dir);

  }

  if (_isGamma1Set && _isGamma2Set && _isBJet1Set && _isBJet2Set){

    Float_t Mbbgg = (_gamma1 + _gamma2 + _bjet1 + _bjet2).M();
    
    hists->fill1DHist(Mbbgg, Form("00_Mbbgg_r1_%s_cut%i", d, num),";m(#gamma#gamma jj), GeV", 100, 100,900, weight, dir);
    hists->fill1DHist(Mbbgg, Form("00_Mbbgg_r2_%s_cut%i", d, num),";m(#gamma#gamma jj), GeV", 100, 0,1500, weight, dir);
  }



  //if (_isGammaSet){
  //  float zep = Zeppenfeld((_lPt1+_lPt2+_gamma),_bjet1,_bjet2);
  //  hists->fill1DHist(zep, Form("052_zep_%s_cut%i", d, num),
  //			";Zeppenfeld: #eta_{ll#gamma} - #frac{1}{2}(#eta_{j1} + #eta_{j2})", 50,-6,6,  weight,dir);
  //  hists->fill1DHist(fabs((_bjet1+_bjet2).DeltaPhi(_lPt1+_lPt2+_gamma)), Form("052_dPhi_diJet_higgs_%s_cut%i", d, num),
  //";#Delta#phi(di-Jet, ll#gamma)", 50,0,3,  weight,dir);


  /*


  }
  if (_isMetSet){
    hists->fill1DHist(_met.Mod(), Form("06_met_Et_%s_cut%i",  d, num),";E_{T}^{miss} (GeV)", 50,0,100,  weight,dir);
    hists->fill1DHist(_met.Phi()-TMath::Pi(), Form("06_met_phi_%s_cut%i", d, num),
		      ";#phi(E_{T}^{miss})", 50,-TMath::Pi(),TMath::Pi(),  weight,dir);

    if (_isJet1Set)
      hists->fill1DHist(fabs(_met.DeltaPhi(_bjet1.Vect().XYvector())), Form("06_met_dPhiJet1_%s_cut%i", d, num),
			";#Delta#phi(E_{T}^{miss}, j_{1})", 50,0,TMath::Pi(),  weight,dir);
    if (_isJet2Set)
      hists->fill1DHist(fabs(_met.DeltaPhi(_bjet2.Vect().XYvector())), Form("06_met_dPhiJet2_%s_cut%i", d, num),
			";#Delta#phi(E_{T}^{miss}, j_{2})", 50,0,TMath::Pi(),  weight,dir);

  }
  */

  //hists->fill1DHist(weight, Form("weight_%s_cut%i", d, num), ";weight", 200, 0,3, 1, dir);
  //hists->fill1DHist(_rhoFactor, Form("rhoFactor_%s_cut%i",     d, num), ";#rho-factor",     100, 0,35, weight, dir);
  //hists->fill1DHist(_nVtx,      Form("vtx_nPV_weight_%s_cut%i",d, num), ";nPV, re-weighted", 40, 0,40, weight, dir);
  //hists->fill1DHist(_nVtx,      Form("vtx_nPV_raw_%s_cut%i",   d, num), ";nPV, Raw",         40, 0,40,      1, dir);
  ////hists->fill1DHist(nVtxTotal, Form("vtx_nPV_tot_%s_cut%i",   d, num), "vtx_nPV_tot",    40, 0,40, weight, dir);

  //hists->fill1DHist(_pv.ndof(), Form("vtx_ndof1_%s_cut%i",  d, num), ";First vertex nDof", 50, 0,200, weight, dir);
}
