#ifndef _FggHistMakerZgamma_H
#define _FggHistMakerZgamma_H

#include <cstring>
#include "FggHistMakerBase.h"
#include "Angles.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class FggHistMakerZgamma : public FggHistMakerBase {

 public:
  FggHistMakerZgamma(HistManager * h);
  virtual ~FggHistMakerZgamma();

  virtual void MakeMainHistos(Int_t n, Double_t w, string s="Main");
  virtual void SetVtx(reco::Vertex v);
  virtual void SetLeptons(TLorentzVector l1, TLorentzVector l2);
  virtual void SetGamma(TLorentzVector );
  virtual void SetGamma2(TLorentzVector );
  virtual void SetJet1(TLorentzVector );
  virtual void SetJet2(TLorentzVector );
  virtual void SetMet(TVector2 );
  virtual void SetDalectron(TLorentzVector );
  virtual void Reset(float r, UInt_t n);

 private:

  TLorentzVector _lPt1, _lPt2, _gamma, _gamma2, _dale, _jet1, _jet2;
  TVector2 _met;
  reco::Vertex _pv;
  Bool_t _isVtxSet;
  Bool_t _isLepSet;
  Bool_t _isGammaSet;
  Bool_t _isGamma2Set;
  Bool_t _isDaleSet;
  Bool_t _isJet1Set;
  Bool_t _isJet2Set;
  Bool_t _isMetSet;
  Float_t _rhoFactor;
  UInt_t _nVtx;
  
  Angles *angles;

};

#endif
