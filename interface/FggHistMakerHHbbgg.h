#ifndef _FggHistMakerHHbbgg_H
#define _FggHistMakerHHbbgg_H

#include <cstring>
#include "FggHistMakerBase.h"
#include "Angles.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class FggHistMakerHHbbgg : public FggHistMakerBase {

 public:
  FggHistMakerHHbbgg(HistManager * h);
  virtual ~FggHistMakerHHbbgg();

  virtual void MakeMainHistos(Int_t n, Double_t w, string s="Main");
  virtual void SetVtx(reco::Vertex v);
  virtual void SetGamma1(TLorentzVector );
  virtual void SetGamma2(TLorentzVector );
  virtual void SetBJet1(TLorentzVector );
  virtual void SetBJet2(TLorentzVector );
  virtual void SetMet(TVector2 );
  virtual void Reset(float r, UInt_t n);

 private:

  TLorentzVector _gamma1, _gamma2, _bjet1, _bjet2;
  TVector2 _met;
  reco::Vertex _pv;
  Bool_t _isVtxSet;
  Bool_t _isGamma1Set;
  Bool_t _isGamma2Set;
  Bool_t _isBJet1Set;
  Bool_t _isBJet2Set;
  Bool_t _isMetSet;
  Float_t _rhoFactor;
  UInt_t _nVtx;
  
  ZGAngles *angles;

};

#endif
