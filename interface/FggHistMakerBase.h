#ifndef _FggHistMakerBase_H
#define _FggHistMakerBase_H

#include <cstring>
#include "HistManager.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Jet.h"

class FggHistMakerBase {
 public:
  FggHistMakerBase(HistManager * h);
  virtual ~FggHistMakerBase();

  void SetEventNumber(ULong_t e);
  template <class EGammaObj>
    void MakeEGammaCommonPlots(const EGammaObj& e, TString n);
  virtual void MakeElectronPlots(const flashgg::Electron& el, string s="Ele");
  virtual void MakePhotonPlots(const flashgg::Photon& ph, string s="Photon");
  
  virtual void MakeJetPlots(const flashgg::Jet& j, string s="Jets");
  virtual void MakeMuonPlots(const flashgg::Muon& mu, string s="Muons");
  virtual void MakeZeePlots(const flashgg::Photon& , const flashgg::Photon& );
  //virtual void MakePhotonEnergyCorrPlots(const flashgg::Photon& p, Float_t , Float_t );
  virtual void MakeNPlots(Int_t , UInt_t , UInt_t , UInt_t , Double_t w);
  float Zeppenfeld(const TLorentzVector& p, const TLorentzVector& pj1, const TLorentzVector& pj2);
  
 protected:
  HistManager * hists;
  ULong_t _eventNumber;
};

#endif
