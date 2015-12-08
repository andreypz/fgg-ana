#ifndef _DUMMY_ANA
#define _DUMMY_ANA

#include <map>
#include <string>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "HistManager.h"

#define nC 13

class DummyAnalyzer : public edm::BasicAnalyzer {
 public:
  /// default constructor
  DummyAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs);
  /// default destructor
  virtual ~DummyAnalyzer(){
    delete hists;
  };
  void beginJob();
  void endJob(UInt_t n);
  void analyze(const edm::EventBase& event);

  virtual std::string open_temp(std::string p, std::string n, std::ofstream& s);
  virtual void CountEvents(Int_t , std::string , Double_t , std::ofstream& s);
  virtual void FillHistoCounts(Int_t n, Double_t w);

  static bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());}

 protected:

  UInt_t nEvents[nC];
  UInt_t totEvents;
  UInt_t count_neg, count_pos;
  

  Double_t nWeights[nC];
  Double_t totWeights;
  Double_t W0;

  std::ofstream fout, fcuts;
  std::string yields_txt_path;
  std::string lep_;
  double lumiWeight_;

  edm::InputTag muons_, electrons_, photons_, jets_;
  edm::InputTag myGen_;

  //edm::EDGetTokenT<std::vector<reco::Muon> > muonsToken_;

  Bool_t isRealData;
  ULong_t eventNumber;
  UInt_t runNumber;
  std::map<std::string, TH1*> hists_;
  //TFile *outFile;
  //TFile *histoFile;
  HistManager *hists;
};

#endif /* _DUMMY_ANA */
