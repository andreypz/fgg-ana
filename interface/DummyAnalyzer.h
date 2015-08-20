#include <map>
#include <string>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
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

#define nC 10

class DummyAnalyzer : public edm::BasicAnalyzer {
 public:
  /// default constructor
  DummyAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs);
  /// default destructor
  virtual ~DummyAnalyzer(){
    delete hists;
  };
  void beginJob();
  void endJob();
  void analyze(const edm::EventBase& event);

  virtual void CountEvents(Int_t , std::string , Double_t , std::ofstream& s);
  virtual void FillHistoCounts(Int_t n, Double_t w);

 protected:

  UInt_t nEvents[nC];
  UInt_t totEvents;
  UInt_t count_neg, count_pos;


  Double_t nWeights[nC];
  Double_t totWeights;
  Double_t W0;

  std::ofstream fout, fcuts;
  std::string lep_;
  double lumiWeight_;

  edm::InputTag muons_, electrons_, photons_;
  edm::InputTag myGen_;

  //edm::EDGetTokenT<std::vector<reco::Muon> > muonsToken_;

  bool isRealData;
  std::map<std::string, TH1*> hists_;
  //TFile *outFile;
  //TFile *histoFile;
  HistManager *hists;
};
