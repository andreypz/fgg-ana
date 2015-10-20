#ifndef _GEN_ANA
#define _GEN_ANA

#include "DummyAnalyzer.h"

#include "FggHistMakerHHbbgg.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"

class GenAnalyzer : public DummyAnalyzer {

 public:
  /// default constructor
  GenAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs);
  /// default destructor
  virtual ~GenAnalyzer(){};
  void beginJob();
  void endJob();
  void analyze(const edm::EventBase& event);

 private:
  FggHistMakerHHbbgg *FHM;

};


#endif /* GEN_ANA */
