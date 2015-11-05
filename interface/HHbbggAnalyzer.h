#ifndef _HHbbgg_ANA
#define _HHbbgg_ANA

#include "DummyAnalyzer.h"

#include "FggHistMakerHHbbgg.h"
#include "flashgg/bbggTools/interface/bbggTools.h"
#include "flashgg/DataFormats/interface/Jet.h"

class HHbbggAnalyzer : public DummyAnalyzer {

 public:

  HHbbggAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs);

  virtual ~HHbbggAnalyzer(){};
  void beginJob();
  void endJob();
  void analyze(const edm::EventBase& event);

  typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

 private:
  FggHistMakerHHbbgg *FHM;
  bbggTools tools_;

  // ORDER MATTERS:
  std::vector<edm::InputTag> inputTagJets_;
  std::vector<double> phoIDcutEB_;
  std::vector<double> phoIDcutEE_;
  
};


#endif /* HHbbgg_ANA */
