#include "DummyAnalyzer.h"

#include "FggHistMakerZgamma.h"

class ZGammaAnalyzer : public DummyAnalyzer {

 public:
  /// default constructor
  ZGammaAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs);
  /// default destructor
  virtual ~ZGammaAnalyzer(){};
  void beginJob();
  void endJob();
  void analyze(const edm::EventBase& event);

 private:

  FggHistMakerZgamma *FHM;
  Double_t a;
};
