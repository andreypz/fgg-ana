#ifndef _GEN_ANA
#define _GEN_ANA

#include "DummyAnalyzer.h"

#include "FggHistMakerBase.h"

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
  FggHistMakerBase *FHM;

};


#endif /* GEN_ANA */
