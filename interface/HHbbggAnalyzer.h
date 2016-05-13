#ifndef _HHbbgg_ANA
#define _HHbbgg_ANA

#include "DummyAnalyzer.h"

#include "FggHistMakerHHbbgg.h"
#include "flashgg/bbggTools/interface/bbggTools.h"
#include "flashgg/DataFormats/interface/Jet.h"
//#include "TTree.h"

class HHbbggAnalyzer : public DummyAnalyzer {

 public:

  HHbbggAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs);

  virtual ~HHbbggAnalyzer(){};
  void beginJob();
  void endJob();
  void analyze(const edm::EventBase& event);

  //typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;
  /*
  struct bTagSort{
  bTagSort( const HHbbggAnalyzer& info ) : m_info(info) { }
    // This is in order to sort the vector of jets by b-discriminatr,
    // isng the bTag from the class memeber.
    // Taken from here:
    // http://stackoverflow.com/questions/1902311/problem-sorting-using-member-function-as-comparator
    const HHbbggAnalyzer& m_info;
    bool operator()( const flashgg::Jet& p1, const flashgg::Jet& p2){
      return (p1.bDiscriminator(m_info.bTag) > p2.bDiscriminator(m_info.bTag));
    }
  };
  */

  // This does not work because the static functions cant use the class member variables
  //statis bool bTagSort(const flashgg::Jet& p1, const flashgg::Jet& p2){
  //return (p1.bDiscriminator(bTag) > p2.bDiscriminator(bTag));}


 private:
  FggHistMakerHHbbgg *FHM;
  bbggTools *tools;

  /*TTree *outTree;
  Char_t    o_category;
  UInt_t    o_run;
  ULong64_t o_evt;
  Double_t  o_weight;
  Double_t  o_bbMass;
  Double_t  o_ggMass;
  Double_t  o_bbggMass;
  */
  
  string bTagName;
  edm::InputTag rhoFixedGrid_;
  // ORDER MATTERS:
  //std::vector<edm::InputTag> inputTagJets_;
  std::vector<std::string> myTriggers_;
  std::vector<double> phoIDcutEB_;
  std::vector<double> phoIDcutEE_;
  std::vector<double> phoISOcutEB_;
  std::vector<double> phoISOcutEE_;
  //edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;

  UInt_t cutFlow_;
  Bool_t useDiPhotons_;
  edm::InputTag diPhotons_;
  UInt_t phoIDtype_;
};


#endif /* HHbbgg_ANA */
