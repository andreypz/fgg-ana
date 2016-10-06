#ifndef _HHbbgg_ANA
#define _HHbbgg_ANA

#include "DummyAnalyzer.h"

#include "FggHistMakerHHbbgg.h"
#include "flashgg/bbggTools/interface/bbggTools.h"
#include "flashgg/bbggTools/interface/bbggJetRegression.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/Taggers/interface/GlobalVariablesDumper.h"

#include "TTree.h"
#include "Angles.h"


typedef math::XYZTLorentzVector LorentzVector;

class HHbbggAnalyzer : public DummyAnalyzer {

 public:

  HHbbggAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs);

  virtual ~HHbbggAnalyzer(){};
  void beginJob();
  void endJob();
  void analyze(const edm::EventBase& event);

  LorentzVector  LorentzToLorentz(const TLorentzVector& v);
  TLorentzVector LorentzToLorentz(const LorentzVector& v);


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
  bbggJetRegression *jetReg;
  
  flashgg::GlobalVariablesDumper* globVar_;

  //Tree objects
  TTree *flatTree, *genTree;
  Char_t    o_category;
  UInt_t    o_run;
  ULong64_t o_evt;
  Double_t  o_weight;
  Double_t  o_bbMass, o_ggMass, o_bbggMass;

  LorentzVector leadingPhoton, subleadingPhoton, diphotonCandidate;
  LorentzVector leadingJet, subleadingJet, dijetCandidate;
  LorentzVector leadingJet_KF, subleadingJet_KF, dijetCandidate_KF;
  LorentzVector leadingJet_Reg, subleadingJet_Reg, dijetCandidate_Reg;
  LorentzVector leadingJet_RegKF, subleadingJet_RegKF, dijetCandidate_RegKF;
  LorentzVector diHiggsCandidate, diHiggsCandidate_KF, diHiggsCandidate_Reg,diHiggsCandidate_RegKF;
  vector<int> leadingPhotonID, leadingPhotonISO, subleadingPhotonID, subleadingPhotonISO;
  vector<double> genWeights;
  float leadingJet_bDis, subleadingJet_bDis, jet1PtRes, jet1EtaRes, jet1PhiRes, jet2PtRes, jet2EtaRes, jet2PhiRes;
  float CosThetaStar, leadingPhotonIDMVA, subleadingPhotonIDMVA, DiJetDiPho_DR_2, DiJetDiPho_DR_1, PhoJetMinDr;
  std::map<std::string, int> myTriggerResults;

  double genTotalWeight;
  unsigned int nPromptInDiPhoton;
  int leadingPhotonEVeto, subleadingPhotonEVeto;
  int leadingJet_flavour, subleadingJet_flavour;
  int isSignal, isPhotonCR;
  int nvtx, nvtx2;

  // End of tree objects


  UInt_t nodeFileNum;
  Bool_t nodesOfHH;
  Double_t gen_mHH, gen_ptH1, gen_ptH2, gen_cosTheta, gen_cosTheta2;

  TFile * NRwFile;
  TH2F * NR_Wei_Hists[1507];
  Float_t NRWeights[1507];


  Angles *angles;
  string bTagName;
  edm::InputTag rhoFixedGrid_;

  // Inputs for Constructor:
  // ORDER MATTERS!
  // std::vector<edm::InputTag> inputTagJets_;
  //edm::InputTag vertexes_;
  std::vector<std::string> myTriggers_;
  std::vector<double> phoIDcutEB_;
  std::vector<double> phoIDcutEE_;
  std::vector<double> phoISOcutEB_;
  std::vector<double> phoISOcutEE_;

  std::vector<double> phoISO_nhCorEB_;
  std::vector<double> phoISO_nhCorEE_;
  std::vector<double> phoISO_phCorEB_;
  std::vector<double> phoISO_phCorEE_;
  //edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;

  UInt_t cutFlow_;
  Bool_t useDiPhotons_;
  edm::InputTag diPhotons_;
  UInt_t phoIDtype_;
  Bool_t doBJetRegression_;
  edm::FileInPath bRegFile_;
  Bool_t doNonResWeights_;
};


#endif /* HHbbgg_ANA */
