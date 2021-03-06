#include "../interface/DummyAnalyzer.h"

using namespace std;
using namespace edm;

DummyAnalyzer::DummyAnalyzer(const edm::ParameterSet& cfg, TFileDirectory& fs):
  edm::BasicAnalyzer::BasicAnalyzer(cfg, fs),
  lep_(cfg.getUntrackedParameter<string>("lep")),
  lumiWeight_( cfg.getParameter<double>( "lumiWeight" ) ),
  muons_(cfg.getParameter<edm::InputTag>("muonTag")),
  electrons_(cfg.getParameter<edm::InputTag>("electronTag")),
  photons_(cfg.getParameter<edm::InputTag>("photonTag")),
  jets_(cfg.getParameter<edm::InputTag>("jetTag")),
  myGen_( cfg.getParameter<edm::InputTag>( "genTag" ) ),
  runSample_(cfg.getUntrackedParameter<string>("runSample") )
{
  cout<<"\t DDDDD \t Dummy is constructing you..."<<endl;
  // This is an exmple of how to use TFileService histograms:
  
  hists_["DummyHisto"] = fs.make<TH1F>("DummyPt"  , "pt"  ,  100,  0., 300.);
  // This is the HistoManager class (notice, no underscore):
  hists    = new HistManager(fs.getBareDirectory()->GetFile());
    
  yields_txt_path = open_temp("/tmp", "out_cutlist", fcuts);
  open_temp("/tmp", "synch", fout);

  fout.precision(3); fout.setf(ios::fixed, ios::floatfield);
  red.SetCode(Color::FG_RED);
  green.SetCode(Color::FG_GREEN);
  blue.SetCode(Color::FG_BLUE);
  def.SetCode(Color::FG_DEFAULT);

  
  W0=1;
  totEvents = 0;
  totWeights = 0;
  count_neg = count_pos=0;
  for (Int_t n1=0; n1<nC; n1++){
    nEvents[n1]=0;
    nWeights[n1]=0;
  }

  //Modifier red(FG_RED);
  ///def(Color::FG_DEFAULT);

}

void DummyAnalyzer::beginJob()
{
  //cout<<"\t DDDDD \t Dummy is Beginning the job."<<endl;
}


void DummyAnalyzer::analyze(const edm::EventBase& event)
{
  throw cms::Exception("DUMMY", "I am a dummy analyzer, please don't run me. Inherit from me instead.\
I provide basic function which are useful for all your other analyzers.");
}

void DummyAnalyzer::endJob(UInt_t effBase = 0)
{
  // effBase argument is used to set the base cut for total efficiency calculation (see below)

  cout<<"\t DDDDDDD \t  Dummy is Terminating... Don't sit too close to the monitor - it may explode into your face! **"<<endl;
  fout.close();
  fcuts.close();

  string allCuts[nC];
  string line;
  ifstream myfile(yields_txt_path);
  if (myfile.is_open()){
    while (! myfile.eof() ){
      getline (myfile,line);
      Int_t n = atoi(line.substr(0,2).c_str());
      allCuts[n] = line;
      //cout << n<<"  "<<allCuts[n]<< endl;
    }
    myfile.close();
  }


  cout<<"\t This is "<<green<<runSample_<<def<<" sample"<<endl;

  // Here set the scale factor for the gen weight in order to match the weighted event to a simple sum
  // This is only for visualization (yield table below)
  if (runSample_=="DYJets")
    W0 = 17000;
  else if (runSample_=="QCD")
    W0=1;
 else if (runSample_=="GJets")
    W0=1;
 else if (runSample_=="DiPhoton")
    W0=0.5;
  else if (runSample_=="Graviton")
    W0=1;
  else if (runSample_=="Radion")
    W0=1;
  else W0=1;


  cout<<" ** YIELDS **"<<endl;
  cout<<"n |"<<setw(45)<<" CUT DESCRIPTION \t\t|"<<" events \t"<< "Weight sum |"<<" Tot eff |"<<" cut eff |"<<endl;
  for (UInt_t n=0; n<nC; n++){
    if (n==0)
      cout<<  "0 |"<<setw(45)<<allCuts[n]<<"\t |"<<setw(8)<< nEvents[0]<<"|"<<setw(8)<<ULong64_t(nWeights[0]/W0)<<"|"
	  <<setw(5)<<std::fixed<<std::setprecision(3)<<1.0<<"|"<<1.0<<"|"<<endl;
    else
      cout<<n<<" |"<<setw(45)<<allCuts[n]<<"\t |"<<setw(8)<< nEvents[n]<<"|"<<setw(8)<<ULong64_t(nWeights[n]/W0)<<"|"
	  <<setw(5)<<std::fixed<<std::setprecision(3)<<float(nWeights[n])/nWeights[effBase]<<"|"<<float(nWeights[n])/nWeights[n-1]<<"|"<<endl;
  }



  hists->fill1DHist(-1, "evt_byCut",";cut #;weighted events", nC+1,-1,nC, totWeights, "Counts");
  hists->fill1DHist(-1, "evt_byCut_raw", ";cut #;events",     nC+1,-1,nC, totEvents,  "Counts");

  cout<<"Tot weights = "<<ULong64_t(totWeights)<<endl;
  //cout<<"Tot weights (double) = "<<totWeights<<endl;
  cout<<"Pos weights count "<<count_pos<<"  negatives = "<<count_neg<<"  W0="<<W0<<endl;
  cout<<" ** End of Analyzer **"<<endl;

}


void DummyAnalyzer::CountEvents(Int_t num, string cutName, Double_t wei, ofstream& s)
{
  if (nEvents[num]==0)
    s<<num<<" "<<cutName<<endl;
  nEvents[num]++;
  nWeights[num]+=wei;
  //else edm::LogError("File is not open!");
}


void DummyAnalyzer::FillHistoCounts(Int_t num, Double_t wei)
{
  hists->fill1DHist(num, "evt_byCut",";cut;weighted events", nC+1, -1,nC, wei, "Counts");
  hists->fill1DHist(num, "evt_byCut_raw", ";cut;events",     nC+1, -1,nC, 1, "Counts");
}

std::string DummyAnalyzer::open_temp(std::string path, std::string nm, std::ofstream& f) {
  // This piece of code is taken from:
  // http://stackoverflow.com/questions/499636/how-to-create-a-stdofstream-to-a-temp-file
  // It is usefull for creating temp files in ofstream etc.
  path += "/"+nm+"_XXXXXX";
  std::vector<char> dst_path(path.begin(), path.end());
  dst_path.push_back('\0');

  // See also: http://man7.org/linux/man-pages/man3/mkstemp.3.html
  int fd = mkstemp(&dst_path[0]);
  if(fd != -1) {
    path.assign(dst_path.begin(), dst_path.end() - 1);
    f.open(path.c_str(),
	   std::ios_base::trunc | std::ios_base::out);
    close(fd);

    //cout<<"\t DDDD Open temp OK.  Path="<<path<<endl;
  }

  return path;
}
