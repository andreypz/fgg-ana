#include "../interface/HHbbggAnalyzer.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "PhysicsTools/UtilAlgos/interface/FWLiteAnalyzerWrapper.h"

typedef fwlite::AnalyzerWrapper<HHbbggAnalyzer> WrappedFWLiteHHbbggAnalyzer;

int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  //AutoLibraryLoader::enable();
  
  gSystem->Load( "libDataFormatsFWLite.so" );
  gSystem->Load( "libDataFormatsPatCandidates.so" );
  gSystem->Load( "libflashggDataFormats.so" );

  for (int a=0; a<argc; a++){
    std::cout<<a<<"  "<<argv[a]<<std::endl;
  }
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  PythonProcessDesc builder( argv[1], argc, argv );
  //WrappedFWLiteHHbbggAnalyzer ana( *( builder.processDesc()->getProcessPSet() ), std::string( "HHbbggAnalyzer" ) );
  WrappedFWLiteHHbbggAnalyzer ana( *( builder.processDesc()->getProcessPSet() ), std::string( "HHbbggAnalyzer" ), std::string( "fsDir" ) );

  ana.beginJob();
  ana.analyze();
  ana.endJob();
  return 0;
}
