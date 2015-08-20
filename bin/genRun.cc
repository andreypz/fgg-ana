#include "../interface/GenAnalyzer.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "PhysicsTools/UtilAlgos/interface/FWLiteAnalyzerWrapper.h"

typedef fwlite::AnalyzerWrapper<GenAnalyzer> WrappedFWLiteGenAnalyzer;

int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
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
  WrappedFWLiteGenAnalyzer ana( *( builder.processDesc()->getProcessPSet() ), std::string( "GenAnalyzer" ), std::string( "test" ) );

  ana.beginJob();
  ana.analyze();
  ana.endJob();
  return 0;
}
