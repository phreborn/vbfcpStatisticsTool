// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Perform quick conditional fits as needed for profile likelihood
//               scans

#include "TFile.h"
#include "TH1D.h"
#include "TTime.h"
#include "TSystem.h"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"
#include "RooNLLVar.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "TStopwatch.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooSimultaneous.h"
#include "RooRealSumPdf.h"
#include "TTree.h"
#include "TH2.h"
#include "RooGaussian.h"

#include "RooMinimizer.h"
#include "Math/MinimizerOptions.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooGamma.h"
#include "RooPoisson.h"
#include "RooBifurGauss.h"

#include "ExtendedModel.hxx"
#include "ExtendedMinimizer.hxx"

#include "log.hxx"
#include "utils.hxx"

#include <iomanip>
#include <stdlib.h>
#include <list>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <chrono>
#include <regex>

#include "boost/program_options.hpp"
#include "boost/program_options/cmdline.hpp"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/variables_map.hpp"

using namespace std;
using namespace RooFit;
using namespace RooStats;

// _____________________________________________________________________________
// Declarations of functions used in this file

// _____________________________________________________________________________
// Main routine
int main(int argc, char** argv)
{
  TTime thistime = gSystem->Now();

  // Load custom classes
  loadCustom();

  // Model information
  string inFileName      = "path/to/workspace.root";
  string wsName          = "combined";
  string modelConfigName = "ModelConfig";
  string dataName        = "combData";
  string snapshot        = "nominalNuis";
  string folder          = "test";

  // Parameter settings
  string poiName         = "mu";
  string profileName     = "";
  string fixName         = "";
  int fixAllNP           = 0;
  bool makeParameterSnapshots = true;

  // Fit settings
  string minimizerType   = "Minuit2";
  string minimizerAlgo   = "Migrad";
  int defaultStrategy    = 0;
  int binnedLikelihood   = 1;
  int offsetting         = 1;
  int constOpt           = 2;
  double eps             = 1.0;
  int numCPU             = 1;
  double precision       = 0.001;
  bool setInitialError   = false;
  int calls              = -1;
  int iters           = -1;

  // Misc settings
  int fixCache           = 1;
  int fixMulti           = 1;
  int eigendecomposition = 0;
  string loglevel        = "INFO";

  // Bookkeeping

  using namespace boost;
  namespace po = boost::program_options;
  po::options_description desc( "Program options" );
  desc.add_options()
    ( "help           , h"                                                                         , "Print this help message")
    ( "input"         , po::value<string>( &inFileName )->default_value( inFileName )              , "File to run over." )
    ( "poi"           , po::value<string>( &poiName )->default_value( poiName )                    , "POIs to measure." )
    ( "snapshot"      , po::value<string>( &snapshot )->default_value( snapshot )                  , "Initial snapshot." )
    ( "folder"        , po::value<string>( &folder )->default_value( folder )                      , "Output folder." )
    ( "profile"       , po::value<string>( &profileName )->default_value( profileName )            , "Parameters to profile." )
    ( "fix"           , po::value<string>( &fixName )->default_value( fixName )                    , "Parameters to fix." )
    ( "workspace"     , po::value<string>( &wsName )->default_value( wsName )                      , "WS to grab." )
    ( "modelconfig"   , po::value<string>( &modelConfigName )->default_value( modelConfigName )    , "MC to load." )
    ( "data"          , po::value<string>( &dataName )->default_value( dataName )                  , "Data to use." )
    ( "minimizerType" , po::value<string>( &minimizerType )->default_value( minimizerType )        , "Minimizer type." )
    ( "minimizerAlgo" , po::value<string>( &minimizerAlgo )->default_value( minimizerAlgo )        , "Minimizer algorithm." )
    ( "strategy"      , po::value<int>( &defaultStrategy )->default_value( defaultStrategy )       , "Default strategy." )
    ( "numCPU"        , po::value<int>( &numCPU )->default_value( numCPU )                         , "Number of CPUs." )
    ( "binned"        , po::value<int>( &binnedLikelihood )->default_value( binnedLikelihood )     , "Binned likelihood." )
    ( "starfix"       , po::value<int>( &fixCache )->default_value( fixCache )                     , "Fix StarMomentMorph cache." )
    ( "multifix"      , po::value<int>( &fixMulti )->default_value( fixMulti )                     , "Fix MultiPdf level 2." )
    ( "precision"     , po::value<double>( &precision )->default_value( precision )                , "Precision for scan." )
    ( "eps"           , po::value<double>( &eps )->default_value( eps )                            , "Convergence criterium." )
    ( "eigen"         , po::value<int>( &eigendecomposition )->default_value( eigendecomposition ) , "Eigenvalues and vectors." )
    ( "offset"        , po::value<int>( &offsetting )->default_value( offsetting )                 , "Offset likelihood." )
    ( "optimize"      , po::value<int>( &constOpt )->default_value( constOpt )                     , "Optimize constant terms." )
    ( "loglevel"      , po::value<string>( &loglevel )->default_value( loglevel )                  , "POIs to use." )
    ( "fixAllNP"      , po::value<int>( &fixAllNP )->default_value( fixAllNP )                     , "Fix all NP." )
    ( "calls"         , po::value<int>( &calls )->default_value( calls )                           , "Maximum number of function calls." )
    ( "iters"         , po::value<int>( &iters )->default_value( iters )                           , "Maximum number of Minuit iterations." )
    ;

  po::variables_map vm0;

  try {
    po::store( po::command_line_parser( argc, argv ).options( desc ).run(), vm0 );
    po::notify( vm0 );
  }

  catch ( std::exception& ex ) {
    cerr << "Invalid options: " << ex.what() << endl;
    cout << "Use ./a.out --help to print the allowed program options" << endl;
    return -1;
  }

  catch ( ... ) {
    cerr << "Unidentified error parsing program options." << endl;
    return -1;
  }

  // if help, print help
  if ( vm0.count( "help" ) ) {
    cout << "Usage: ./a.out [PROGRAMOPTIONS]\n";
    cout << desc;
    return 0;
  }

  // DEBUG OUTPUT
  // - ERROR
  // - WARNING
  // - INFO
  // - DEBUG
  LOG::ReportingLevel() = LOG::FromString(loglevel);

  // Configuration of minimizer
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerType.c_str(), minimizerAlgo.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
  if (loglevel == "DEBUG") {
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
  } else {
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  }
  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // Parse options
  vector<string> parsed = parseString(poiName, ",");

  // Load the model
  ExtendedModel* model = new ExtendedModel("model", inFileName, wsName,
                                            modelConfigName, dataName, snapshot,
                                            binnedLikelihood, "pdf_", fixCache,
                                            fixMulti);

  RooWorkspace* ws = model->GetWorkspace();
  ModelConfig* mc = model->GetModelConfig();
  RooAbsPdf* pdf = model->GetPdf();
  RooAbsData* data = model->GetData();
  RooArgSet* nuis = model->GetNuisanceParameters();
  RooArgSet* globs = model->GetGlobalObservables();
  RooArgSet* pois = model->GetParametersOfInterest();
  RooArgSet* obs = model->GetObservables();

  std::system(("mkdir -vp root-files/" + folder).c_str());

  if (fixAllNP) {
    model->fixNuisanceParameters();
  }

  if (setInitialError) {
    model->setInitialErrors();
  }

  if (fixName != "") {
    model->fixNuisanceParameters(fixName);
  }

  model->fixParametersOfInterest();

  // Collect POIs
  vector<RooRealVar*> scan_poi_vector;
  RooArgSet scan_poi_set;
  vector<double> poiVals;
  for (size_t i = 0; i < parsed.size(); i++) {
    TString thisName = parsed[i].c_str();
    TString thisVal;
    if (thisName.Contains("[")) {
      assert(thisName.Contains("]"));
      TObjArray* thisNameArray = thisName.Tokenize("[");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      thisVal = ((TObjString*)thisNameArray->At(1))->GetString();
      thisVal.ReplaceAll("]","");
    }

    RooRealVar* poi = (RooRealVar*)ws->var(thisName.Data());
    if (!poi) {
      LOG(logERROR) << "POI: " << thisName.Data() << " doesn't exist!";
      exit(1);
    }
    double val = poi->getVal();

    if (thisVal.IsFloat()) {
      val = thisVal.Atof();
      poi->setVal(val);
      poi->setConstant(1);
    }

    LOG(logINFO) << "Getting POI " << poi->GetName() << " and set value to " << val;

    scan_poi_vector.push_back(poi);
    scan_poi_set.add(*poi);
    poiVals.push_back(val);
  }

  model->profileParameters(profileName);

  // Perform the final fit
  // for (size_t i_eps = 0; i_eps < parsed_eps.size(); i_eps++) {
  // double thisEps = atof(parsed_eps[i_eps].c_str());

  MyTimer timer;
  LOG(logINFO) << "Making ExtendedMinimizer";
  ExtendedMinimizer minimizer("minimizer", pdf, data);
  LOG(logINFO) << "Starting minimization";

  minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()),
                               Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                               Constrain(*nuis), GlobalObservables(*globs),
                               NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                               Precision(precision), Save(), PrintLevel(printLevel), Timer(0),
                               ExtendedMinimizer::MaxFunctionCalls(calls), ExtendedMinimizer::MaxIterations(iters));
  double time = timer.elapsed();
  LOG(logINFO) << "Fitting time: " << setprecision(9) << time << " seconds";

  double minNll = minimizer.GetMinNll();
  LOG(logINFO) << "NLL after minimisation with eps " << eps << ": " << setprecision(15) << minNll;

  // if (i_eps == parsed_eps.size()-1) {
  RooArgSet profiledParameters = minimizer.GetFitResult()->floatParsFinal();
  int status = minimizer.GetFitResult()->status();
  double edm = minimizer.GetFitResult()->edm();

  stringstream filename;
  filename << "root-files/" << folder << "/scan/";
  for (size_t i = 0; i < scan_poi_vector.size(); ++i) {
    filename << scan_poi_vector[i]->GetName() << scan_poi_vector[i]->getVal();
  }
  filename << ".root";
  TFile fout(filename.str().c_str(), "recreate");

  TTree* resultTree = new TTree("result", "result");
  resultTree->SetDirectory(0);

  resultTree->Branch("status", &status);
  resultTree->Branch("edm", &edm);
  resultTree->Branch("nll", &minNll);

  vector<double> fill_poi_vals;
  vector<string> fill_poi_names;
  for (size_t i = 0; i < scan_poi_vector.size(); ++i) {
    double thisVal = scan_poi_vector[i]->getVal();
    string thisName = scan_poi_vector[i]->GetName();
    fill_poi_vals.push_back(thisVal);
    fill_poi_names.push_back(thisName);
  }

  for (size_t i = 0; i < fill_poi_names.size(); ++i) {
    resultTree->Branch(fill_poi_names[i].c_str(), &fill_poi_vals[i]);
  }

  vector<double> fill_nuis_vals;
  vector<string> fill_nuis_names;
  for (RooLinkedListIter it = profiledParameters.iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    double thisVal = v->getVal();
    string thisName = v->GetName();
    fill_nuis_vals.push_back(thisVal);
    fill_nuis_names.push_back(thisName);
  }

  for (size_t i = 0; i < fill_nuis_names.size(); ++i) {
    resultTree->Branch(fill_nuis_names[i].c_str(), &fill_nuis_vals[i]);
  }

  resultTree->Fill();
  resultTree->ResetBranchAddresses();
  resultTree->Write("",TObject::kOverwrite);
  fout.Close();
  // }
  // }

  PrintResourcesUsed(thistime);
}
