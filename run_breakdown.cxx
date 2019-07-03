// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2018-05-02
// Description : Perform uncertainty breakdowns

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
#include "cpp-text-table/TextTable.h"

#include "atlasrootstyle/AtlasUtils.h"
#include "atlasrootstyle/AtlasLabels.h"
#include "atlasrootstyle/AtlasStyle.h"

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

#include "yaml-cpp/yaml.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

// _____________________________________________________________________________
// Declarations of functions used in this file
vector<string> getClassification(string name, YAML::Node classes);

// _____________________________________________________________________________
// Main routine
int main(int argc, char** argv)
{
  TTime thistime = gSystem->Now();

  // Load custom classes
  loadCustom();

  // ATLAS style
  SetAtlasStyle();

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

  // Misc settings
  int fixCache           = 1;
  int fixMulti           = 1;
  string loglevel        = "INFO";

  // Bookkeeping
  string classification  = "classification.yaml";
  string category        = "";
  bool subtractFromTotal = false;
  bool doIndividual = false;

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
    ( "offset"        , po::value<int>( &offsetting )->default_value( offsetting )                 , "Offset likelihood." )
    ( "optimize"      , po::value<int>( &constOpt )->default_value( constOpt )                     , "Optimize constant terms." )
    ( "loglevel"      , po::value<string>( &loglevel )->default_value( loglevel )                  , "Control verbosity." )
    ( "classification", po::value<string>( &classification )->default_value( classification )      , "Definition of uncertainty categories." )
    ( "category"      , po::value<string>( &category )->default_value( category )                        , "Specific category which should be submitted" )
    ( "subtractFromTotal", po::bool_switch( &subtractFromTotal )                                   , "Subtract uncertainties from total." )
    ( "doIndividual"  , po::bool_switch( &doIndividual )                                           , "Compute uncertainty for individual sources." )
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
    TString thisName = parsed[i];
    TString range;
    TString boundary;
    int sign = 0;

    bool useRange    = kFALSE;
    bool useBoundary = kFALSE;

    // Get ranges
    if (thisName.Contains("[")) {
      assert(thisName.Contains("]"));
      TObjArray* thisNameArray = thisName.Tokenize("[");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      range = ((TObjString*)thisNameArray->At(1))->GetString();
      range.ReplaceAll("]","");
      assert(range.Contains(":"));
      useRange = kTRUE;
    }

    // Get sign
    if (thisName.Contains("+")) {
      thisName.ReplaceAll("+",">0");
    } else if (thisName.Contains("-")) {
      thisName.ReplaceAll("-","<0");
    }

    // Get boundaries
    if (thisName.Contains(">")) {
      TObjArray* thisNameArray = thisName.Tokenize(">");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      boundary = ((TObjString*)thisNameArray->At(1))->GetString();
      sign = +1;
      useBoundary = kTRUE;
    } else if (thisName.Contains("<")) {
      TObjArray* thisNameArray = thisName.Tokenize("<");
      thisName = ((TObjString*)thisNameArray->At(0))->GetString();
      boundary = ((TObjString*)thisNameArray->At(1))->GetString();
      sign = -1;
      useBoundary = kTRUE;
    }

    RooRealVar* thisPoi = (RooRealVar*)ws->var(thisName);
    if (!thisPoi) {
      LOG(logERROR) << "POI: " << thisPoi->GetName() << " doesn't exist!";
      exit(-1);
    }

    thisPoi->removeRange();

    if (useRange) {
      double origVal = thisPoi->getVal();
      TObjArray* rangeArray = range.Tokenize(":");
      TString s_lo = ((TObjString*)rangeArray->At(0))->GetString();
      TString s_hi = ((TObjString*)rangeArray->At(1))->GetString();
      double lo = atof(s_lo.Data());
      double hi = atof(s_hi.Data());
      thisPoi->setRange(lo, hi);
      if ((origVal < lo) || (origVal > hi)) {
        double newVal = (hi - lo) / 2;
        thisPoi->setVal(newVal);
        LOG(logERROR) << "Setting: " << thisPoi->GetName() << " to value " << newVal;
      }
    }

    if (useBoundary) {
      double tmpBoundary = atof(boundary.Data());
      double origVal = thisPoi->getVal();
      double forigVal = fabs(thisPoi->getVal());
      bool boundaryIsZero = AlmostEqualUlpsAndAbs(tmpBoundary, 0.0, 0.0001, 4);

      if (sign > 0) {
        thisPoi->setMin(tmpBoundary);
        if (origVal < tmpBoundary) {
          thisPoi->setVal(tmpBoundary);
        }
        if (boundaryIsZero && origVal < 0) {
          thisPoi->setVal(forigVal);
        }
      } else if (sign < 0) {
        thisPoi->setMax(tmpBoundary);
        if (origVal > tmpBoundary) {
          thisPoi->setVal(tmpBoundary);
        }
        if (boundaryIsZero && origVal > 0) {
          thisPoi->setVal(-forigVal);
        }
      }
    }

    thisPoi->setError(0.2);
    thisPoi->setConstant(0);
    double val = thisPoi->getVal();
    LOG(logINFO) << "Getting POI " << thisPoi->GetName() << " and set value to " << val;

    scan_poi_vector.push_back(thisPoi);
    scan_poi_set.add(*thisPoi);
    poiVals.push_back(val);
  }

  model->profileParameters(profileName);

  // Add other POIs to nuisance set
  RooArgSet* nuisOtherPoi = new RooArgSet();
  nuisOtherPoi->add(*nuis);

  for (RooLinkedListIter it = pois->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    if ( std::find(scan_poi_vector.begin(), scan_poi_vector.end(), v) == scan_poi_vector.end() ) {
      nuisOtherPoi->add(*v);
    }
  }

  nuisOtherPoi->Print("v");

  // Load the classification file
  LOG(logINFO) << "Load classification from " << classification;
  YAML::Node classes = YAML::LoadFile(classification);

  // Assign classes
  LOG(logINFO) << "Assign nuisance parameter classes";

  map<string, vector<string>> nuisance_assignments;

  for (RooLinkedListIter it = nuisOtherPoi->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    string name = v->GetName();
    vector<string> nuisanceClasses = getClassification(name, classes);

    if (nuisanceClasses.size() > 0) {
      for (auto nuisance_class : nuisanceClasses) {
        nuisance_assignments[nuisance_class].push_back(name);
      }
    } else {
      nuisance_assignments["None"].push_back(name);
      LOG(logWARNING) << "Parameter " << name << " assigned to no categories!";
    }
  }

  ofstream myfile ("class_details.yaml");
  for (auto nuisance_class : nuisance_assignments) {
    myfile << "\n";
    myfile << nuisance_class.first << ":\n";
    for (auto name : nuisance_class.second) {
      myfile << "  - " << name << "\n";
    }
  }
  myfile.close();

  // Unconditional fit
  LOG(logINFO) << "Run unconditional fit";

  MyTimer timer;
  ExtendedMinimizer minimizer("minimizer", pdf, data);

  vector<RooCmdArg> opt;
  opt.push_back(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()));
  opt.push_back(Strategy(defaultStrategy));
  opt.push_back(ExtendedMinimizer::Eps(eps));
  opt.push_back(Constrain(*nuis));
  opt.push_back(GlobalObservables(*globs));
  opt.push_back(NumCPU(numCPU, 3));
  opt.push_back(Offset(offsetting));
  opt.push_back(Optimize(constOpt));
  opt.push_back(Precision(precision));

  RooLinkedList* cmdList_init = new RooLinkedList();
  for (auto &o : opt) {
    cmdList_init->Add((TObject*)&o);
  }

  minimizer.minimize(*cmdList_init);
  double time = timer.elapsed();
  LOG(logINFO) << "Fitting time: " << setprecision(9) << time << " seconds";
  double mle = scan_poi_vector[0]->getVal();

  ws->saveSnapshot("tmp_shot", *mc->GetPdf()->getParameters(data));

  // Find total uncertainty
  LOG(logINFO) << "Find total uncertainty";

  opt.push_back(ExtendedMinimizer::ReuseNLL(1));
  opt.push_back(ExtendedMinimizer::Scan(scan_poi_set));

  RooLinkedList* cmdList = new RooLinkedList();
  for (auto &o : opt) {
    cmdList->Add((TObject*)&o);
  }

  ws->loadSnapshot("tmp_shot");
  minimizer.minimize(*cmdList);
  double err_hi = scan_poi_vector[0]->getErrorHi(), err_lo = scan_poi_vector[0]->getErrorLo();

  // Use pulls to define sign of the uncertainty
  LOG(logINFO) << "Define sign of the uncertainty";

  map<string, string> signs_up;
  map<string, string> signs_down;

  for (RooLinkedListIter it = nuisOtherPoi->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    string name = v->GetName();
    double pull = v->getVal();
    if (pull < 0) {
      signs_up[name] = "+";
      signs_down[name] = "-";
    } else {
      signs_up[name] = "-";
      signs_down[name] = "+";
    }
  }

  // Find the statistical uncertainty
  LOG(logINFO) << "Find statistical uncertainty";

  for (RooLinkedListIter it = nuisOtherPoi->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    string name = v->GetName();

    if (subtractFromTotal) {
      vector<string> vec = nuisance_assignments["Normalisation"];
      if (std::find(vec.begin(), vec.end(), name) == vec.end()) {
        LOG(logINFO) << "Fixing parameter " << name;
        v->setConstant(true);
      }
    } else {
      vector<string> vec = nuisance_assignments["Normalisation"];
      if (std::find(vec.begin(), vec.end(), name) == vec.end()) {
        LOG(logINFO) << "Fixing parameter " << name;
        v->setConstant(true);

      }
    }
  }

  ws->saveSnapshot("tmp_shot2", *mc->GetPdf()->getParameters(data));

  minimizer.minimize(*cmdList);
  double stat_err_hi = scan_poi_vector[0]->getErrorHi(), stat_err_lo = scan_poi_vector[0]->getErrorLo();

  // Individual nuisance parameters
  double quad_hi = 0, quad_lo = 0;
  map<string, double> ind_err_hi;
  map<string, double> ind_err_lo;
  set<pair<double, string>>  set_sys;

  if (doIndividual) {
    LOG(logINFO) << "Evaluate individual uncertainties";

    if (subtractFromTotal) {
      ws->loadSnapshot("tmp_shot");
    } else {
      ws->loadSnapshot("tmp_shot2");
    }

    for (RooLinkedListIter it = nuisOtherPoi->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      string name = v->GetName();
      LOG(logINFO) << "Getting impact from " << name;

      if (subtractFromTotal) {
        ws->loadSnapshot("tmp_shot");
        v->setConstant(true);
      } else {
        ws->loadSnapshot("tmp_shot2");
        v->setConstant(false);
      }

      minimizer.minimize(*cmdList);
      ind_err_hi[name] = scan_poi_vector[0]->getErrorHi();
      ind_err_lo[name] = scan_poi_vector[0]->getErrorLo();

      ind_err_hi[name] = subtract_error(ind_err_hi[name], stat_err_hi);
      ind_err_lo[name] = subtract_error(ind_err_lo[name], stat_err_lo);

      if (ind_err_hi[name] != ind_err_hi[name] || ind_err_hi[name] < 0.001) {
        ind_err_hi[name] = 0;
      }
      if (ind_err_lo[name] != ind_err_lo[name] || ind_err_lo[name] < 0.001) {
        ind_err_lo[name] = 0;
      }

      set_sys.insert(make_pair(sqrt((1 + ind_err_hi[name]) * (1 + ind_err_lo[name])) - 1, name));

      quad_hi += ind_err_hi[name] * ind_err_hi[name];
      quad_lo += ind_err_lo[name] * ind_err_lo[name];

      if (subtractFromTotal) {
        v->setConstant(false);
      } else {
        v->setConstant(true);
      }
    }

    quad_hi = sqrt(quad_hi);
    quad_lo = sqrt(quad_lo);
  }

  // Getting uncertainties for each group
  LOG(logINFO) << "Get uncertainty for each group";

  map<string, double> all_err_hi;
  map<string, double> all_err_lo;

  for (auto class_itr : nuisance_assignments) {
    string nuisance_class = class_itr.first;
    if (category != "" && !(category == nuisance_class || nuisance_class == "Normalisation" || nuisance_class == "TemplateStatistics")){
         continue;
    }
    vector<string> vec = class_itr.second;
    LOG(logINFO) << "On class " << nuisance_class;

    ws->loadSnapshot("tmp_shot");

    for (RooLinkedListIter it = nuisOtherPoi->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      string name = v->GetName();


      if (subtractFromTotal) {
        if (std::find(vec.begin(), vec.end(), name) != vec.end()) {
          LOG(logINFO) << "Fixing parameter " << name;
          v->setConstant(true);
        }
      } else {
        // normalisation should always float
        vector<string> vec_stat = nuisance_assignments["Normalisation"];
        if (std::find(vec_stat.begin(), vec_stat.end(), name) != vec_stat.end()) {
          LOG(logINFO) << "Floating parameter " << name;
          continue;
        }

        if (std::find(vec.begin(), vec.end(), name) != vec.end()) {
          LOG(logINFO) << "Floating parameter " << name;
          continue;
        }

        LOG(logINFO) << "Fixing parameter " << name;
        v->setConstant(true);
      }
    }

    minimizer.minimize(*cmdList);
    all_err_hi[nuisance_class] = scan_poi_vector[0]->getErrorHi();
    all_err_lo[nuisance_class] = scan_poi_vector[0]->getErrorLo();
  }

  // Find data statistics and control region statistics uncertainties
  LOG(logINFO) << "Find data statistics and control region statistics uncertainties";

  ws->loadSnapshot("tmp_shot");
  for (RooLinkedListIter it = nuisOtherPoi->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    string name = v->GetName();
    LOG(logINFO) << "Fixing parameter " << name;
    v->setConstant(true);
  }

  ws->saveSnapshot("tmp_shot3", *mc->GetPdf()->getParameters(data));

  minimizer.minimize(*cmdList);
  double sr_err_hi = scan_poi_vector[0]->getErrorHi();
  double sr_err_lo = scan_poi_vector[0]->getErrorLo();

  set<pair<double, string>>  set_stat;
  map<string, double> cr_err_hi;
  map<string, double> cr_err_lo;
  double cr_sig_err_hi = 0;
  double cr_sig_err_lo = 0;

  vector<string> poi_names;
  for (auto v : scan_poi_vector) {
    poi_names.push_back(v->GetName());
  }

  for (RooLinkedListIter it = nuisOtherPoi->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    string name = v->GetName();

    vector<string> vec_stat = nuisance_assignments["Normalisation"];
    if (std::find(vec_stat.begin(), vec_stat.end(), name) == vec_stat.end()) {
      LOG(logINFO) << "Skipping parameter " << name;
      continue;
    }

    ws->loadSnapshot("tmp_shot3");
    LOG(logINFO) << "Floating parameter " << name;
    v->setConstant(false);

    minimizer.minimize(*cmdList);
    cr_err_hi[name] = scan_poi_vector[0]->getErrorHi();
    cr_err_lo[name] = scan_poi_vector[0]->getErrorLo();

    cr_err_hi[name] = subtract_error(cr_err_hi[name], sr_err_hi);
    cr_err_lo[name] = subtract_error(cr_err_lo[name], sr_err_lo);

    if (cr_err_hi[name] != cr_err_hi[name] || cr_err_hi[name] < 0.001) {
      cr_err_hi[name] = 0;
    }
    if (cr_err_lo[name] != cr_err_lo[name] || cr_err_lo[name] < 0.001) {
      cr_err_lo[name] = 0;
    }

    if (std::find(poi_names.begin(), poi_names.end(), name) != poi_names.end()) {
      cr_sig_err_hi = cr_err_hi[name];
      cr_sig_err_lo = cr_err_lo[name];
    }

    set_stat.insert(make_pair(sqrt((1 + cr_err_hi[name]) * (1 + cr_err_lo[name])) - 1, name));

    v->setConstant(true);
  }

  // Compute components
  LOG(logINFO) << "Compute component uncertainties";

  double sys_err_hi = subtract_error(err_hi, stat_err_hi);
  double sys_err_lo = subtract_error(err_lo, stat_err_lo);

  double mc_stat_err_hi = 0;
  double mc_stat_err_lo = 0;

  if (subtractFromTotal) {
    mc_stat_err_hi = subtract_error(err_hi, all_err_hi["TemplateStatistics"]);
    mc_stat_err_lo = subtract_error(err_lo, all_err_lo["TemplateStatistics"]);
  } else {
    mc_stat_err_hi = subtract_error(all_err_hi["TemplateStatistics"], stat_err_hi);
    mc_stat_err_lo = subtract_error(all_err_lo["TemplateStatistics"], stat_err_lo);
  }

  double sys_nomc_err_hi = subtract_error(sys_err_hi, mc_stat_err_hi);
  double sys_nomc_err_lo = subtract_error(sys_err_lo, mc_stat_err_lo);

  string hi, lo;

  set<string> set_sysNames;
  if (doIndividual) {
    TextTable t('-', '|', '+');

    t.add("nuisance");
    t.add("hi");
    t.add("lo");
    t.endOfRow();

    for (auto it : set_sys) {
      string name = it.second;

      t.add(name);
      hi = to_string(ind_err_hi[name]);
      lo = to_string(ind_err_lo[name]);
      hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
      lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
      t.add(signs_up[name] + hi);
      t.add(signs_down[name] + lo);
      t.endOfRow();

      set_sysNames.insert(name);
    }

    t.setAlignment(1, TextTable::Alignment::RIGHT);
    t.setAlignment(2, TextTable::Alignment::RIGHT);
    std::cout << t;
  }

  TextTable t_stat('-', '|', '+');

  t_stat.add("nuisance (stat)");
  t_stat.add("hi");
  t_stat.add("lo");
  t_stat.endOfRow();

  for (auto it : set_stat) {
    string name = it.second;

    t_stat.add(name);
    hi = to_string(cr_err_hi[name]);
    lo = to_string(cr_err_lo[name]);
    hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
    lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
    t_stat.add("+" + hi);
    t_stat.add("-" + lo);
    t_stat.endOfRow();

    set_sysNames.insert(name);
  }

  t_stat.setAlignment(1, TextTable::Alignment::RIGHT);
  t_stat.setAlignment(2, TextTable::Alignment::RIGHT);
  std::cout << t_stat;

  map<string, double> all_err_hi_comp;
  map<string, double> all_err_lo_comp;

  for (auto class_itr : nuisance_assignments) {
    string nuisance_class = class_itr.first;
    if (category != "" && !(category == nuisance_class || nuisance_class == "Normalisation" || nuisance_class == "TemplateStatistics")){
         continue;
    }

    if (subtractFromTotal) {
      all_err_hi_comp[nuisance_class] = subtract_error(err_hi, all_err_hi[nuisance_class]);
      all_err_lo_comp[nuisance_class] = subtract_error(err_lo, all_err_lo[nuisance_class]);
    } else {
      all_err_hi_comp[nuisance_class] = subtract_error(all_err_hi[nuisance_class], stat_err_hi);
      all_err_lo_comp[nuisance_class] = subtract_error(all_err_lo[nuisance_class], stat_err_lo);
    }
  }

  TextTable t_breakdown('-', '|', '+');

  t_breakdown.add("group");
  t_breakdown.add("hi");
  t_breakdown.add("lo");
  t_breakdown.endOfRow();

  t_breakdown.add("total");
  hi = to_string(PDGrounding(mle, err_hi, 1).second);
  lo = to_string(PDGrounding(mle, err_lo, 1).second);
  hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
  lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
  t_breakdown.add("+" + hi);
  t_breakdown.add(      lo);
  t_breakdown.endOfRow();

  t_breakdown.add("data statistics");
  hi = to_string(PDGrounding(mle, stat_err_hi, 1).second);
  lo = to_string(PDGrounding(mle, stat_err_lo, 1).second);
  hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
  lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
  t_breakdown.add("+" + hi);
  t_breakdown.add(      lo);
  t_breakdown.endOfRow();

  t_breakdown.add("template statistics");
  hi = to_string(PDGrounding(mle, mc_stat_err_hi, 1).second);
  lo = to_string(PDGrounding(mle, mc_stat_err_lo, 1).second);
  hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
  lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
  t_breakdown.add("+" + hi);
  t_breakdown.add(      lo);
  t_breakdown.endOfRow();

  t_breakdown.add("systematics");
  hi = to_string(PDGrounding(mle, sys_err_hi, 1).second);
  lo = to_string(PDGrounding(mle, sys_err_lo, 1).second);
  hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
  lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
  t_breakdown.add("+" + hi);
  t_breakdown.add(      lo);
  t_breakdown.endOfRow();

  t_breakdown.add("systematics w/o template statistics");
  hi = to_string(PDGrounding(mle, sys_nomc_err_hi, 1).second);
  lo = to_string(PDGrounding(mle, sys_nomc_err_lo, 1).second);
  hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
  lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
  t_breakdown.add("+" + hi);
  t_breakdown.add(      lo);
  t_breakdown.endOfRow();

  for (auto class_itr : nuisance_assignments) {
    string nuisance_class = class_itr.first;
    if (category != "" && !(category == nuisance_class || nuisance_class == "Normalisation" || nuisance_class == "TemplateStatistics")){
         continue;
    }

    if (nuisance_class == "Normalisation") {
      continue;
    }

    if (nuisance_class == "TemplateStatistics") {
      continue;
    }

    t_breakdown.add(nuisance_class);
    hi = to_string(PDGrounding(mle, all_err_hi_comp[nuisance_class], 1).second);
    lo = to_string(PDGrounding(mle, all_err_lo_comp[nuisance_class], 1).second);
    hi.erase(hi.find_last_not_of('0') + 1, std::string::npos);
    lo.erase(lo.find_last_not_of('0') + 1, std::string::npos);
    t_breakdown.add("+" + hi);
    t_breakdown.add(      lo);
    t_breakdown.endOfRow();
  }

  t_breakdown.setAlignment(1, TextTable::Alignment::RIGHT);
  t_breakdown.setAlignment(2, TextTable::Alignment::RIGHT);
  std::cout << t_breakdown;

  PrintResourcesUsed(thistime);
}

// _____________________________________________________________________________
vector<string> getClassification(string name, YAML::Node classes)
{
  vector<string> all_classes;
  map<string, vector<string>> veto_classes;

  for (const auto node : classes) {
    YAML::Node identifier = node.first;
    YAML::Node negators = node.second["negators"];
    YAML::Node selectors = node.second["selectors"];

    for (const auto expression : negators) {
      veto_classes[identifier.as<string>()].push_back(expression.as<string>());
    }

    for (const auto expression : selectors) {
      std::regex select_regex(expression.as<string>());
      if (std::regex_match(name, select_regex)) {
        all_classes.push_back(identifier.as<string>());
      }
    }
  }

  vector<string> matched_classes = all_classes;
  set<string> s(matched_classes.begin(), matched_classes.end());
  matched_classes.assign(s.begin(), s.end());

  for (auto nuisance_class : all_classes) {
    for (auto check_class : veto_classes[nuisance_class]) {
      if ((std::find(all_classes.begin(), all_classes.end(), check_class) != all_classes.end()) && (std::find(matched_classes.begin(), matched_classes.end(), nuisance_class) != matched_classes.end())) {
        matched_classes.erase(std::remove(matched_classes.begin(), matched_classes.end(), nuisance_class), matched_classes.end());
      }
    }
  }

  return matched_classes;
}
