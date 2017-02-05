// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Compute the impact (read correlation) of a given NP to a set of
//               POIs

#include "TFile.h"
#include "TH1D.h"
#include "TTime.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "Math/MinimizerOptions.h"

#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBifurGauss.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooPoisson.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"

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
  string variable        = "";

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
    ( "parameter"     , po::value<string>( &variable )->default_value( variable )                  , "Parameter to rank." )
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
  vector<string> parsed_profile = parseString(profileName, ",");
  vector<string> parsed_fix = parseString(fixName, ",");

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

  system(("mkdir -vp root-files/" + folder + "/pulls").c_str());

  // Fix nuisance parameters at initial values
  if (fixName != "") {
    for (size_t i = 0; i < parsed_fix.size(); i++) {
      TString thisName = parsed_fix[i].c_str();
      TString thisVal;
      if (thisName.Contains("[")) {
        assert(thisName.Contains("]"));
        TObjArray* thisNameArray = thisName.Tokenize("[");
        thisName = ((TObjString*)thisNameArray->At(0))->GetString();
        thisVal = ((TObjString*)thisNameArray->At(1))->GetString();
        thisVal.ReplaceAll("]","");
      }

      RooRealVar* par = (RooRealVar*)ws->var(thisName.Data());

      if (!par) {
        LOG(logWARNING) << "Parameter " << thisName << " does not exist";
        continue;
      }

      double val = par->getVal();
      if (thisVal.IsFloat()) {
        val = thisVal.Atof();
        par->setVal(val);
      }

      LOG(logINFO) << "Fixing nuisance " << thisName << " at value " << par->getVal();
      par->setConstant(1);
    }
  }

  // By default fix all POIs before floating
  for (RooLinkedListIter it = pois->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    LOG(logINFO) << "Fixing POI " << v->GetName() << " at value " << v->getVal();
    v->removeRange();
    v->setError(0.15);
    v->setConstant(1);
  }

  // Collect POIs
  vector<RooRealVar*> rank_poi_vector;
  for (size_t i = 0; i < parsed.size(); i++) {
    TString thisName = parsed[i].c_str();
    TString thisVal;
    TString range;
    bool useRange    = kFALSE;

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

    RooRealVar* thisPoi = (RooRealVar*)ws->var(thisName.Data());
    if (!thisPoi) {
      LOG(logERROR) << "POI: " << thisName.Data() << " doesn't exist!";
      exit(1);
    }

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
      }
    }

    thisPoi->setError(0.15);
    thisPoi->setConstant(0);

    LOG(logINFO) << "Getting POI " << thisName.Data() << " = " << thisPoi->getVal() << " in [" << thisPoi->getMin() << "," << thisPoi->getMax() << "]";

    rank_poi_vector.push_back(thisPoi);
  }

  // Explicitly float parameters to profile and define ranges and boundaries
  LOG(logINFO) << "Getting POIs to profile";
  for (size_t i = 0; i < parsed_profile.size(); i++) {
    TString thisName = parsed_profile[i];
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

    LOG(logDEBUG) << "Profiling POI: " << thisName.Data();
    RooRealVar* thisPoi = (RooRealVar*)ws->var(thisName);
    if (!thisPoi) {
      LOG(logERROR) << "POI: " << thisPoi->GetName() << " doesn't exist!";
      exit(-1);
    }

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
    thisPoi->setConstant(0);
    LOG(logINFO) << "Profiling POI " << thisName.Data() << " = " << thisPoi->getVal() << " in [" << thisPoi->getMin() << "," << thisPoi->getMax() << "]";
  }

  MyTimer timer_total;

  RooRealVar* nuip = (RooRealVar*)ws->var(variable.c_str());
  nuip->setConstant(0);
  LOG(logINFO) << "Computing error for parameter " << nuip->GetName() << " at " << nuip->getVal();

  MyTimer timer_unconditional;
  LOG(logINFO) << "Making ExtendedMinimizer for unconditional fit";
  ExtendedMinimizer minimizer("minimizer", pdf, data);
  LOG(logINFO) << "Starting minimization";
  minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()),
                               Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                               Constrain(*nuis), GlobalObservables(*globs),
                               NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                               Precision(precision), ExtendedMinimizer::Scan(RooArgSet(*nuip)),
                               PrintLevel(printLevel));
  LOG(logINFO) << "Fitting time: " << setprecision(9) << timer_unconditional.elapsed() << " seconds";
  PrintResourcesUsed(thistime);

  vector<double> pois_hat;
  for (size_t i = 0; i < rank_poi_vector.size(); i++) {
    pois_hat.push_back(rank_poi_vector[i]->getVal());
    LOG(logINFO) << rank_poi_vector[i]->GetName() << " " << rank_poi_vector[i]->getVal();
  }

  // save snapshot
  ws->saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));
  LOG(logINFO) << "Made unconditional snapshot with name tmp_snapshot";

  MyTimer timer_prefit_uncertainty;

  double nuip_hat = nuip->getVal();
  double nuip_errup = nuip->getErrorHi(), nuip_errdown = nuip->getErrorLo();
  double prefitvariation = 1.0;
  RooArgSet* AllConstraints = new RooArgSet();
  // RooArgSet* tmpAllGlobalObservables = (RooArgSet*)mc->GetGlobalObservables();

  if (ws->set(Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", pdf->GetName(), RooNameSet(*data->get()).content()))) {

    // retrieve from cache
    const RooArgSet* constr = ws->set(Form("CACHE_CONSTR_OF_PDF_%s_FOR_OBS_%s", pdf->GetName(), RooNameSet(*data->get()).content()));
    AllConstraints->add(*constr);
    delete constr;
  } else {
    // Load information needed to determine attributes from ModelConfig
    RooAbsPdf* tmpPdf = (RooAbsPdf*)mc->GetPdf();
    RooArgSet* tmpAllNuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
    RooArgSet* tmpAllObservables = (RooArgSet*)mc->GetObservables();

    // Copies, to keep original sets in place after getAllconstraints call
    RooArgSet tmpAllNuisanceParameters2 = *tmpAllNuisanceParameters;
    RooArgSet tmpAllObservables2 = *tmpAllObservables;
    AllConstraints = tmpPdf->getAllConstraints(tmpAllObservables2, tmpAllNuisanceParameters2, kFALSE);
  }

  // Take care of the case where we have a product of constraint terms
  TIterator* ConstraintItrAll = AllConstraints->createIterator();
  RooAbsArg* nextConstraint;
  RooArgSet* tmpAllConstraints = new RooArgSet(AllConstraints->GetName());
  while ((nextConstraint = (RooAbsArg*)ConstraintItrAll->Next())) {
    if (nextConstraint->IsA() == RooProdPdf::Class()) {
      RooArgSet thisComponents;
      FindUniqueProdComponents((RooProdPdf*)nextConstraint, thisComponents);
      tmpAllConstraints->add(thisComponents);
    } else {
      tmpAllConstraints->add(*nextConstraint);
    }
  }

  TIterator* ConstraintItr = tmpAllConstraints->createIterator();
  bool foundConstraint = kFALSE;
  while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next()) && !foundConstraint) {
    if (nextConstraint->dependsOn(*nuip)) {
      foundConstraint = kTRUE;

      // loop over global observables to match nuisance parameter and
      // global observable in case of a constrained nuisnace parameter
      TIterator* GlobsItr = globs->createIterator();
      RooRealVar* nextGlobalObservable;
      bool foundGlobalObservable = kFALSE;
      while ((nextGlobalObservable = (RooRealVar*)GlobsItr->Next()) && !foundGlobalObservable) {
        if (nextConstraint->dependsOn(*nextGlobalObservable)) {
          foundGlobalObservable = kTRUE;

          // find constraint width in case of a Gaussian
          if (nextConstraint->IsA() == RooGaussian::Class()) {
            double oldSigmaVal = 1.0;
            TIterator* ServerItr = nextConstraint->serverIterator();
            RooRealVar* nextServer;
            bool foundSigma = kFALSE;
            while ((nextServer = (RooRealVar*)ServerItr->Next()) && !foundSigma) {
              if (nextServer != nextGlobalObservable && nextServer != nuip) {
                oldSigmaVal = nextServer->getVal();
                foundSigma = kTRUE;
              }
            }

            if (AlmostEqualUlpsAndAbs(oldSigmaVal, 1.0, 0.001, 4)) {
              oldSigmaVal = 1.0;
            }

            if (!foundSigma) {
              LOG(logINFO) << "Sigma for pdf " << nextConstraint->GetName() << " not found. Using 1.0.";
            } else {
              LOG(logINFO) << "Using " << oldSigmaVal << " for sigma of pdf " << nextConstraint->GetName();
            }

            prefitvariation = oldSigmaVal;
          } else if (nextConstraint->IsA() == RooPoisson::Class()) {
            double tau = nextGlobalObservable->getVal();
            LOG(logINFO) << "Found tau " << tau << " of pdf " << nextConstraint->GetName();

            prefitvariation = 1. / sqrt(tau);

            LOG(logINFO) << "Prefit variation is " << prefitvariation;
          }
        }
      }
      delete GlobsItr;
    }
  }
  delete ConstraintItr;

  if (!foundConstraint) {
    LOG(logWARNING) << "Not a constrained parameter. No prefit impact can be determined. Use postfit impact instead.";
  }

  LOG(logINFO) << "Time to find prefit variation: " << setprecision(9) << timer_prefit_uncertainty.elapsed() << " seconds";
  PrintResourcesUsed(thistime);

  // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
  MyTimer timer_postfit_up_impact;
  LOG(logINFO) << "Evaluating effect of moving " << nuip->GetName() << " up by " << nuip_hat << " + " << fabs(nuip_errup);
  ws->loadSnapshot("tmp_snapshot");
  nuip->setVal(nuip_hat + fabs(nuip_errup));
  nuip->setConstant(1);

  minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()),
                               Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                               Constrain(*nuis), GlobalObservables(*globs),
                               NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                               ExtendedMinimizer::ReuseNLL(1), Precision(precision));

  vector<double> pois_up;
  for (size_t i = 0; i < rank_poi_vector.size(); i++) {
    pois_up.push_back(rank_poi_vector[i]->getVal());
  }

  LOG(logINFO) << "Time to find postfit up impact: " << setprecision(9) << timer_postfit_up_impact.elapsed() << " seconds";
  PrintResourcesUsed(thistime);

  MyTimer timer_postfit_down_impact;
  LOG(logINFO) << "Evaluating effect of moving " << nuip->GetName() << " down by " << nuip_hat << " - " << fabs(nuip_errdown);
  ws->loadSnapshot("tmp_snapshot");
  nuip->setVal(nuip_hat - fabs(nuip_errdown));
  nuip->setConstant(1);

  minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()),
                               Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                               Constrain(*nuis), GlobalObservables(*globs),
                               NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                               ExtendedMinimizer::ReuseNLL(1), Precision(precision));

  vector<double> pois_down;
  for (size_t i = 0; i < rank_poi_vector.size(); i++) {
    pois_down.push_back(rank_poi_vector[i]->getVal());
  }

  LOG(logINFO) << "Time to find postfit down impact: " << setprecision(9) << timer_postfit_down_impact.elapsed() << " seconds";
  PrintResourcesUsed(thistime);

  // fix theta at the MLE value +/- prefit uncertainty and minimize again to estimate the change in the POI
  vector<double> pois_nom_up;
  vector<double> pois_nom_down;

  if (foundConstraint) {

    MyTimer timer_prefit_up_impact;
    LOG(logINFO) << "Evaluating effect of moving " << nuip->GetName() << " up by " << nuip_hat << " + " << prefitvariation;
    ws->loadSnapshot("tmp_snapshot");
    nuip->setVal(nuip_hat + prefitvariation);
    nuip->setConstant(1);

    minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()),
                                 Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                                 Constrain(*nuis), GlobalObservables(*globs),
                                 NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                                 ExtendedMinimizer::ReuseNLL(1), Precision(precision));

    for (size_t i = 0; i < rank_poi_vector.size(); i++) {
      pois_nom_up.push_back(rank_poi_vector[i]->getVal());
    }

    LOG(logINFO) << "Time to find prefit up impact: " << setprecision(9) << timer_prefit_up_impact.elapsed() << " seconds";
    PrintResourcesUsed(thistime);

    MyTimer timer_prefit_down_impact;
    LOG(logINFO) << "Evaluating effect of moving " << nuip->GetName() << " down by " << nuip_hat << " - " << prefitvariation;
    ws->loadSnapshot("tmp_snapshot");
    nuip->setVal(nuip_hat-prefitvariation);
    nuip->setConstant(1);

    // ExtendedMinimizer minimizer_prefit_down("minimizer_prefit_down", pdf, data);
    minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()),
                                 Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                                 Constrain(*nuis), GlobalObservables(*globs),
                                 NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                                 ExtendedMinimizer::ReuseNLL(1), Precision(precision));

    for (size_t i = 0; i < rank_poi_vector.size(); i++) {
      pois_nom_down.push_back(rank_poi_vector[i]->getVal());
    }

    LOG(logINFO) << "Time to find prefit down impact: " << setprecision(9) << timer_prefit_down_impact.elapsed() << " seconds";
    PrintResourcesUsed(thistime);
  } else {
    LOG(logWARNING) << "Prefit impact not estimated, instead postfit impact is cloned.";
    for (size_t i = 0; i < rank_poi_vector.size(); i++) {
      pois_nom_up.push_back(pois_up[i]);
    }

    for (size_t i = 0; i < rank_poi_vector.size(); i++) {
      pois_nom_down.push_back(pois_down[i]);
    }
  }

  cout << endl;
  LOG(logINFO) << "Time to perform all fits: " << setprecision(9) << timer_total.elapsed() << " seconds";
  cout << endl;
  LOG(logINFO) << "Unconditional minimum of NP " << variable << ": " << nuip_hat << " +" << fabs(nuip_errup) << " -" << fabs(nuip_errdown);
  LOG(logINFO) << "Prefit uncertainty of NP " << variable << ": " << nuip_hat << " +/- " << prefitvariation;

  for (size_t i = 0; i < rank_poi_vector.size(); i++) {
    cout << endl;
    LOG(logINFO) << "Unconditional minimum of POI " << rank_poi_vector[i]->GetName() << ": " << pois_hat[i];
    LOG(logINFO) << "POI when varying NP up by 1 sigma postfit (prefit): " << pois_up[i] << " (" << pois_nom_up[i] << ")";
    LOG(logINFO) << "POI when varying NP down by 1 sigma postfit (prefit): " << pois_down[i] << " (" << pois_nom_down[i] << ")";
  }

  // store result in root file
  stringstream fileName;
  fileName << "root-files/" << folder << "/pulls/" << variable << ".root";
  TFile fout(fileName.str().c_str(), "recreate");

  TTree* resultTree = new TTree("result", "result");
  resultTree->SetDirectory(0);

  resultTree->Branch("nuisance", &variable);

  resultTree->Branch("nuis_hat", &nuip_hat);
  resultTree->Branch("nuis_hi", &nuip_errup);
  resultTree->Branch("nuis_lo", &nuip_errdown);
  resultTree->Branch("nuis_prefit", &prefitvariation);

  vector<double> fill_poi_vals;
  vector<string> fill_poi_names;

  for (size_t i = 0; i < rank_poi_vector.size(); ++i) {
    string thisName = rank_poi_vector[i]->GetName();

    double val_pois_hat = pois_hat[i];
    double val_pois_up = pois_up[i];
    double val_pois_down = pois_down[i];
    double val_pois_nom_up = pois_nom_up[i];
    double val_pois_nom_down = pois_nom_down[i];

    fill_poi_names.push_back(thisName + "_hat");
    fill_poi_names.push_back(thisName + "_up");
    fill_poi_names.push_back(thisName + "_down");
    fill_poi_names.push_back(thisName + "_up_nom");
    fill_poi_names.push_back(thisName + "_down_nom");

    fill_poi_vals.push_back(val_pois_hat);
    fill_poi_vals.push_back(val_pois_up);
    fill_poi_vals.push_back(val_pois_down);
    fill_poi_vals.push_back(val_pois_nom_up);
    fill_poi_vals.push_back(val_pois_nom_down);
  }

  for (size_t i = 0; i < fill_poi_names.size(); ++i) {
    resultTree->Branch(fill_poi_names[i].c_str(), &fill_poi_vals[i]);
  }

  resultTree->Fill();
  resultTree->ResetBranchAddresses();
  resultTree->Write("",TObject::kOverwrite);
  fout.Close();

  PrintResourcesUsed(thistime);

  return 0;
}
