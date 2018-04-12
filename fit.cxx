// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Perform quick unconditional fits

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

#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasStyle.h"

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
  bool fixAllNP          = false;
  bool makeParameterSnapshots = false;

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
  vector<double> scanRange = {};
  bool minosScan         = false;
  bool computeZ          = false;
  double altPOI          = 0.0;

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
    ( "fixAllNP"      , po::bool_switch( &fixAllNP )                                               , "Fix all NP." )
    ( "minos"         , po::bool_switch( &minosScan )                                              , "Minos confidence intervals for POI." )
    ( "scan"          , po::value<vector<double>> ( &scanRange )->multitoken()                     , "Range for PLR scan od POI." )
    ( "makeSnapshots" , po::bool_switch( &makeParameterSnapshots )                                 , "Make nominal paramter snapshots." )
    ( "significance"  , po::bool_switch( &computeZ )                                               , "Compute significance wrt alternative hypothesis." )
    ( "alternative"   , po::value<double>( &altPOI )->default_value( altPOI )                      , "POI values for alternative hypothesis." )
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

  if (makeParameterSnapshots) {
    // Save the snapshots of nominal parameters
    LOG(logINFO) << "MakeSnapshots() Saving nominal snapshots.";
    ws->saveSnapshot("nominalGlobs", *mc->GetGlobalObservables());
    ws->saveSnapshot("nominalNuis", *mc->GetNuisanceParameters());
    ws->saveSnapshot("nominalPois", *mc->GetParametersOfInterest());
  }

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
  opt.push_back(Hesse());
  opt.push_back(Save());
  if (minosScan) opt.push_back(ExtendedMinimizer::Scan(scan_poi_set));

  RooLinkedList* cmdList = new RooLinkedList();
  for (auto &o : opt) {
    cmdList->Add((TObject*)&o);
  }

  minimizer.minimize(*cmdList);
  double time = timer.elapsed();
  LOG(logINFO) << "Fitting time: " << setprecision(9) << time << " seconds";

  double minNll = minimizer.GetMinNll();
  LOG(logINFO) << "NLL after minimisation: " << setprecision(15) << minNll;

  // Save fit result
  RooFitResult* result = (minimizer.GetFitResult());

  stringstream filename;
  filename <<  inFileName << ".fitresult";
  TFile fout(filename.str().c_str(), "recreate");
  result->Write("", TObject::kOverwrite);

  if (minosScan) {
    TTree* resultTree = new TTree("result", "result");
    resultTree->SetDirectory(0);
    vector<double> fill_poi_vals;
    vector<string> fill_poi_names;

    for (size_t i = 0; i < scan_poi_vector.size(); ++i) {
      string thisName = scan_poi_vector[i]->GetName();

      double val_pois_hat = scan_poi_vector[i]->getVal();
      double val_pois_hi = scan_poi_vector[i]->getErrorHi();
      double val_pois_lo = scan_poi_vector[i]->getErrorLo();

      fill_poi_names.push_back(thisName + "_hat");
      fill_poi_names.push_back(thisName + "_hi");
      fill_poi_names.push_back(thisName + "_lo");

      fill_poi_vals.push_back(val_pois_hat);
      fill_poi_vals.push_back(val_pois_hi);
      fill_poi_vals.push_back(val_pois_lo);
    }

    for (size_t i = 0; i < fill_poi_names.size(); ++i) {
      resultTree->Branch(fill_poi_names[i].c_str(), &fill_poi_vals[i]);
    }

    resultTree->Fill();
    resultTree->ResetBranchAddresses();
    resultTree->Write("",TObject::kOverwrite);
  }
  fout.Close();
  LOG(logINFO) << "Saved fitresult as " << filename.str();

  // Compute significance (two-sided)
  if (computeZ) {
    for (RooLinkedListIter it = scan_poi_set.iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      v->setVal(altPOI);
      v->setConstant(1);
      LOG(logINFO) << "Fixing POI " << v->GetName() << " to " << altPOI;
    }

    if (minosScan) cmdList->Remove((TObject*)&opt.back());

    minimizer.minimize(*cmdList);

    double minNllAlt = minimizer.GetMinNll();
    LOG(logINFO) << "NLL after minimisation for alternative hypothesis: " << setprecision(15) << minNllAlt;

    int ndf = scan_poi_vector.size();
    double tmu = 2. * (minNllAlt - minNll);
    double p = TMath::Prob(tmu, ndf);
    double z = ROOT::Math::gaussian_quantile(1. - p / 2., 1);

    LOG(logINFO) << "Significance: " << setprecision(15) << z << " sigma";

    // Store significance
    stringstream filenameZ;
    filenameZ <<  inFileName << ".significance";
    TFile foutZ(filenameZ.str().c_str(), "recreate");

    TTree* sigTree = new TTree("significance", "significance");
    sigTree->SetDirectory(0);

    sigTree->Branch("ndf", &ndf);
    sigTree->Branch("tmu", &tmu);
    sigTree->Branch("p", &p);
    sigTree->Branch("z", &z);
    sigTree->Branch("minNllAlt", &minNllAlt);
    sigTree->Branch("minNll", &minNll);

    sigTree->Fill();
    sigTree->ResetBranchAddresses();
    sigTree->Write("",TObject::kOverwrite);
    foutZ.Close();
  }

  if (makeParameterSnapshots) {
    TString filename_with_snapshot(inFileName);
    filename_with_snapshot.ReplaceAll(".root", "_snapshots.root");

    RooArgSet floatingParameters(result->floatParsFinal());
    ws->saveSnapshot("ucmles", floatingParameters);

    // Bring us back to nominal for exporting
    LOG(logINFO) << "MakeSnapshots() Return to nominal parameter values.";
    ws->loadSnapshot("nominalNuis");
    ws->loadSnapshot("nominalGlobs");
    ws->loadSnapshot("nominalPois");

    ws->writeToFile(filename_with_snapshot.Data());
  }

  // Perform profile likelihood scans
  for (auto thisPoi : scan_poi_vector) {
    TString name = thisPoi->GetName();

    if (scanRange.size() != 3) {
      break;
    }

    double lo = scanRange[0];
    double hi = scanRange[1];
    int nbins = scanRange[2];

    LOG(logINFO) << "Perform profile likelihood scan of " << name.Data() << " in the range [" << lo << "," << hi << "] in " << nbins << " bins";

    pair<TGraph*, TGraph*> scan = minimizer.createProfile(ws->var(name.Data()), lo, hi, nbins);

    (scan.first)->SetLineColor(kBlack);
    (scan.first)->SetMarkerColor(kBlack);
    (scan.second)->SetLineColor(kBlack);
    (scan.second)->SetMarkerColor(kBlack);

    LOG(logINFO) << "Plotting profile likelihood scan";

    TCanvas *c = new TCanvas();
    c->SetCanvasSize(800, 800);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0, 0, 0, 0);
    pad->Draw();
    pad->cd();

    double max = (scan.first)->GetYaxis()->GetXmax();
    TH2F *boundary = new TH2F("boundary", "boundary", 1, lo, hi, 1, 0, max);

    boundary->GetXaxis()->SetLabelSize(28);
    boundary->GetYaxis()->SetLabelSize(28);

    boundary->GetXaxis()->SetTitleSize(28);
    boundary->GetYaxis()->SetTitleSize(28);

    boundary->GetXaxis()->SetTitle(name.Data());
    boundary->GetYaxis()->SetTitle("-2 log #Lambda");

    boundary->Draw();
    scan.second->Draw("L SAME");
    scan.first->Draw("P SAME");

    TString scanName = "scan_" + name + ".pdf";
    c->SaveAs(scanName.Data());
  }

  PrintResourcesUsed(thistime);
}
