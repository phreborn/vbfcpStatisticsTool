// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2013-04-13
// Description : Draw pulls of nuisance parameters, rank by importance

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <list>
#include <map>
#include <math.h>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTime.h"
#include "TTree.h"
#include "TTreeFormula.h"

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

#include "log.hxx"
#include "utils.hxx"

// #include "atlasstyle-00-03-05/AtlasUtils.h"
// #include "atlasstyle-00-03-05/AtlasLabels.h"
// #include "atlasstyle-00-03-05/AtlasStyle.h"

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
int main(int argc, char **argv) {
  TTime thistime = gSystem->Now();
  // SetAtlasStyle();

  string cardName = "";
  string overlayCard = "";
  string poiname = "mu";

  int firstParameter = 1;
  int showTopParameters = -1;
  double showHighImpact = 0.0;

  double scale_poi = 10.0;
  double scale_theta = 10.0;

  bool drawParamNames = true;
  bool drawPostfitImpactBand = true;
  bool drawPrefitImpactBand = true;
  bool useRelativeImpact = false;
  bool visualizeCorrelation = true;
  bool showOneSided = true;

  bool rankNuis = true;
  string loglevel = "INFO";

  string hex_color_standardband = "#ff0000";
  string hex_color_standardband_ol = "#950098";
  string hex_color_pulls = "#9d9d9d";
  string hex_color_prefit = "#000155";
  string hex_color_prefit_ol = "#027734";
  string hex_color_postfit = "#0094ff";
  string hex_color_postfit_ol = "#44ff94";

  // Bookkeeping

  using namespace boost;
  namespace po = boost::program_options;
  po::options_description desc( "Program options" );
  desc.add_options()
    ( "help           , h"                                                                         , "Print this help message")
    ( "input"         , po::value<string>( &cardName )->default_value( cardName )                  , "Folder to run over." )
    ( "overlay"       , po::value<string>( &overlayCard )->default_value( overlayCard )            , "Folder to overlay." )
    ( "poi"           , po::value<string>( &poiname )->default_value( poiname )                    , "POI in ranking." )
    ( "first"         , po::value<int>( &firstParameter )->default_value( firstParameter )         , "First parameter." )
    ( "top"           , po::value<int>( &showTopParameters )->default_value( showTopParameters )   , "Top X parameters." )
    ( "threshold"     , po::value<double>( &showHighImpact )->default_value( showHighImpact )      , "Threshold on impact." )
    ( "scale_poi"     , po::value<double>( &scale_poi )->default_value( scale_poi )                , "Scaling of impact axis." )
    ( "scale_theta"   , po::value<double>( &scale_theta )->default_value( scale_theta )            , "Scaling of pull axis." )
    ( "names"         , po::value<bool>( &drawParamNames )->default_value( drawParamNames )        , "Nuisance names." )
    ( "postfit"       , po::value<bool>( &drawPostfitImpactBand )->default_value( drawPostfitImpactBand ) , "Postfit impact." )
    ( "prefit"        , po::value<bool>( &drawPrefitImpactBand )->default_value( drawPrefitImpactBand )   , "Prefit impact." )
    ( "relative"      , po::value<bool>( &useRelativeImpact )->default_value( useRelativeImpact )         , "Show relative variation." )
    ( "correlation"   , po::value<bool>( &visualizeCorrelation )->default_value( visualizeCorrelation )   , "Visualize sign of correlation." )
    ( "onesided"      , po::value<bool>( &showOneSided )->default_value( showOneSided )                   , "Show one sided variations." )
    ( "rank"          , po::value<bool>( &rankNuis )->default_value( rankNuis )                           , "Order nuisances by impact." )
    ( "loglevel"      , po::value<string>( &loglevel )->default_value( loglevel )                         , "POIs to use." )
    ;

  po::variables_map vm0;

  try {
    po::store( po::command_line_parser( argc, argv ).options( desc ).run(), vm0 );
    po::notify( vm0 );
  }

  catch (std::exception &ex) {
    cerr << "Invalid options: " << ex.what() << endl;
    cout << "Use ./a.out --help to print the allowed program options" << endl;
    return -1;
  }

  catch (...) {
    cerr << "Unidentified error parsing program options." << endl;
    return -1;
  }

  // if help, print help
  if (vm0.count("help")) {
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

  if (useRelativeImpact) {
    LOG(logWARNING) << "Can't use relative impact at the moment. Switch to absolute variation";
    useRelativeImpact = false;
  }

  // Modifications of inputs
  firstParameter -= 1;
  Color_t color_standardband = TColor::GetColor(hex_color_standardband.c_str());
  Color_t color_standardband_ol = TColor::GetColor(hex_color_standardband_ol.c_str());
  Color_t color_pulls = TColor::GetColor(hex_color_pulls.c_str());
  Color_t color_prefit = TColor::GetColor(hex_color_prefit.c_str());
  Color_t color_prefit_ol = TColor::GetColor(hex_color_prefit_ol.c_str());
  Color_t color_postfit = TColor::GetColor(hex_color_postfit.c_str());
  Color_t color_postfit_ol = TColor::GetColor(hex_color_postfit_ol.c_str());

  // Make the canvas for the ranking plot
  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 2000);
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.0, 1.0, 1.0, 0);

  if (drawParamNames) {
    pad1->SetLeftMargin(0.35);
  } else {
    pad1->SetLeftMargin(0.05);
  }

  pad1->SetRightMargin(0.05);
  pad1->SetBottomMargin(0.09);
  pad1->SetTopMargin(0.09);

  pad1->Draw();
  pad1->cd();

  // Run the actual plotting
  vector<double> points_nuis;
  vector<double> val;
  vector<double> up;
  vector<double> down;
  vector<double> poi_hat;
  vector<double> poi_up;
  vector<double> poi_down;
  vector<double> poi_nom_up;
  vector<double> poi_nom_down;
  vector<string> labels;
  map<string, double> prefit_variations_from_file;

  system(("hadd -f " + cardName + "/tmp.root " + cardName + "/*.root").c_str());

  TFile *file = TFile::Open((cardName + "/tmp.root").c_str());
  TTree *tree = (TTree *)(file->Get("result"));

  TTreeFormula *formula_nuis_hat = new TTreeFormula("nuis_hat", "nuis_hat", tree);
  TTreeFormula *formula_nuis_hi = new TTreeFormula("nuis_hi", "nuis_hi", tree);
  TTreeFormula *formula_nuis_lo = new TTreeFormula("nuis_lo", "nuis_lo", tree);
  TTreeFormula *formula_nuis_prefit = new TTreeFormula("nuis_prefit", "nuis_prefit", tree);

  TTreeFormula *formula_poi_hat = new TTreeFormula((poiname + "_hat").c_str(), (poiname + "_hat").c_str(), tree);
  TTreeFormula *formula_poi_up = new TTreeFormula((poiname + "_up").c_str(), (poiname + "_up").c_str(), tree);
  TTreeFormula *formula_poi_down = new TTreeFormula((poiname + "_down").c_str(), (poiname + "_down").c_str(), tree);
  TTreeFormula *formula_poi_nom_up = new TTreeFormula((poiname + "_up_nom").c_str(), (poiname + "_up_nom").c_str(), tree);
  TTreeFormula *formula_poi_nom_down = new TTreeFormula((poiname + "_down_nom").c_str(), (poiname + "_down_nom").c_str(), tree);

  string *label_ptr = new string();
  tree->SetBranchAddress("nuisance", &label_ptr);

  int n_events = tree->GetEntries();

  for (int j = 0; j < n_events; ++j) {
    tree->GetEntry(j);

    double val_nuis_hat = 1.0 * formula_nuis_hat->EvalInstance();
    double val_nuis_hi = 1.0 * formula_nuis_hi->EvalInstance();
    double val_nuis_lo = 1.0 * formula_nuis_lo->EvalInstance();
    double val_nuis_prefit = 1.0 * formula_nuis_prefit->EvalInstance();

    points_nuis.push_back(j);
    val.push_back(val_nuis_hat);
    up.push_back(val_nuis_hi);
    down.push_back(val_nuis_lo);
    prefit_variations_from_file[*label_ptr] = val_nuis_prefit;

    double val_poi_hat = 1.0 * formula_poi_hat->EvalInstance();
    double val_poi_up = 1.0 * formula_poi_up->EvalInstance();
    double val_poi_down = 1.0 * formula_poi_down->EvalInstance();
    double val_poi_nom_up = 1.0 * formula_poi_nom_up->EvalInstance();
    double val_poi_nom_down = 1.0 * formula_poi_nom_down->EvalInstance();

    poi_hat.push_back(val_poi_hat);
    poi_up.push_back(val_poi_up);
    poi_down.push_back(val_poi_down);
    poi_nom_up.push_back(val_poi_nom_up);
    poi_nom_down.push_back(val_poi_nom_down);

    labels.push_back(*label_ptr);
  }

  int nrNuis = points_nuis.size();
  vector<double> points_nuis2 = points_nuis;
  points_nuis2.resize(2 * nrNuis);

  for (int i = 0; i < nrNuis; i++) {
    points_nuis[i] = i + 0.5;
    points_nuis2[2 * i] = i + 0.25;
    points_nuis2[2 * i + 1] = i + 0.75;
  }

  vector<double> poi_hat2 = poi_hat;
  vector<double> poi_up2 = poi_up;
  vector<double> poi_down2 = poi_down;
  vector<double> poi_nom_up2 = poi_nom_up;
  vector<double> poi_nom_down2 = poi_nom_down;

  poi_hat.resize(2 * nrNuis);
  poi_up.resize(2 * nrNuis);
  poi_down.resize(2 * nrNuis);
  poi_nom_up.resize(2 * nrNuis);
  poi_nom_down.resize(2 * nrNuis);

  vector<double> points_nuis_ol;
  vector<double> val_ol;
  vector<double> up_ol;
  vector<double> down_ol;
  vector<double> poi_hat_ol;
  vector<double> poi_up_ol;
  vector<double> poi_down_ol;
  vector<double> poi_nom_up_ol;
  vector<double> poi_nom_down_ol;
  vector<string> labels_ol;
  map<string, double> prefit_variations_from_file_ol;

  if (overlayCard != "") {
    system(("hadd -f " + overlayCard + "/tmp_ol.root " + overlayCard + "/*.root").c_str());

    TFile *file_ol = TFile::Open((overlayCard + "/tmp_ol.root").c_str());
    TTree *tree_ol = (TTree *)(file_ol->Get("result"));

    TTreeFormula *formula_nuis_hat_ol = new TTreeFormula("nuis_hat", "nuis_hat", tree_ol);
    TTreeFormula *formula_nuis_hi_ol = new TTreeFormula("nuis_hi", "nuis_hi", tree_ol);
    TTreeFormula *formula_nuis_lo_ol = new TTreeFormula("nuis_lo", "nuis_lo", tree_ol);
    TTreeFormula *formula_nuis_prefit_ol = new TTreeFormula("nuis_prefit", "nuis_prefit", tree_ol);

    TTreeFormula *formula_poi_hat_ol = new TTreeFormula((poiname + "_hat").c_str(), (poiname + "_hat").c_str(), tree_ol);
    TTreeFormula *formula_poi_up_ol = new TTreeFormula((poiname + "_up").c_str(), (poiname + "_up").c_str(), tree_ol);
    TTreeFormula *formula_poi_down_ol = new TTreeFormula((poiname + "_down").c_str(), (poiname + "_down").c_str(), tree_ol);
    TTreeFormula *formula_poi_nom_up_ol = new TTreeFormula((poiname + "_up_nom").c_str(), (poiname + "_up_nom").c_str(), tree_ol);
    TTreeFormula *formula_poi_nom_down_ol = new TTreeFormula((poiname + "_down_nom").c_str(), (poiname + "_down_nom").c_str(), tree_ol);

    string *label_ptr_ol = new string();
    tree_ol->SetBranchAddress("nuisance", &label_ptr_ol);

    int n_events_ol = tree_ol->GetEntries();

    for (int j = 0; j < n_events_ol; ++j) {
      tree_ol->GetEntry(j);

      double val_nuis_hat = 1.0 * formula_nuis_hat_ol->EvalInstance();
      double val_nuis_hi = 1.0 * formula_nuis_hi_ol->EvalInstance();
      double val_nuis_lo = 1.0 * formula_nuis_lo_ol->EvalInstance();
      double val_nuis_prefit = 1.0 * formula_nuis_prefit_ol->EvalInstance();

      points_nuis_ol.push_back(j);
      val_ol.push_back(val_nuis_hat);
      up_ol.push_back(val_nuis_hi);
      down_ol.push_back(val_nuis_lo);
      prefit_variations_from_file_ol[*label_ptr] = val_nuis_prefit;

      double val_poi_hat = 1.0 * formula_poi_hat_ol->EvalInstance();
      double val_poi_up = 1.0 * formula_poi_up_ol->EvalInstance();
      double val_poi_down = 1.0 * formula_poi_down_ol->EvalInstance();
      double val_poi_nom_up = 1.0 * formula_poi_nom_up_ol->EvalInstance();
      double val_poi_nom_down = 1.0 * formula_poi_nom_down_ol->EvalInstance();

      poi_hat_ol.push_back(val_poi_hat);
      poi_up_ol.push_back(val_poi_up);
      poi_down_ol.push_back(val_poi_down);
      poi_nom_up_ol.push_back(val_poi_nom_up);
      poi_nom_down_ol.push_back(val_poi_nom_down);

      labels_ol.push_back(*label_ptr);
    }
  }

  int nrNuis_ol = points_nuis_ol.size();
  vector<double> points_nuis2_ol = points_nuis_ol;
  points_nuis2_ol.resize(2 * nrNuis_ol);

  for (int i = 0; i < nrNuis_ol; i++) {
    points_nuis_ol[i] = i + 0.25;
    points_nuis2_ol[2 * i] = i + 0.125;
    points_nuis2_ol[2 * i + 1] = i + 0.375;
  }

  vector<double> poi_hat2_ol = poi_hat_ol;
  vector<double> poi_up2_ol = poi_up_ol;
  vector<double> poi_down2_ol = poi_down_ol;
  vector<double> poi_nom_up2_ol = poi_nom_up_ol;
  vector<double> poi_nom_down2_ol = poi_nom_down_ol;

  poi_hat_ol.resize(2 * nrNuis_ol);
  poi_up_ol.resize(2 * nrNuis_ol);
  poi_down_ol.resize(2 * nrNuis_ol);
  poi_nom_up_ol.resize(2 * nrNuis_ol);
  poi_nom_down_ol.resize(2 * nrNuis_ol);

  if (overlayCard != "") {
    for (int i = 0; i < nrNuis; i++) {
      points_nuis[i] = i + 0.75;
      points_nuis2[2 * i] = i + 0.625;
      points_nuis2[2 * i + 1] = i + 0.875;
    }
  }

  // set correct values for the poi
  vector<string> antiCorrelated;
  vector<string> antiCorrelated_ol;
  vector<string> antiCorrelatedNom;
  vector<string> antiCorrelatedNom_ol;

  for (int i = 0; i < nrNuis; i++) {
    val[i] *= scale_theta;

    poi_up[2 * i] = poi_up2[i] - poi_hat2[i];
    poi_up[2 * i + 1] = poi_up2[i] - poi_hat2[i];
    poi_down[2 * i] = poi_down2[i] - poi_hat2[i];
    poi_down[2 * i + 1] = poi_down2[i] - poi_hat2[i];

    poi_nom_up[2 * i] = poi_nom_up2[i] - poi_hat2[i];
    poi_nom_up[2 * i + 1] = poi_nom_up2[i] - poi_hat2[i];
    poi_nom_down[2 * i] = poi_nom_down2[i] - poi_hat2[i];
    poi_nom_down[2 * i + 1] = poi_nom_down2[i] - poi_hat2[i];

    if (poi_up[2 * i] < 0 && poi_down[2 * i] > 0) {
      antiCorrelated.push_back(labels[i]);
      if (visualizeCorrelation) {
        swap(poi_up[2 * i], poi_down[2 * i]);
        swap(poi_up[2 * i + 1], poi_down[2 * i + 1]);
      }
    }

    if (poi_nom_up[2 * i] < 0 && poi_nom_down[2 * i] > 0) {
      antiCorrelatedNom.push_back(labels[i]);
      if (visualizeCorrelation) {
        swap(poi_nom_up[2 * i], poi_nom_down[2 * i]);
        swap(poi_nom_up[2 * i + 1], poi_nom_down[2 * i + 1]);
      }
    }

    if (!visualizeCorrelation || !showOneSided) {
      if (poi_up[2 * i] < 0)
        swap(poi_up[2 * i], poi_down[2 * i]);
      if (poi_up[2 * i + 1] < 0)
        swap(poi_up[2 * i + 1], poi_down[2 * i + 1]);
      if (poi_nom_up[2 * i] < 0)
        swap(poi_nom_up[2 * i], poi_nom_down[2 * i]);
      if (poi_nom_up[2 * i + 1] < 0)
        swap(poi_nom_up[2 * i + 1], poi_nom_down[2 * i + 1]);
    }

    if (showOneSided) {
      // both positive or negative
      if (poi_up[2 * i] > 0 && poi_down[2 * i] > 0) {
        poi_up[2 * i] = poi_down[2 * i + 1];
        poi_down[2 * i] = 0;
        poi_down[2 * i + 1] = 0;
      } else if (poi_up[2 * i] < 0 && poi_down[2 * i] < 0) {
        poi_down[2 * i + 1] = poi_up[2 * i];
        poi_up[2 * i] = 0;
        poi_up[2 * i + 1] = 0;
      }

      if (poi_nom_up[2 * i] > 0 && poi_nom_down[2 * i] > 0) {
        poi_nom_up[2 * i] = poi_nom_down[2 * i + 1];
        poi_nom_down[2 * i] = 0;
        poi_nom_down[2 * i + 1] = 0;
      } else if (poi_nom_up[2 * i] < 0 && poi_nom_down[2 * i] < 0) {
        poi_nom_down[2 * i + 1] = poi_nom_up[2 * i];
        poi_nom_up[2 * i] = 0;
        poi_nom_up[2 * i + 1] = 0;
      }
    }

    poi_up[2 * i] = fabs(poi_up[2 * i]);
    poi_up[2 * i + 1] = fabs(poi_up[2 * i + 1]);
    poi_down[2 * i] = fabs(poi_down[2 * i]);
    poi_down[2 * i + 1] = fabs(poi_down[2 * i + 1]);

    poi_nom_up[2 * i] = fabs(poi_nom_up[2 * i]);
    poi_nom_up[2 * i + 1] = fabs(poi_nom_up[2 * i + 1]);
    poi_nom_down[2 * i] = fabs(poi_nom_down[2 * i]);
    poi_nom_down[2 * i + 1] = fabs(poi_nom_down[2 * i + 1]);

    poi_hat[2 * i] = 0;
    poi_hat[2 * i + 1] = 0;
  }

  for (int i = 0; i < nrNuis_ol; i++) {
    val_ol[i] *= scale_theta;

    poi_up_ol[2 * i] = poi_up2_ol[i] - poi_hat2_ol[i];
    poi_up_ol[2 * i + 1] = poi_up2_ol[i] - poi_hat2_ol[i];
    poi_down_ol[2 * i] = poi_down2_ol[i] - poi_hat2_ol[i];
    poi_down_ol[2 * i + 1] = poi_down2_ol[i] - poi_hat2_ol[i];

    poi_nom_up_ol[2 * i] = poi_nom_up2_ol[i] - poi_hat2_ol[i];
    poi_nom_up_ol[2 * i + 1] = poi_nom_up2_ol[i] - poi_hat2_ol[i];
    poi_nom_down_ol[2 * i] = poi_nom_down2_ol[i] - poi_hat2_ol[i];
    poi_nom_down_ol[2 * i + 1] = poi_nom_down2_ol[i] - poi_hat2_ol[i];

    if (poi_up_ol[2 * i] < 0 && poi_down_ol[2 * i] > 0) {
      antiCorrelated_ol.push_back(labels_ol[i]);
      if (visualizeCorrelation) {
        swap(poi_up_ol[2 * i], poi_down_ol[2 * i]);
        swap(poi_up_ol[2 * i + 1], poi_down_ol[2 * i + 1]);
      }
    }

    if (poi_nom_up_ol[2 * i] < 0 && poi_nom_down_ol[2 * i] > 0) {
      antiCorrelatedNom_ol.push_back(labels_ol[i]);
      if (visualizeCorrelation) {
        swap(poi_nom_up_ol[2 * i], poi_nom_down_ol[2 * i]);
        swap(poi_nom_up_ol[2 * i + 1], poi_nom_down_ol[2 * i + 1]);
      }
    }

    if (!visualizeCorrelation || !showOneSided) {
      if (poi_up_ol[2 * i] < 0)
        swap(poi_up_ol[2 * i], poi_down_ol[2 * i]);
      if (poi_up_ol[2 * i + 1] < 0)
        swap(poi_up_ol[2 * i + 1], poi_down_ol[2 * i + 1]);
      if (poi_nom_up_ol[2 * i] < 0)
        swap(poi_nom_up_ol[2 * i], poi_nom_down_ol[2 * i]);
      if (poi_nom_up_ol[2 * i + 1] < 0)
        swap(poi_nom_up_ol[2 * i + 1], poi_nom_down_ol[2 * i + 1]);
    }

    if (showOneSided) {
      // both positive or negative
      if (poi_up_ol[2 * i] > 0 && poi_down_ol[2 * i] > 0) {
        poi_up_ol[2 * i] = poi_down_ol[2 * i + 1];
        poi_down_ol[2 * i] = 0;
        poi_down_ol[2 * i + 1] = 0;
      } else if (poi_up_ol[2 * i] < 0 && poi_down_ol[2 * i] < 0) {
        poi_down_ol[2 * i + 1] = poi_up_ol[2 * i];
        poi_up_ol[2 * i] = 0;
        poi_up_ol[2 * i + 1] = 0;
      }

      if (poi_nom_up_ol[2 * i] > 0 && poi_nom_down_ol[2 * i] > 0) {
        poi_nom_up_ol[2 * i] = poi_nom_down_ol[2 * i + 1];
        poi_nom_down_ol[2 * i] = 0;
        poi_nom_down_ol[2 * i + 1] = 0;
      } else if (poi_nom_up_ol[2 * i] < 0 && poi_nom_down_ol[2 * i] < 0) {
        poi_nom_down_ol[2 * i + 1] = poi_nom_up_ol[2 * i];
        poi_nom_up_ol[2 * i] = 0;
        poi_nom_up_ol[2 * i + 1] = 0;
      }
    }

    poi_up_ol[2 * i] = fabs(poi_up_ol[2 * i]);
    poi_up_ol[2 * i + 1] = fabs(poi_up_ol[2 * i + 1]);
    poi_down_ol[2 * i] = fabs(poi_down_ol[2 * i]);
    poi_down_ol[2 * i + 1] = fabs(poi_down_ol[2 * i + 1]);

    poi_nom_up_ol[2 * i] = fabs(poi_nom_up_ol[2 * i]);
    poi_nom_up_ol[2 * i + 1] = fabs(poi_nom_up_ol[2 * i + 1]);
    poi_nom_down_ol[2 * i] = fabs(poi_nom_down_ol[2 * i]);
    poi_nom_down_ol[2 * i + 1] = fabs(poi_nom_down_ol[2 * i + 1]);

    poi_hat_ol[2 * i] = 0;
    poi_hat_ol[2 * i + 1] = 0;
  }

  // find maximal error due to a single nuisance parameter
  // TODO: for overlay as well? maybe get rid of this...
  double max_poi = 0.;
  for (int i = 0; i < nrNuis * 2; ++i) {
    if (poi_up[i] > max_poi)
      max_poi = poi_up[i];
    if (poi_down[i] > max_poi)
      max_poi = poi_down[i];
  }

  // Total uncertainty of POI
  // TODO hard coded at the moment
  double sigma_tot_hi = 1.0;
  double sigma_tot_lo = 1.0;
  double sigma_tot_ol_hi;
  double sigma_tot_ol_lo;

  if (overlayCard != "") {
    sigma_tot_ol_hi = 1.0;
    sigma_tot_ol_lo = 1.0;
  }

  // dump everything in maps
  map<string, vector<double>> nuis_map;
  for (int i = 0; i < nrNuis; i++) {
    string index = labels[i];
    nuis_map[index].push_back(val[i]);
    nuis_map[index].push_back(up[i]);
    nuis_map[index].push_back(down[i]);
    nuis_map[index].push_back(poi_hat[2 * i]);
    nuis_map[index].push_back(poi_hat[2 * i + 1]);
    nuis_map[index].push_back(poi_up[2 * i]);
    nuis_map[index].push_back(poi_up[2 * i + 1]);
    nuis_map[index].push_back(poi_down[2 * i]);
    nuis_map[index].push_back(poi_down[2 * i + 1]);
    nuis_map[index].push_back(poi_nom_up[2 * i]);
    nuis_map[index].push_back(poi_nom_up[2 * i + 1]);
    nuis_map[index].push_back(poi_nom_down[2 * i]);
    nuis_map[index].push_back(poi_nom_down[2 * i + 1]);
  }

  map<string, vector<double>> nuis_map_ol;
  if (overlayCard != "") {
    for (int i = 0; i < nrNuis_ol; i++) {
      string index = labels_ol[i];
      nuis_map_ol[index].push_back(val_ol[i]);
      nuis_map_ol[index].push_back(up_ol[i]);
      nuis_map_ol[index].push_back(down_ol[i]);
      nuis_map_ol[index].push_back(poi_hat_ol[2 * i]);
      nuis_map_ol[index].push_back(poi_hat_ol[2 * i + 1]);
      nuis_map_ol[index].push_back(poi_up_ol[2 * i]);
      nuis_map_ol[index].push_back(poi_up_ol[2 * i + 1]);
      nuis_map_ol[index].push_back(poi_down_ol[2 * i]);
      nuis_map_ol[index].push_back(poi_down_ol[2 * i + 1]);
      nuis_map_ol[index].push_back(poi_nom_up_ol[2 * i]);
      nuis_map_ol[index].push_back(poi_nom_up_ol[2 * i + 1]);
      nuis_map_ol[index].push_back(poi_nom_down_ol[2 * i]);
      nuis_map_ol[index].push_back(poi_nom_down_ol[2 * i + 1]);
    }
  }

  // check that we have in both maps the same keys
  if (overlayCard != "") {
    std::vector<string> tmpstring;

    int counter = nrNuis;
    int counter_ol = nrNuis_ol;

    for (int i = 0; i < nrNuis; i++) {
      bool found = 0;
      for (int ii = 0; ii < nrNuis_ol; ii++) {
        if (labels[i] == labels_ol[ii])
          found = 1;
      }
      if (!found) {
        nuis_map_ol[labels[i]].push_back(-999);
        val_ol.push_back(-999);
        nuis_map_ol[labels[i]].push_back(0);
        up_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        down_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_hat_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_hat_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_up_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_up_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_down_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_down_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_nom_up_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_nom_up_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_nom_down_ol.push_back(0);
        nuis_map_ol[labels[i]].push_back(0);
        poi_nom_down_ol.push_back(0);
        labels_ol.push_back(labels[i]);
        points_nuis_ol.push_back(counter_ol + 0.25);
        points_nuis2_ol.push_back(counter_ol + 0.125);
        points_nuis2_ol.push_back(counter_ol + 0.375);
        counter_ol++;
      }
    }

    for (int i = 0; i < nrNuis_ol; i++) {
      bool found = 0;
      for (int ii = 0; ii < nrNuis; ii++) {
        if (labels_ol[i] == labels[ii])
          found = 1;
      }
      if (!found) {
        nuis_map[labels_ol[i]].push_back(-999);
        val.push_back(-999);
        nuis_map[labels_ol[i]].push_back(0);
        up.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        down.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_hat.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_hat.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_up.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_up.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_down.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_down.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_nom_up.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_nom_up.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_nom_down.push_back(0);
        nuis_map[labels_ol[i]].push_back(0);
        poi_nom_down.push_back(0);
        labels.push_back(labels_ol[i]);
        points_nuis.push_back(counter + 0.75);
        points_nuis2.push_back(counter + 0.625);
        points_nuis2.push_back(counter + 0.875);
        counter++;
      }
    }
  }

  // Getting the vectors back
  nrNuis = labels.size();
  if (overlayCard != "")
    nrNuis_ol = labels_ol.size();

  for (int i = 0; i < nrNuis - 1; i++) {
    for (int j = 0; j < nrNuis - 1 - i; j++) {
      if (strcmp(labels[i].c_str(), labels[i + 1].c_str())) {
        swap(labels[j], labels[j + 1]);
        if (overlayCard != "")
          labels_ol[j + 1] = labels[j + 1];
        if (overlayCard != "")
          labels_ol[j] = labels[j];
      }
    }
  }

  for (int i = 0; i < nrNuis; i++) {
    val[i] = nuis_map[labels[i]][0];
    up[i] = nuis_map[labels[i]][1];
    down[i] = nuis_map[labels[i]][2];
    poi_hat[2 * i] = nuis_map[labels[i]][3];
    poi_hat[2 * i + 1] = nuis_map[labels[i]][4];
    poi_up[2 * i] = nuis_map[labels[i]][5];
    poi_up[2 * i + 1] = nuis_map[labels[i]][6];
    poi_down[2 * i] = nuis_map[labels[i]][7];
    poi_down[2 * i + 1] = nuis_map[labels[i]][8];
    poi_nom_up[2 * i] = nuis_map[labels[i]][9];
    poi_nom_up[2 * i + 1] = nuis_map[labels[i]][10];
    poi_nom_down[2 * i] = nuis_map[labels[i]][11];
    poi_nom_down[2 * i + 1] = nuis_map[labels[i]][12];
  }

  for (int i = 0; i < nrNuis_ol; i++) {
    val_ol[i] = nuis_map_ol[labels_ol[i]][0];
    up_ol[i] = nuis_map_ol[labels_ol[i]][1];
    down_ol[i] = nuis_map_ol[labels_ol[i]][2];
    poi_hat_ol[2 * i] = nuis_map_ol[labels_ol[i]][3];
    poi_hat_ol[2 * i + 1] = nuis_map_ol[labels_ol[i]][4];
    poi_up_ol[2 * i] = nuis_map_ol[labels_ol[i]][5];
    poi_up_ol[2 * i + 1] = nuis_map_ol[labels_ol[i]][6];
    poi_down_ol[2 * i] = nuis_map_ol[labels_ol[i]][7];
    poi_down_ol[2 * i + 1] = nuis_map_ol[labels_ol[i]][8];
    poi_nom_up_ol[2 * i] = nuis_map_ol[labels_ol[i]][9];
    poi_nom_up_ol[2 * i + 1] = nuis_map_ol[labels_ol[i]][10];
    poi_nom_down_ol[2 * i] = nuis_map_ol[labels_ol[i]][11];
    poi_nom_down_ol[2 * i + 1] = nuis_map_ol[labels_ol[i]][12];
  }

  if (overlayCard != "") {
    for (int i = 0; i < nrNuis; i++) {
      cout << "Matches: " << labels[i] << " " << labels_ol[i] << endl;
    }
  }

  // sort poi values by variation size
  if (rankNuis) {
    for (int i = 0; i < nrNuis - 1; i++) {
      for (int j = 0; j < nrNuis - 1 - i; j++) {
        double sumi = 0.0;
        double sumii = 0.0;
        int counti = 0;
        int countii = 0;

        if (poi_up[2 * j] > 0) {
          sumi += poi_up[2 * j];
          counti++;
        }
        if (poi_down[2 * j] > 0) {
          sumi += poi_down[2 * j];
          counti++;
        }
        if (poi_up[2 * j + 1] > 0) {
          sumi += poi_up[2 * j + 1];
          counti++;
        }
        if (poi_down[2 * j + 1] > 0) {
          sumi += poi_down[2 * j + 1];
          counti++;
        }

        if (poi_up[2 * j + 2] > 0) {
          sumii += poi_up[2 * j + 2];
          countii++;
        }
        if (poi_down[2 * j + 2] > 0) {
          sumii += poi_down[2 * j + 2];
          countii++;
        }
        if (poi_up[2 * j + 3] > 0) {
          sumii += poi_up[2 * j + 3];
          countii++;
        }
        if (poi_down[2 * j + 3] > 0) {
          sumii += poi_down[2 * j + 3];
          countii++;
        }

        if (counti > 0)
          sumi /= counti;
        if (countii > 0)
          sumii /= countii;

        if (sumi > sumii) {

          // swap postfit poi
          swap(poi_up[2 * j], poi_up[2 * j + 2]);
          swap(poi_up[2 * j + 1], poi_up[2 * j + 3]);
          swap(poi_down[2 * j], poi_down[2 * j + 2]);
          swap(poi_down[2 * j + 1], poi_down[2 * j + 3]);
          if (overlayCard != "") {
            swap(poi_up_ol[2 * j], poi_up_ol[2 * j + 2]);
            swap(poi_up_ol[2 * j + 1], poi_up_ol[2 * j + 3]);
            swap(poi_down_ol[2 * j], poi_down_ol[2 * j + 2]);
            swap(poi_down_ol[2 * j + 1], poi_down_ol[2 * j + 3]);
          }

          // swap prefit poi
          swap(poi_nom_up[2 * j], poi_nom_up[2 * j + 2]);
          swap(poi_nom_up[2 * j + 1], poi_nom_up[2 * j + 3]);
          swap(poi_nom_down[2 * j], poi_nom_down[2 * j + 2]);
          swap(poi_nom_down[2 * j + 1], poi_nom_down[2 * j + 3]);
          if (overlayCard != "") {
            swap(poi_nom_up_ol[2 * j], poi_nom_up_ol[2 * j + 2]);
            swap(poi_nom_up_ol[2 * j + 1], poi_nom_up_ol[2 * j + 3]);
            swap(poi_nom_down_ol[2 * j], poi_nom_down_ol[2 * j + 2]);
            swap(poi_nom_down_ol[2 * j + 1], poi_nom_down_ol[2 * j + 3]);
          }

          // swap pulls
          swap(up[j], up[j + 1]);
          swap(down[j], down[j + 1]);
          swap(val[j], val[j + 1]);
          if (overlayCard != "") {
            swap(up_ol[j], up_ol[j + 1]);
            swap(down_ol[j], down_ol[j + 1]);
            swap(val_ol[j], val_ol[j + 1]);
          }

          // swap names
          swap(labels[j], labels[j + 1]);
          if (overlayCard != "") {
            swap(labels_ol[j], labels_ol[j + 1]);
          }
        }
      }
    }
  }

  // make the 1 sigma boxes
  vector<double> boxup;
  vector<double> boxdown;
  vector<double> cenup;
  vector<double> cendown;

  vector<double> boxup_ol;
  vector<double> boxdown_ol;
  vector<double> cenup_ol;
  vector<double> cendown_ol;

  for (int i = 0; i < nrNuis; i++) {
    double prefitVariation = prefit_variations_from_file[labels[i]];

    up[i] /= prefitVariation;
    down[i] /= prefitVariation;

    boxup.push_back(prefitVariation * scale_theta / prefitVariation);
    boxdown.push_back(prefitVariation * scale_theta / prefitVariation);
    double height = 0.25;

    if (overlayCard != "")
      height = 0.125;

    cenup.push_back(height);
    cenup.push_back(height);
    cendown.push_back(height);
    cendown.push_back(height);
  }

  for (int i = 0; i < nrNuis_ol; i++) {
    double prefitVariation = prefit_variations_from_file_ol[labels[i]];
    boxup_ol.push_back(prefitVariation * scale_theta);
    boxdown_ol.push_back(prefitVariation * scale_theta);
    double height = 0.125;
    cenup_ol.push_back(height);
    cenup_ol.push_back(height);
    cendown_ol.push_back(height);
    cendown_ol.push_back(height);
  }

  // make the final arrays for plotting, in particular remove parameters
  int nrNuis2remove = 0;
  for (int i = 0; i < nrNuis; i++) {
    cout << "DEBUG::Checking " << labels[i] << " " << fabs(poi_down[2 * i] - poi_hat[2 * i]) << " " << fabs(poi_up[2 * i] - poi_hat[2 * i]) << endl;

    if ((fabs(poi_down[2 * i]) + fabs(poi_up[2 * i])) / (sigma_tot_lo + sigma_tot_hi) < showHighImpact) {
      cout << "WARNING::Removing " << labels[i] << ". Below threshold." << endl;
      nrNuis2remove++;
    }
  }

  cout << "firstParameter is " << firstParameter << endl;

  reverse(boxdown_ol.begin(), boxdown_ol.end());
  reverse(boxup_ol.begin(), boxup_ol.end());
  reverse(cendown_ol.begin(), cendown_ol.end());
  reverse(cenup_ol.begin(), cenup_ol.end());
  reverse(down_ol.begin(), down_ol.end());
  reverse(labels_ol.begin(), labels_ol.end());
  reverse(poi_down_ol.begin(), poi_down_ol.end());
  reverse(poi_hat_ol.begin(), poi_hat_ol.end());
  reverse(poi_nom_down_ol.begin(), poi_nom_down_ol.end());
  reverse(poi_nom_up_ol.begin(), poi_nom_up_ol.end());
  reverse(poi_up_ol.begin(), poi_up_ol.end());
  reverse(points_nuis_ol.begin(), points_nuis_ol.end());
  reverse(up_ol.begin(), up_ol.end());
  reverse(val_ol.begin(), val_ol.end());
  reverse(boxdown.begin(), boxdown.end());
  reverse(boxup.begin(), boxup.end());
  reverse(cendown.begin(), cendown.end());
  reverse(cenup.begin(), cenup.end());
  reverse(down.begin(), down.end());
  reverse(labels.begin(), labels.end());
  reverse(poi_down.begin(), poi_down.end());
  reverse(poi_hat.begin(), poi_hat.end());
  reverse(poi_nom_down.begin(), poi_nom_down.end());
  reverse(poi_nom_up.begin(), poi_nom_up.end());
  reverse(poi_up.begin(), poi_up.end());
  reverse(points_nuis.begin(), points_nuis.end());
  reverse(up.begin(), up.end());
  reverse(val.begin(), val.end());

  labels.erase(labels.begin(), labels.begin() + firstParameter);
  // points_nuis.erase(points_nuis.end() - firstParameter, points_nuis.end());

  if (overlayCard != "") {
    labels_ol.erase(labels_ol.begin(), labels_ol.begin() + firstParameter);
    // points_nuis_ol.erase(points_nuis_ol.end() - firstParameter, points_nuis_ol.end());
  }

  val.erase(val.begin(), val.begin() + firstParameter);
  down.erase(down.begin(), down.begin() + firstParameter);
  up.erase(up.begin(), up.begin() + firstParameter);

  if (overlayCard != "") {
    val_ol.erase(val_ol.begin(), val_ol.begin() + firstParameter);
    down_ol.erase(down_ol.begin(), down_ol.begin() + firstParameter);
    up_ol.erase(up_ol.begin(), up_ol.begin() + firstParameter);
  }

  poi_hat.erase(poi_hat.begin(), poi_hat.begin() + 2 * firstParameter);
  poi_down.erase(poi_down.begin(), poi_down.begin() + 2 * firstParameter);
  poi_up.erase(poi_up.begin(), poi_up.begin() + 2 * firstParameter);

  if (overlayCard != "") {
    poi_hat_ol.erase(poi_hat_ol.begin(), poi_hat_ol.begin() + 2 * firstParameter);
    poi_down_ol.erase(poi_down_ol.begin(), poi_down_ol.begin() + 2 * firstParameter);
    poi_up_ol.erase(poi_up_ol.begin(), poi_up_ol.begin() + 2 * firstParameter);
  }

  poi_nom_down.erase(poi_nom_down.begin(), poi_nom_down.begin() + 2 * firstParameter);
  poi_nom_up.erase(poi_nom_up.begin(), poi_nom_up.begin() + 2 * firstParameter);

  if (overlayCard != "") {
    poi_nom_down_ol.erase(poi_nom_down_ol.begin(), poi_nom_down_ol.begin() + 2 * firstParameter);
    poi_nom_up_ol.erase(poi_nom_up_ol.begin(), poi_nom_up_ol.begin() + 2 * firstParameter);
  }

  boxdown.erase(boxdown.begin(), boxdown.begin() + firstParameter);
  boxup.erase(boxup.begin(), boxup.begin() + firstParameter);
  cendown.erase(cendown.begin(), cendown.begin() + 2 * firstParameter);
  cenup.erase(cenup.begin(), cenup.begin() + 2 * firstParameter);

  if (overlayCard != "") {
    boxdown_ol.erase(boxdown_ol.begin(), boxdown_ol.begin() + firstParameter);
    boxup_ol.erase(boxup_ol.begin(), boxup_ol.begin() + firstParameter);
    cendown_ol.erase(cendown_ol.begin(), cendown_ol.begin() + 2 * firstParameter);
    cenup_ol.erase(cenup_ol.begin(), cenup_ol.begin() + 2 * firstParameter);
  }

  reverse(boxdown_ol.begin(), boxdown_ol.end());
  reverse(boxup_ol.begin(), boxup_ol.end());
  reverse(cendown_ol.begin(), cendown_ol.end());
  reverse(cenup_ol.begin(), cenup_ol.end());
  reverse(down_ol.begin(), down_ol.end());
  reverse(labels_ol.begin(), labels_ol.end());
  reverse(poi_down_ol.begin(), poi_down_ol.end());
  reverse(poi_hat_ol.begin(), poi_hat_ol.end());
  reverse(poi_nom_down_ol.begin(), poi_nom_down_ol.end());
  reverse(poi_nom_up_ol.begin(), poi_nom_up_ol.end());
  reverse(poi_up_ol.begin(), poi_up_ol.end());
  reverse(points_nuis_ol.begin(), points_nuis_ol.end());
  reverse(up_ol.begin(), up_ol.end());
  reverse(val_ol.begin(), val_ol.end());
  reverse(boxdown.begin(), boxdown.end());
  reverse(boxup.begin(), boxup.end());
  reverse(cendown.begin(), cendown.end());
  reverse(cenup.begin(), cenup.end());
  reverse(down.begin(), down.end());
  reverse(labels.begin(), labels.end());
  reverse(poi_down.begin(), poi_down.end());
  reverse(poi_hat.begin(), poi_hat.end());
  reverse(poi_nom_down.begin(), poi_nom_down.end());
  reverse(poi_nom_up.begin(), poi_nom_up.end());
  reverse(poi_up.begin(), poi_up.end());
  reverse(points_nuis.begin(), points_nuis.end());
  reverse(up.begin(), up.end());
  reverse(val.begin(), val.end());

  nrNuis -= firstParameter;

  if (showTopParameters != -1)
    nrNuis2remove = nrNuis - showTopParameters;

  labels.erase(labels.begin(), labels.begin() + nrNuis2remove);
  points_nuis.erase(points_nuis.end() - nrNuis2remove, points_nuis.end());

  if (overlayCard != "") {
    labels_ol.erase(labels_ol.begin(), labels_ol.begin() + nrNuis2remove);
    points_nuis_ol.erase(points_nuis_ol.end() - nrNuis2remove, points_nuis_ol.end());
  }

  val.erase(val.begin(), val.begin() + nrNuis2remove);
  down.erase(down.begin(), down.begin() + nrNuis2remove);
  up.erase(up.begin(), up.begin() + nrNuis2remove);

  if (overlayCard != "") {
    val_ol.erase(val_ol.begin(), val_ol.begin() + nrNuis2remove);
    down_ol.erase(down_ol.begin(), down_ol.begin() + nrNuis2remove);
    up_ol.erase(up_ol.begin(), up_ol.begin() + nrNuis2remove);
  }

  poi_hat.erase(poi_hat.begin(), poi_hat.begin() + 2 * nrNuis2remove);
  poi_down.erase(poi_down.begin(), poi_down.begin() + 2 * nrNuis2remove);
  poi_up.erase(poi_up.begin(), poi_up.begin() + 2 * nrNuis2remove);

  if (overlayCard != "") {
    poi_hat_ol.erase(poi_hat_ol.begin(), poi_hat_ol.begin() + 2 * nrNuis2remove);
    poi_down_ol.erase(poi_down_ol.begin(), poi_down_ol.begin() + 2 * nrNuis2remove);
    poi_up_ol.erase(poi_up_ol.begin(), poi_up_ol.begin() + 2 * nrNuis2remove);
  }

  poi_nom_down.erase(poi_nom_down.begin(), poi_nom_down.begin() + 2 * nrNuis2remove);
  poi_nom_up.erase(poi_nom_up.begin(), poi_nom_up.begin() + 2 * nrNuis2remove);

  if (overlayCard != "") {
    poi_nom_down_ol.erase(poi_nom_down_ol.begin(), poi_nom_down_ol.begin() + 2 * nrNuis2remove);
    poi_nom_up_ol.erase(poi_nom_up_ol.begin(), poi_nom_up_ol.begin() + 2 * nrNuis2remove);
  }

  boxdown.erase(boxdown.begin(), boxdown.begin() + nrNuis2remove);
  boxup.erase(boxup.begin(), boxup.begin() + nrNuis2remove);
  cendown.erase(cendown.begin(), cendown.begin() + 2 * nrNuis2remove);
  cenup.erase(cenup.begin(), cenup.begin() + 2 * nrNuis2remove);

  if (overlayCard != "") {
    boxdown_ol.erase(boxdown_ol.begin(), boxdown_ol.begin() + nrNuis2remove);
    boxup_ol.erase(boxup_ol.begin(), boxup_ol.begin() + nrNuis2remove);
    cendown_ol.erase(cendown_ol.begin(), cendown_ol.begin() + 2 * nrNuis2remove);
    cenup_ol.erase(cenup_ol.begin(), cenup_ol.begin() + 2 * nrNuis2remove);
  }

  nrNuis -= nrNuis2remove;
  nrNuis_ol -= nrNuis2remove;
  cout << "INFO::" << nrNuis << " (" << nrNuis_ol << ") nuisance paramters remaining." << endl;

  int offset = ceil(2 * nrNuis / 10); // used for space to plot the labels and legend

  for (int i = 0; i < nrNuis; i++) {
    poi_up[2 * i] = fabs(poi_up[2 * i]) * scale_poi / max_poi;
    poi_up[2 * i + 1] = fabs(poi_up[2 * i + 1]) * scale_poi / max_poi;
    poi_down[2 * i] = fabs(poi_down[2 * i]) * scale_poi / max_poi;
    poi_down[2 * i + 1] = fabs(poi_down[2 * i + 1]) * scale_poi / max_poi;

    if (overlayCard != "") {
      poi_up_ol[2 * i] = fabs(poi_up_ol[2 * i]) * scale_poi / max_poi;
      poi_up_ol[2 * i + 1] = fabs(poi_up_ol[2 * i + 1]) * scale_poi / max_poi;
      poi_down_ol[2 * i] = fabs(poi_down_ol[2 * i]) * scale_poi / max_poi;
      poi_down_ol[2 * i + 1] = fabs(poi_down_ol[2 * i + 1]) * scale_poi / max_poi;
    }

    poi_nom_up[2 * i] = fabs(poi_nom_up[2 * i]) * scale_poi / max_poi;
    poi_nom_up[2 * i + 1] = fabs(poi_nom_up[2 * i + 1]) * scale_poi / max_poi;
    poi_nom_down[2 * i] = fabs(poi_nom_down[2 * i]) * scale_poi / max_poi;
    poi_nom_down[2 * i + 1] = fabs(poi_nom_down[2 * i + 1]) * scale_poi / max_poi;

    if (overlayCard != "") {
      poi_nom_up_ol[2 * i] = fabs(poi_nom_up_ol[2 * i]) * scale_poi / max_poi;
      poi_nom_up_ol[2 * i + 1] = fabs(poi_nom_up_ol[2 * i + 1]) * scale_poi / max_poi;
      poi_nom_down_ol[2 * i] = fabs(poi_nom_down_ol[2 * i]) * scale_poi / max_poi;
      poi_nom_down_ol[2 * i + 1] = fabs(poi_nom_down_ol[2 * i + 1]) * scale_poi / max_poi;
    }

    if (useRelativeImpact) {
      poi_up[2 * i] /= sigma_tot_hi;
      poi_up[2 * i + 1] /= sigma_tot_hi;
      poi_down[2 * i] /= sigma_tot_lo;
      poi_down[2 * i + 1] /= sigma_tot_lo;

      if (overlayCard != "") {
        poi_up_ol[2 * i] /= sigma_tot_ol_hi;
        poi_up_ol[2 * i + 1] /= sigma_tot_ol_hi;
        poi_down_ol[2 * i] /= sigma_tot_ol_lo;
        poi_down_ol[2 * i + 1] /= sigma_tot_ol_lo;
      }

      poi_nom_up[2 * i] /= sigma_tot_hi;
      poi_nom_up[2 * i + 1] /= sigma_tot_hi;
      poi_nom_down[2 * i] /= sigma_tot_lo;
      poi_nom_down[2 * i + 1] /= sigma_tot_lo;

      if (overlayCard != "") {
        poi_nom_up_ol[2 * i] /= sigma_tot_ol_hi;
        poi_nom_up_ol[2 * i + 1] /= sigma_tot_ol_hi;
        poi_nom_down_ol[2 * i] /= sigma_tot_ol_lo;
        poi_nom_down_ol[2 * i + 1] /= sigma_tot_ol_lo;
      }
    }

    up[i] = fabs(up[i]) * scale_theta;
    down[i] = fabs(down[i]) * scale_theta;

    if (overlayCard != "") {
      up_ol[i] = fabs(up_ol[i]) * scale_theta;
      down_ol[i] = fabs(down_ol[i]) * scale_theta;
    }
  }

  // Finally do the plotting
  pad1->cd();

  // make plot of pulls for nuisance parameters
  int markerSize = 1;
  TGraphAsymmErrors *gr = makeGraphErr("", nrNuis, getAry(val), getAry(points_nuis), getAry(down), getAry(up), NULL, NULL);
  gr->SetLineColor(kBlack);
  gr->SetMarkerColor(kBlack);
  gr->SetMarkerStyle(20);
  gr->SetLineStyle(1);
  gr->SetLineWidth(1);
  gr->SetMarkerSize(markerSize);
  gr->GetXaxis()->SetTitleOffset(1.2);

  TGraphAsymmErrors *gr_ol;
  if (overlayCard != "") {
    gr_ol = makeGraphErr("", nrNuis_ol, getAry(val_ol), getAry(points_nuis_ol), getAry(down_ol), getAry(up_ol), NULL, NULL);
    gr_ol->SetLineColor(kBlack);
    gr_ol->SetMarkerColor(kBlack);
    gr_ol->SetMarkerStyle(24);
    gr_ol->SetLineStyle(1);
    gr_ol->SetLineWidth(2);
    gr_ol->SetMarkerSize(markerSize);
    gr_ol->GetXaxis()->SetTitleOffset(1.2);
  }

  // make plot of 1 sigma boxes
  TGraphAsymmErrors *gr1s = makeGraphErr("", nrNuis, getAry(val), getAry(points_nuis), getAry(boxdown), getAry(boxup), NULL, NULL);
  gr1s->SetLineColor(color_standardband);
  gr1s->SetMarkerColor(color_standardband);
  gr1s->SetLineStyle(1);
  gr1s->SetLineWidth(1);
  gr1s->SetMarkerSize(markerSize * 1.25);
  gr1s->GetXaxis()->SetTitleOffset(1.2);

  TGraphAsymmErrors *gr1s_ol;
  if (overlayCard != "") {
    gr1s_ol = makeGraphErr("", nrNuis_ol, getAry(val_ol), getAry(points_nuis_ol), getAry(boxdown_ol), getAry(boxup_ol), NULL, NULL);
    gr1s_ol->SetLineColor(color_standardband_ol);
    gr1s_ol->SetMarkerColor(color_standardband_ol);
    gr1s_ol->SetLineStyle(1);
    gr1s_ol->SetLineWidth(3);
    gr1s_ol->SetMarkerSize(markerSize * 1.25);
    gr1s_ol->GetXaxis()->SetTitleOffset(1.2);
  }

  vector<double> plot_poi_up_corr;
  vector<double> plot_poi_down_corr;
  vector<double> plot_poi_up_ol_corr;
  vector<double> plot_poi_down_ol_corr;
  vector<double> plot_poi_nom_up_corr;
  vector<double> plot_poi_nom_down_corr;
  vector<double> plot_poi_nom_up_ol_corr;
  vector<double> plot_poi_nom_down_ol_corr;

  vector<double> plot_poi_up_anticorr;
  vector<double> plot_poi_down_anticorr;
  vector<double> plot_poi_up_ol_anticorr;
  vector<double> plot_poi_down_ol_anticorr;
  vector<double> plot_poi_nom_up_anticorr;
  vector<double> plot_poi_nom_down_anticorr;
  vector<double> plot_poi_nom_up_ol_anticorr;
  vector<double> plot_poi_nom_down_ol_anticorr;

  for (int i = 0; i < nrNuis; i++) {
    if (std::find(antiCorrelated.begin(), antiCorrelated.end(), labels[i]) != antiCorrelated.end()) {
      plot_poi_up_corr.push_back(0.0);
      plot_poi_up_corr.push_back(0.0);
      plot_poi_up_anticorr.push_back(poi_up[2 * i]);
      plot_poi_up_anticorr.push_back(poi_up[2 * i + 1]);
      plot_poi_down_corr.push_back(0.0);
      plot_poi_down_corr.push_back(0.0);
      plot_poi_down_anticorr.push_back(poi_down[2 * i]);
      plot_poi_down_anticorr.push_back(poi_down[2 * i + 1]);
    } else {
      plot_poi_up_anticorr.push_back(0.0);
      plot_poi_up_anticorr.push_back(0.0);
      plot_poi_up_corr.push_back(poi_up[2 * i]);
      plot_poi_up_corr.push_back(poi_up[2 * i + 1]);
      plot_poi_down_anticorr.push_back(0.0);
      plot_poi_down_anticorr.push_back(0.0);
      plot_poi_down_corr.push_back(poi_down[2 * i]);
      plot_poi_down_corr.push_back(poi_down[2 * i + 1]);
    }

    if (std::find(antiCorrelatedNom.begin(), antiCorrelatedNom.end(), labels[i]) != antiCorrelatedNom.end()) {
      plot_poi_nom_up_corr.push_back(0.0);
      plot_poi_nom_up_corr.push_back(0.0);
      plot_poi_nom_up_anticorr.push_back(poi_nom_up[2 * i]);
      plot_poi_nom_up_anticorr.push_back(poi_nom_up[2 * i + 1]);
      plot_poi_nom_down_corr.push_back(0.0);
      plot_poi_nom_down_corr.push_back(0.0);
      plot_poi_nom_down_anticorr.push_back(poi_nom_down[2 * i]);
      plot_poi_nom_down_anticorr.push_back(poi_nom_down[2 * i + 1]);
    } else {
      plot_poi_nom_up_anticorr.push_back(0.0);
      plot_poi_nom_up_anticorr.push_back(0.0);
      plot_poi_nom_up_corr.push_back(poi_nom_up[2 * i]);
      plot_poi_nom_up_corr.push_back(poi_nom_up[2 * i + 1]);
      plot_poi_nom_down_anticorr.push_back(0.0);
      plot_poi_nom_down_anticorr.push_back(0.0);
      plot_poi_nom_down_corr.push_back(poi_nom_down[2 * i]);
      plot_poi_nom_down_corr.push_back(poi_nom_down[2 * i + 1]);
    }

    if (overlayCard != "") {
      if (std::find(antiCorrelated_ol.begin(), antiCorrelated_ol.end(), labels_ol[i]) != antiCorrelated_ol.end()) {
        plot_poi_up_ol_corr.push_back(0.0);
        plot_poi_up_ol_corr.push_back(0.0);
        plot_poi_up_ol_anticorr.push_back(poi_up_ol[2 * i]);
        plot_poi_up_ol_anticorr.push_back(poi_up_ol[2 * i + 1]);
        plot_poi_down_ol_corr.push_back(0.0);
        plot_poi_down_ol_corr.push_back(0.0);
        plot_poi_down_ol_anticorr.push_back(poi_down_ol[2 * i]);
        plot_poi_down_ol_anticorr.push_back(poi_down_ol[2 * i + 1]);
      } else {
        plot_poi_up_ol_anticorr.push_back(0.0);
        plot_poi_up_ol_anticorr.push_back(0.0);
        plot_poi_up_ol_corr.push_back(poi_up_ol[2 * i]);
        plot_poi_up_ol_corr.push_back(poi_up_ol[2 * i + 1]);
        plot_poi_down_ol_anticorr.push_back(0.0);
        plot_poi_down_ol_anticorr.push_back(0.0);
        plot_poi_down_ol_corr.push_back(poi_down_ol[2 * i]);
        plot_poi_down_ol_corr.push_back(poi_down_ol[2 * i + 1]);
      }

      if (std::find(antiCorrelatedNom_ol.begin(), antiCorrelatedNom_ol.end(), labels_ol[i]) != antiCorrelatedNom_ol.end()) {
        plot_poi_nom_up_ol_corr.push_back(0.0);
        plot_poi_nom_up_ol_corr.push_back(0.0);
        plot_poi_nom_up_ol_anticorr.push_back(poi_nom_up_ol[2 * i]);
        plot_poi_nom_up_ol_anticorr.push_back(poi_nom_up_ol[2 * i + 1]);
        plot_poi_nom_down_ol_corr.push_back(0.0);
        plot_poi_nom_down_ol_corr.push_back(0.0);
        plot_poi_nom_down_ol_anticorr.push_back(poi_nom_down_ol[2 * i]);
        plot_poi_nom_down_ol_anticorr.push_back(poi_nom_down_ol[2 * i + 1]);
      } else {
        plot_poi_nom_up_ol_anticorr.push_back(0.0);
        plot_poi_nom_up_ol_anticorr.push_back(0.0);
        plot_poi_nom_up_ol_corr.push_back(poi_nom_up_ol[2 * i]);
        plot_poi_nom_up_ol_corr.push_back(poi_nom_up_ol[2 * i + 1]);
        plot_poi_nom_down_ol_anticorr.push_back(0.0);
        plot_poi_nom_down_ol_anticorr.push_back(0.0);
        plot_poi_nom_down_ol_corr.push_back(poi_nom_down_ol[2 * i]);
        plot_poi_nom_down_ol_corr.push_back(poi_nom_down_ol[2 * i + 1]);
      }
    }
  }

  // make plot for the POI change for postfit uncertainties
  TGraphAsymmErrors *gr_poi = makeGraphErr("", 2 * nrNuis, getAry(poi_hat), getAry(points_nuis2), getAry(plot_poi_down_corr), getAry(plot_poi_up_corr), getAry(cenup), getAry(cendown));
  gr_poi->SetLineColor(color_postfit);
  gr_poi->SetFillColor(color_postfit);
  // gr_poi->SetFillStyle(3354);
  gr_poi->SetLineWidth(0);
  gr_poi->SetMarkerSize(0);

  TGraphAsymmErrors *gr_poi_ol;
  if (overlayCard != "") {
    gr_poi_ol = makeGraphErr("", 2 * nrNuis_ol, getAry(poi_hat_ol), getAry(points_nuis2_ol), getAry(plot_poi_down_ol_corr), getAry(plot_poi_up_ol_corr), getAry(cenup_ol), getAry(cendown_ol));
    gr_poi_ol->SetLineColor(color_postfit_ol);
    gr_poi_ol->SetFillColor(color_postfit_ol);
    // gr_poi_ol->SetFillStyle(3345);
    gr_poi_ol->SetLineWidth(0);
    gr_poi_ol->SetMarkerSize(0);
  }

  TGraphAsymmErrors *gr_poi_anticorr = makeGraphErr("", 2 * nrNuis, getAry(poi_hat), getAry(points_nuis2), getAry(plot_poi_down_anticorr), getAry(plot_poi_up_anticorr), getAry(cenup), getAry(cendown));
  gr_poi_anticorr->SetLineColor(visualizeCorrelation ? color_postfit_ol : color_postfit);
  gr_poi_anticorr->SetFillColor(visualizeCorrelation ? color_postfit_ol : color_postfit);
  // gr_poi_anticorr->SetFillStyle(3354);
  gr_poi_anticorr->SetLineWidth(0);
  gr_poi_anticorr->SetMarkerSize(0);

  TGraphAsymmErrors *gr_poi_anticorr_ol;
  if (overlayCard != "") {
    gr_poi_anticorr_ol = makeGraphErr("", 2 * nrNuis_ol, getAry(poi_hat_ol), getAry(points_nuis2_ol), getAry(plot_poi_down_ol_anticorr), getAry(plot_poi_up_ol_anticorr), getAry(cenup_ol), getAry(cendown_ol));
    gr_poi_anticorr_ol->SetLineColor(visualizeCorrelation ? color_postfit : color_postfit_ol);
    gr_poi_anticorr_ol->SetFillColor(visualizeCorrelation ? color_postfit : color_postfit_ol);
    // gr_poi_anticorr_ol->SetFillStyle(3345);
    gr_poi_anticorr_ol->SetLineWidth(0);
    gr_poi_anticorr_ol->SetMarkerSize(0);
  }

  // make plot for the POI change for prefit uncertainties
  TGraphAsymmErrors *gr_poi_nom = makeGraphErr("", 2 * nrNuis, getAry(poi_hat), getAry(points_nuis2), getAry(plot_poi_nom_down_corr), getAry(plot_poi_nom_up_corr), getAry(cenup), getAry(cendown));
  gr_poi_nom->SetLineColor(color_prefit);
  gr_poi_nom->SetFillColor(color_prefit);
  gr_poi_nom->SetFillStyle(3354);
  gr_poi_nom->SetLineWidth(1);
  gr_poi_nom->SetMarkerSize(0);

  TGraphAsymmErrors *gr_poi_nom_ol;
  if (overlayCard != "") {
    gr_poi_nom_ol = makeGraphErr("", 2 * nrNuis_ol, getAry(poi_hat_ol), getAry(points_nuis2_ol), getAry(plot_poi_nom_down_ol_corr), getAry(plot_poi_nom_up_ol_corr), getAry(cenup_ol), getAry(cendown_ol));
    gr_poi_nom_ol->SetLineColor(color_prefit_ol);
    gr_poi_nom_ol->SetFillColor(color_prefit_ol);
    gr_poi_nom_ol->SetFillStyle(3345);
    gr_poi_nom_ol->SetLineWidth(1);
    gr_poi_nom_ol->SetMarkerSize(0);
  }

  TGraphAsymmErrors *gr_poi_anticorr_nom = makeGraphErr("", 2 * nrNuis, getAry(poi_hat), getAry(points_nuis2), getAry(plot_poi_nom_down_anticorr), getAry(plot_poi_nom_up_anticorr), getAry(cenup), getAry(cendown));
  gr_poi_anticorr_nom->SetLineColor(visualizeCorrelation ? color_prefit_ol : color_prefit);
  gr_poi_anticorr_nom->SetFillColor(visualizeCorrelation ? color_prefit_ol : color_prefit);
  gr_poi_anticorr_nom->SetFillStyle(3354);
  gr_poi_anticorr_nom->SetLineWidth(1);
  gr_poi_anticorr_nom->SetMarkerSize(0);

  TGraphAsymmErrors *gr_poi_anticorr_nom_ol;
  if (overlayCard != "") {
    gr_poi_anticorr_nom_ol = makeGraphErr("", 2 * nrNuis_ol, getAry(poi_hat_ol), getAry(points_nuis2_ol), getAry(plot_poi_nom_down_ol_anticorr), getAry(plot_poi_nom_up_ol_anticorr), getAry(cenup_ol), getAry(cendown_ol));
    gr_poi_anticorr_nom_ol->SetLineColor(visualizeCorrelation ? color_prefit : color_prefit_ol);
    gr_poi_anticorr_nom_ol->SetFillColor(visualizeCorrelation ? color_prefit : color_prefit_ol);
    gr_poi_anticorr_nom_ol->SetFillStyle(3345);
    gr_poi_anticorr_nom_ol->SetLineWidth(1);
    gr_poi_anticorr_nom_ol->SetMarkerSize(0);
  }

  double border_lo = -sigma_tot_lo / max_poi;
  double border_hi = sigma_tot_hi / max_poi;

  // different shades for better readability
  int nrShades = ceil((nrNuis + 1) / 2);
  std::vector<double> shadeCenter;
  std::vector<double> shadePoints;
  std::vector<double> shadeWidth;
  std::vector<double> shadeHeight;
  for (int ishade = 0; ishade < nrShades; ishade++) {
    shadeCenter.push_back(0.0);
    shadePoints.push_back(2.0 * ishade + 0.5);
    shadeWidth.push_back(100); // TODO: should not be hardcoded
    shadeHeight.push_back(0.5);
  }

  TGraphAsymmErrors *gr_shades = makeGraphErr("", nrShades, getAry(shadeCenter), getAry(shadePoints), getAry(shadeWidth), getAry(shadeWidth), getAry(shadeHeight), getAry(shadeHeight));
  gr_shades->SetLineColor(18);
  gr_shades->SetFillColor(18);
  gr_shades->SetFillStyle(3001);
  gr_shades->SetLineWidth(1);
  gr_shades->SetMarkerSize(0);

  // histogram to get the nuisance parameter labels correct
  TH2F *h = new TH2F("h", "", 1, border_lo, border_hi, nrNuis + offset + 1, -offset, nrNuis + 1);
  for (int i = offset; i < nrNuis + offset; i++)
    h->GetYaxis()->SetBinLabel(i + 1, drawParamNames ? labels[i - offset].c_str() : "");
  // for (int i = offset; i < nrNuis + offset; i++) h->GetYaxis()->SetBinLabel(i + 1, "");
  h->LabelsOption("h");
  double labelSize = 1. / nrNuis;
  h->SetLabelSize(labelSize > 0.02 ? 0.02 : labelSize, "Y");
  h->GetXaxis()->SetLabelColor(kWhite);
  h->GetXaxis()->SetAxisColor(kWhite);
  h->GetYaxis()->SetLabelColor(kBlack);
  h->GetYaxis()->SetAxisColor(kBlack);
  h->GetYaxis()->SetTickLength(0.);
  h->SetStats(0);
  // h->LabelsDeflate();
  h->Draw("h");

  // TODO coloring depending on tags

  // axis for the POI correlation
  TGaxis *axis_poi = new TGaxis(border_lo, nrNuis + 1, border_hi, nrNuis + 1, (-sigma_tot_lo) / scale_poi, (sigma_tot_hi) / scale_poi, 510, "-L");
  axis_poi->ImportAxisAttributes(h->GetXaxis());
  axis_poi->SetName("axis_poi");
  stringstream ss_poi;
  ss_poi << "#Delta#hat{" << poiname << "}";
  if (useRelativeImpact)
    axis_poi->SetTitle("#Delta#hat{#mu}/#Delta#hat{#mu}_{tot}");
  // else axis_poi->SetTitle("#Delta#hat{m}_{H} [GeV]");
  else
    axis_poi->SetTitle(ss_poi.str().c_str());
  axis_poi->SetTitleOffset(1.1);
  axis_poi->SetLineColor(kBlack);
  axis_poi->SetLabelColor(kBlack);
  axis_poi->SetTitleColor(kBlack);
  axis_poi->SetLabelSize(0.034);
  axis_poi->SetTitleSize(0.034);

  // axis for the nuisance parameter pull
  TGaxis *axis_theta = new TGaxis(border_lo, -offset, border_hi, -offset, (-sigma_tot_lo / max_poi) / scale_theta, (sigma_tot_hi / max_poi) / scale_theta, 510, "+R");
  axis_theta->ImportAxisAttributes(h->GetXaxis());
  axis_theta->SetName("axis_theta");
  axis_theta->SetTitle("(#hat{#theta} - #theta_{0})/#Delta#theta");
  axis_theta->SetTitleOffset(1.1);
  axis_theta->SetLineColor(kBlack);
  axis_theta->SetLabelColor(kBlack);
  axis_theta->SetTitleColor(kBlack);
  axis_theta->SetLabelSize(0.034);
  axis_theta->SetTitleSize(0.034);

  // axis for the nuisance parameter labels
  TGaxis *axis_label = new TGaxis(border_lo, 0, border_lo, nrNuis + 1, 0, nrNuis + 1, 0, "-R");
  axis_label->SetLineColor(kBlack);
  axis_label->SetTitleColor(kWhite);
  axis_label->SetLabelSize(0);
  axis_label->SetNdivisions(0);

  // some line definitions
  TLine l;
  l.SetLineWidth(2);
  l.SetLineColor(color_pulls);
  l.SetLineStyle(2);

  // draw the nuisance parameter pulls including error bands and impact on poi
  gr_shades->Draw("p2");

  if (drawPostfitImpactBand) {
    gr_poi->Draw("p2");
    gr_poi_anticorr->Draw("p2");
    if (overlayCard != "") {
      gr_poi_ol->Draw("p2");
      gr_poi_anticorr_ol->Draw("p2");
    }
  }
  if (drawPrefitImpactBand) {
    gr_poi_nom->Draw("p2");
    gr_poi_anticorr_nom->Draw("p2");
    if (overlayCard != "") {
      gr_poi_nom_ol->Draw("p2");
      gr_poi_anticorr_nom_ol->Draw("p2");
    }
  }

  // draw axes
  if (drawPrefitImpactBand || drawPostfitImpactBand)
    axis_poi->Draw();
  axis_theta->Draw();
  axis_label->Draw();

  // draw +-1 and 0 sigma lines for pulls
  l.DrawLine(0., 0., 0., nrNuis);
  l.DrawLine(1. * scale_theta, 0., 1. * scale_theta, nrNuis);
  l.DrawLine(-1. * scale_theta, 0., -1. * scale_theta, nrNuis);

  gStyle->SetEndErrorSize(5.0);
  gr1s->Draw("p");
  if (overlayCard != "")
    gr1s_ol->Draw("p");

  gr->Draw("p");
  if (overlayCard != "")
    gr_ol->Draw("p");

  pad1->SetTicks(0, 0);

  c1->SetLogy(0);

  double channelPosX = 0.36;
  double channelPosY = 0.19;
  TLegend *leg = new TLegend(channelPosX + 0.28, channelPosY - 0.0775, channelPosX + 0.68, channelPosY + 0.02, "", "NDC");
  leg->SetFillStyle(0);
  leg->SetTextSize(0.0225);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  leg->AddEntry(gr, "Pull", "lp");
  if (nrNuis_ol > 0)
    leg->AddEntry(gr_ol, "Alt pull", "lp");
  leg->AddEntry(gr1s, "1 standard deviation", "l");
  if (nrNuis_ol > 0)
    leg->AddEntry(gr1s_ol, "Alt 1 standard deviation", "l");
  if (drawPrefitImpactBand) {
    stringstream ss_poi;
    ss_poi << "Prefit Impact on #hat{" << poiname << "}";
    leg->AddEntry(gr_poi_nom, ss_poi.str().c_str(), "f");
    if (nrNuis_ol > 0)
      leg->AddEntry(gr_poi_nom_ol, "Alt Prefit Impact on #hat{#mu}", "f");
  }
  if (drawPostfitImpactBand) {
    stringstream ss_poi;
    ss_poi << "Postfit Impact on #hat{" << poiname << "}";
    leg->AddEntry(gr_poi, ss_poi.str().c_str(), "f");
    if (nrNuis_ol > 0)
      leg->AddEntry(gr_poi_ol, "Alt Postfit Impact on #hat{#mu}", "f");
  }

  leg->Draw();

  // Labels
  double position_label_x = 0.36;
  double position_label_y = 0.185;
  double position_label_delta_y = 0.0275;
  double position_label_delta_x1 = 0.1225;
  double position_label_delta_x2 = 0.185;

  TLatex label_tex_1;
  label_tex_1.SetNDC();
  label_tex_1.SetTextFont(72);
  label_tex_1.SetTextColor(kBlack);
  label_tex_1.SetTextSize(0.70 * label_tex_1.GetTextSize());

  TLatex label_tex_2;
  label_tex_2.SetNDC();
  label_tex_2.SetTextFont(42);
  label_tex_2.SetTextColor(1);
  label_tex_2.SetTextSize(0.70 * label_tex_2.GetTextSize());

  label_tex_1.DrawLatex(position_label_x, position_label_y, "ATLAS");
  label_tex_2.DrawLatex(position_label_x + position_label_delta_x1, position_label_y, "Internal");

  stringstream lumiLatex;
  lumiLatex << "#sqrt{s} = " << 13 << " TeV, " << 36.5 << " fb^{-1}";

  stringstream ranklabel;
  ranklabel << "Rank " << firstParameter + 1 << " to " << (firstParameter + showTopParameters == -1 ? labels.size() : firstParameter + showTopParameters);

  label_tex_2.DrawLatex(position_label_x, position_label_y - 1 * position_label_delta_y, "LHC Run 1");
  label_tex_2.DrawLatex(position_label_x, position_label_y - 2 * position_label_delta_y, ranklabel.str().c_str());

  // Cleanup
  system(("rm -f " + cardName + "/tmp.root").c_str());

  // Save the plot
  stringstream baseName;
  // baseName << cardName << "_" << poiname << "_rank_" << setfill('0') << setw(4) << firstParameter + 1 << "_to_" << setfill('0') << setw(4) << firstParameter + showTopParameters;
  baseName << "ranking_" << poiname << "_rank_" << setfill('0') << setw(4) << firstParameter + 1 << "_to_" << setfill('0') << setw(4) << firstParameter + showTopParameters;

  string type = "pdf";
  system(("mkdir -vp " + type + "-files").c_str());

  stringstream saveName;
  saveName << type << "-files/" << baseName.str() << "." << type;

  c1->SaveAs(saveName.str().c_str());

  // Finish plotting
  PrintResourcesUsed(thistime);
}
