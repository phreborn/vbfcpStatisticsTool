// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2018-04-10
// Description : Draw profile likelihood scans

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
#include <numeric>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TMultiGraph.h"
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
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooGaussian.h"
#include "RooHistFunc.h"
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

#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasStyle.h"

#include "boost/program_options.hpp"
#include "boost/program_options/cmdline.hpp"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"

using namespace std;
using namespace RooFit;
using namespace RooStats;

// _____________________________________________________________________________
// Declarations of functions used in this file
void save(string baseName, string type, TCanvas* c1);
void readFiles(vector<string> filenames, vector<string> poinames, vector<string> nuis,
               vector<pair<string, map<string, map< vector<double>, double > > > > &map_overlay2folder2poi2nll,
               vector<pair<string, map<string, map< vector<double>, int > > > > &map_overlay2folder2poi2status,
               vector<pair<string, map<string, map< vector<double>, map<string, double> > > > > &map_overlay2folder2poi2nuis);
void plot1D(vector<string> filenames, string poi, vector<string> nuis = vector<string>(),
            vector<pair<string, map<string, map< vector<double>, double > > > > map_overlay2folder2poi2nll = vector<pair<string, map<string, map< vector<double>, double > > > >(),
            vector<pair<string, map<string, map< vector<double>, int > > > > map_overlay2folder2poi2status = vector<pair<string, map<string, map< vector<double>, int > > > >(),
            vector<pair<string, map<string, map< vector<double>, map<string, double> > > > > map_overlay2folder2poi2nuis = vector<pair<string, map<string, map< vector<double>, map<string, double> > > > >(),
            double x_lo = -3.0, double x_hi = 3.0, double y_lo = 0.0, double y_hi = 10.0,
            vector<string> color = vector<string>(), vector<int> style = vector<int>(), vector<string> legend = vector<string>(), vector<string> axis_label = vector<string>(),
            string label = "", string luminosity = "", bool twoStepInterpolation=false);
void plot2D(vector<string> filenames, vector<string> poi, vector<string> nuis = vector<string>(),
            vector<pair<string, map<string, map< vector<double>, double > > > > map_overlay2folder2poi2nll = vector<pair<string, map<string, map< vector<double>, double > > > >(),
            vector<pair<string, map<string, map< vector<double>, int > > > > map_overlay2folder2poi2status = vector<pair<string, map<string, map< vector<double>, int > > > >(),
            vector<pair<string, map<string, map< vector<double>, map<string, double> > > > > map_overlay2folder2poi2nuis = vector<pair<string, map<string, map< vector<double>, map<string, double> > > > >(),
            double x_lo = -3.0, double x_hi = 3.0, double y_lo = -3.0, double y_hi = 3.0,
            vector<string> color = vector<string>(), vector<int> style = vector<int>(), vector<string> legend = vector<string>(), vector<string> axis_label = vector<string>(),
            string label = "", string luminosity = "");
TH2D* IncreaseResolutionAndSmooth(TH2D* h, int xbins, double xlo, double xhi, int ybins, double ylo, double yhi);
map<int, deque<TGraph*> > GetContours(TH2D* h);
TGraph* GetBestFit(TH2D* h);

// _____________________________________________________________________________
// Main routine
int main(int argc, char **argv) {
  TTime thistime = gSystem->Now();
  SetAtlasStyle();

  vector<string> input      = {};
  vector<string> color      = {};
  vector<int> style         = {};
  vector<string> legend     = {};
  vector<string> poi        = {};
  vector<string> axis_label = {};
  vector<string> nuis       = {};
  vector<double> x_range    = {};
  vector<double> y_range    = {};
  string label              = "Internal";
  string luminosity         = "#sqrt{s} = 13 TeV, 36.1 fb^{-1}";
  string loglevel           = "INFO";

  using namespace boost;
  namespace po = boost::program_options;
  po::options_description desc( "Program options" );
  desc.add_options()
    ( "help        , h"                                                            , "Print this help message")
    ( "input"      , po::value<vector<string>> ( &input )->multitoken()            , "Path to input." )
    ( "color"      , po::value<vector<string>> ( &color )->multitoken()            , "Line color." )
    ( "style"      , po::value<vector<int>> ( &style )->multitoken()               , "Line style." )
    ( "legend"     , po::value<vector<string>> ( &legend )->multitoken()           , "Legend text." )
    ( "poi"        , po::value<vector<string>> ( &poi )->multitoken()              , "Parameter of interest." )
    ( "axis_label" , po::value<vector<string>> ( &axis_label )->multitoken()       , "Axis label." )
    ( "nuis"       , po::value<vector<string>> ( &nuis )->multitoken()             , "Nuisance parameters." )
    ( "x"          , po::value<vector<double>> ( &x_range )->multitoken()          , "x-axis range." )
    ( "y"          , po::value<vector<double>> ( &y_range )->multitoken()          , "y-axis range." )
    ( "label"      , po::value<string>( &label )->default_value( label )           , "Internal, Preliminary, etc." )
    ( "luminosity" , po::value<string>( &luminosity )->default_value( luminosity ) , "Luminosity." )
    ( "loglevel"   , po::value<string>( &loglevel )->default_value( loglevel )     , "Loglevel." )
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

  // Bookkeeping
  vector<pair<string, map<string, map< vector<double>, double > > > > map_overlay2folder2poi2nll;
  vector<pair<string, map<string, map< vector<double>, int > > > > map_overlay2folder2poi2status;
  vector<pair<string, map<string, map< vector<double>, map<string, double> > > > > map_overlay2folder2poi2nuis;

  // Read input
  readFiles(input, poi, nuis, map_overlay2folder2poi2nll, map_overlay2folder2poi2status, map_overlay2folder2poi2nuis);

  // Plotting
  if (poi.size() == 1) {
    plot1D(input, poi[0], nuis, map_overlay2folder2poi2nll, map_overlay2folder2poi2status, map_overlay2folder2poi2nuis, x_range[0], x_range[1], y_range[0], y_range[1], color, style, legend, axis_label, label, luminosity);
  } else if (poi.size() == 2) {
    plot2D(input, poi, nuis, map_overlay2folder2poi2nll, map_overlay2folder2poi2status, map_overlay2folder2poi2nuis, x_range[0], x_range[1], y_range[0], y_range[1], color, style, legend, axis_label, label, luminosity);
  } else {
    LOG(logERROR) << "Plotting for multiple pois not yet implemented.";
  }

  // Finish plotting
  PrintResourcesUsed(thistime);
}

// _____________________________________________________________________________
// Save canvas
void save(string baseName, string type, TCanvas* c1) {
  system(("mkdir -vp " + type + "-files").c_str());
  stringstream saveName;
  saveName << type << "-files/" << baseName << "." << type;
  c1->SaveAs(saveName.str().c_str());
}

// _____________________________________________________________________________
// Reads ROOT files with NLL value per scan point
void readFiles(vector<string> filenames, vector<string> poinames, vector<string> nuis,
               vector<pair<string, map<string, map< vector<double>, double > > > > &map_overlay2folder2poi2nll,
               vector<pair<string, map<string, map< vector<double>, int > > > > &map_overlay2folder2poi2status,
               vector<pair<string, map<string, map< vector<double>, map<string, double> > > > > &map_overlay2folder2poi2nuis) {
  LOG(logINFO) << "Reading files.";

  for (auto filename : filenames) {
    vector<string> folders;
    boost::split(folders, filename, boost::is_any_of(";"));

    map<string, map< vector<double>, double > > map_folder2poi2nll;
    map<string, map< vector<double>, int > > map_folder2poi2status;
    map<string, map< vector<double>, map<string, double> > > map_folder2poi2nuis;

    for (auto folder : folders) {
      vector<string> files;
      boost::split(files, folder, boost::is_any_of(","));

      map< vector<double>, double > map_poi2nll;
      map< vector<double>, int > map_poi2status;
      map< vector<double>, map<string, double> > map_poi2nuis;

      for (auto file : files) {
        boost::filesystem::path targetDir(file);
        boost::filesystem::directory_iterator it(targetDir), eod;

        BOOST_FOREACH(boost::filesystem::path const &fname, std::make_pair(it, eod))
        {
          if(boost::filesystem::is_regular_file(fname))
          {
            LOG(logDEBUG) << "On file " << fname;

            TFile* f = TFile::Open(fname.c_str());
            if (!f) {
              LOG(logERROR) << "Could not open file " << fname;
              exit(-1);
            }

            TTree* resultTree = NULL;
            f->GetObject("result", resultTree);
            if (!resultTree) exit(-1);

            int status;
            double nll;
            double edm;

            TBranch *b_status = NULL;
            TBranch *b_nll = NULL;
            TBranch *b_edm = NULL;

            resultTree->SetBranchAddress("status", &status, &b_status);
            resultTree->SetBranchAddress("nll", &nll, &b_nll);
            resultTree->SetBranchAddress("edm", &edm, &b_edm);

            map<string, TBranch*> b_pois;
            map<string, double> pois_val;
            map<string, TBranch*> b_nuis;
            map<string, double> nuis_val;

            for (auto n : nuis) {
              b_nuis[n] = NULL;
              nuis_val[n] = 0.0;
              resultTree->SetBranchAddress(n.c_str(), &nuis_val[n], &b_nuis[n]);
            }

            for (auto p : poinames) {
              b_pois[p] = NULL;
              pois_val[p] = 0.0;
              resultTree->SetBranchAddress(p.c_str(), &pois_val[p], &b_pois[p]);
            }

            Long64_t tentry = resultTree->LoadTree(0);

            b_status->GetEntry(tentry);
            b_nll->GetEntry(tentry);
            b_edm->GetEntry(tentry);

            for (auto n : nuis) {
              b_nuis[n]->GetEntry(tentry);
            }

            for (auto p : poinames) {
              b_pois[p]->GetEntry(tentry);
            }

            if (nll == nll && status == 0) {
              double two_nll = 2.0 * nll;
              vector<double> coordinates;

              for (auto p : poinames) {
                coordinates.push_back(pois_val[p]);
              }

              map_poi2nll[coordinates] = two_nll;
              map_poi2status[coordinates] = status;
              map_poi2nuis[coordinates] = nuis_val;
            }

            f->Close();
            delete f;
          }
        }
      }

      map_folder2poi2nll[folder] = map_poi2nll;
      map_folder2poi2status[folder] = map_poi2status;
      map_folder2poi2nuis[folder] = map_poi2nuis;
    }

    map_overlay2folder2poi2nll.push_back(make_pair(filename, map_folder2poi2nll));
    map_overlay2folder2poi2status.push_back(make_pair(filename, map_folder2poi2status));
    map_overlay2folder2poi2nuis.push_back(make_pair(filename, map_folder2poi2nuis));
  }
}

// ____________________________________________________________________________|__________
// Plot 1D profile likelihood scan
void plot1D(vector<string> filenames, string poi, vector<string> nuis,
            vector<pair<string, map<string, map< vector<double>, double > > > > map_overlay2folder2poi2nll,
            vector<pair<string, map<string, map< vector<double>, int > > > > map_overlay2folder2poi2status,
            vector<pair<string, map<string, map< vector<double>, map<string, double> > > > > map_overlay2folder2poi2nuis,
            double x_lo, double x_hi, double y_lo, double y_hi,
            vector<string> color, vector<int> style, vector<string> legend, vector<string> axis_label,
            string label, string luminosity, bool twoStepInterpolation)
{
  LOG(logINFO) << "Plotting 1D likelihood scan.";

  Color_t color_1sigma  = kRed+1;
  Color_t color_2sigma  = kGreen+2;
  Color_t color_3sigma  = kCyan-6;
  Color_t color_4sigma  = kPink+4;
  Color_t color_5sigma  = kAzure-3;
  Color_t color_bestfit = kOrange+7;

  map<string, map<string, TGraph* > > map_overlay2folder2graph;
  map<string, map<string, TGraph* > > map_overlay2folder2interpolatedgraph;
  map<string, map<string, map<string, TGraph* > > > map_overlay2folder2nuis2graph;
  map<string, map<string, map< string, TGraph* > > > map_overlay2folder2nuis2interpolatedgraph;

  for (auto it_overlay : map_overlay2folder2poi2nll) {
    int i_overlay = &it_overlay - &map_overlay2folder2poi2nll[0];
    string thisOverlay = it_overlay.first;
    map<string, map< vector<double>, double > > map_folder2poi2nll = it_overlay.second;

    map<string, TGraph* > map_folder2graph;
    map<string, TGraph* > map_folder2interpolatedgraph;
    map<string, map<string, TGraph* > > map_folder2nuis2graph;
    map<string, map<string, TGraph* > > map_folder2nuis2interpolatedgraph;

    double xlo =  numeric_limits<double>::infinity();
    double xhi = -numeric_limits<double>::infinity();

    // Order coordinates for every folder and make graph, find poi range
    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      map< vector<double>, double > map_poi2nll = it_folder->second;

      vector<double> x, y;

      for(map< vector<double>, double>::iterator it_poi = map_poi2nll.begin(); it_poi != map_poi2nll.end(); ++it_poi) {
        double nll = it_poi->second;

        if (nll == nll) {
          x.push_back(it_poi->first[0]);
          y.push_back(it_poi->second);
        }
      }

      int nrPoints = x.size();

      for (int i = 0; i < nrPoints-1; i++) {
        for (int j = 0; j < nrPoints-1-i; j++) {
          if (x[j] > x[j+1]) {
            swap(x[j], x[j+1]);
            swap(y[j], y[j+1]);
          }
        }
      }

      if (x[0] < xlo) xlo = x[0];
      if (x[nrPoints-1] > xhi) xhi = x[nrPoints-1];

      TGraph* g = new TGraph(nrPoints, getAry(x), getAry(y));
      map_folder2graph[thisFolder] = g;

      // Make plots of the NPs as function of the POI
      if (nuis.size() > 0) {
        map< vector<double>, map<string, double> > map_poi2nuis = ((map_overlay2folder2poi2nuis[i_overlay]).second)[thisFolder];
        for (size_t inuis = 0; inuis < nuis.size(); inuis++) {
          string thisNuis = nuis[inuis];
          vector<double> xnuis, ynuis;

          for(map< vector<double>, map<string, double> >::iterator it_poi2nuis = map_poi2nuis.begin(); it_poi2nuis != map_poi2nuis.end(); ++it_poi2nuis) {
            double nuisval = it_poi2nuis->second[thisNuis];

            if (nuisval == nuisval) {
              xnuis.push_back(it_poi2nuis->first[0]);
              ynuis.push_back(nuisval);
            }
          }

          int nrPointsNuis = xnuis.size();

          for (int iinuis = 0; iinuis < nrPointsNuis-1; iinuis++) {
            for (int jjnuis = 0; jjnuis < nrPointsNuis-1-iinuis; jjnuis++) {
              if (xnuis[jjnuis] > xnuis[jjnuis+1]) {
                swap(x[jjnuis], x[jjnuis+1]);
                swap(y[jjnuis], y[jjnuis+1]);
              }
            }
          }

          TGraph* gnuis = new TGraph(nrPointsNuis, getAry(xnuis), getAry(ynuis));
          map_folder2nuis2graph[thisFolder][thisNuis] = gnuis;
        }
      }
    }

    double minNll = TMath::Infinity();
    double minNll_x = 0.0;

    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TGraph* g = map_folder2graph[thisFolder];

      for (int i_point = 0; i_point < g->GetN(); ++i_point) {
        double x, y = 0;
        g->GetPoint(i_point, x, y);
        if (y < minNll) {
          minNll = y;
          minNll_x = x;
        }
      }
    }

    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TGraph* g = map_folder2graph[thisFolder];

      for (int i_point = 0; i_point < g->GetN(); ++i_point) {
        double x, y = 0;
        g->GetPoint(i_point, x, y);
        y -= minNll;
        g->SetPoint(i_point, x, y);
      }

      map_folder2graph[thisFolder] = g;
    }

    minNll = TMath::Infinity();
    minNll_x = 0.0;

    // Make smooth interpolated graph for every folder in poi range, find minimum nll
    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TGraph* g = map_folder2graph[thisFolder];

      vector<double> x_interpolated_coarse, y_interpolated_coarse;

      double stepsize_coarse = fabs(xhi - xlo) / 100.0;
      for (double thisX = xlo; thisX <= xhi; thisX += stepsize_coarse) {
        double thisY = g->Eval(thisX, 0);
        x_interpolated_coarse.push_back(thisX);
        y_interpolated_coarse.push_back(thisY);
      }

      int nrPoints_interpolated_coarse = x_interpolated_coarse.size();
      TGraph* g_interpolated_coarse = new TGraph(nrPoints_interpolated_coarse, getAry(x_interpolated_coarse), getAry(y_interpolated_coarse));

      vector<double> x_interpolated, y_interpolated;

      double stepsize = fabs(xhi - xlo) / 500.0;
      for (double thisX = xlo; thisX <= xhi; thisX += stepsize) {
        double thisY = 0.0;
        if (twoStepInterpolation) thisY = g_interpolated_coarse->Eval(thisX, 0, "S");
        else thisY = g->Eval(thisX, 0, "S");
        x_interpolated.push_back(thisX);
        y_interpolated.push_back(thisY);
      }

      int nrPoints_interpolated = x_interpolated.size();
      TGraph* g_interpolated = new TGraph(nrPoints_interpolated, getAry(x_interpolated), getAry(y_interpolated));

      for (int i_point = 0; i_point < g_interpolated->GetN(); ++i_point) {
        double x, y = 0;
        g_interpolated->GetPoint(i_point, x, y);
        if (y < minNll) {
          minNll = y;
          minNll_x = x;
        }
      }

      map_folder2interpolatedgraph[thisFolder] = g_interpolated;

      for (size_t inuis = 0; inuis < nuis.size(); inuis++) {
        string thisNuis = nuis[inuis];
        TGraph* gnuis = map_folder2nuis2graph[thisFolder][thisNuis];

        vector<double> xnuis_interpolated_coarse, ynuis_interpolated_coarse;

        for (double thisX = xlo; thisX <= xhi; thisX += stepsize_coarse) {
          double thisY = gnuis->Eval(thisX, 0);
          xnuis_interpolated_coarse.push_back(thisX);
          ynuis_interpolated_coarse.push_back(thisY);
        }

        int nrPointsNuis_interpolated_coarse = xnuis_interpolated_coarse.size();
        TGraph* gnuis_interpolated_coarse = new TGraph(nrPointsNuis_interpolated_coarse, getAry(xnuis_interpolated_coarse), getAry(ynuis_interpolated_coarse));

        vector<double> xnuis_interpolated, ynuis_interpolated;

        for (double thisX = xlo; thisX <= xhi; thisX += stepsize) {
          double thisY = 0.0;
          if (twoStepInterpolation) thisY = gnuis_interpolated_coarse->Eval(thisX, 0, "S");
          else thisY = gnuis->Eval(thisX, 0, "S");
          xnuis_interpolated.push_back(thisX);
          ynuis_interpolated.push_back(thisY);
        }

        int nrPointsNuis_interpolated = xnuis_interpolated.size();
        TGraph* gnuis_interpolated = new TGraph(nrPointsNuis_interpolated, getAry(xnuis_interpolated), getAry(ynuis_interpolated));

        map_folder2nuis2interpolatedgraph[thisFolder][thisNuis] = gnuis_interpolated;
      }
    }

    // Subtract minimum nll in every point of coarse and smooth graph
    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TGraph* g = map_folder2graph[thisFolder];
      TGraph* g_interpolated = map_folder2interpolatedgraph[thisFolder];

      for (int i_point = 0; i_point < g->GetN(); ++i_point) {
        double x, y = 0;
        g->GetPoint(i_point, x, y);
        y -= minNll;
        g->SetPoint(i_point, x, y);
      }

      for (int i_point = 0; i_point < g_interpolated->GetN(); ++i_point) {
        double x, y = 0;
        g_interpolated->GetPoint(i_point, x, y);
        y -= minNll;
        g_interpolated->SetPoint(i_point, x, y);
      }

      map_folder2graph[thisFolder] = g;
      map_folder2interpolatedgraph[thisFolder] = g_interpolated;
    }

    map_overlay2folder2graph[thisOverlay] = map_folder2graph;
    map_overlay2folder2interpolatedgraph[thisOverlay] = map_folder2interpolatedgraph;
    map_overlay2folder2nuis2graph[thisOverlay] = map_folder2nuis2graph;
    map_overlay2folder2nuis2interpolatedgraph[thisOverlay] = map_folder2nuis2interpolatedgraph;
  }

  map<string, TGraph* > map_overlay2bestinterpolatedgraph;
  map<string, map<string, TGraph* > > map_overlay2nuis2bestinterpolatedgraph;

  for(map<string, map<string, TGraph* > >::iterator it_overlay = map_overlay2folder2interpolatedgraph.begin(); it_overlay != map_overlay2folder2interpolatedgraph.end(); ++it_overlay) {
    string thisOverlay = it_overlay->first;
    map<string, TGraph* > map_folder2interpolatedgraph = it_overlay->second;
    map<string, map<string, TGraph* > > map_folder2nuis2interpolatedgraph = map_overlay2folder2nuis2interpolatedgraph[thisOverlay];

    vector<double> xmin, ymin;
    map<string, vector<double> > xminnuis, yminnuis;
    bool isFirst = true;

    for (map<string, TGraph* >::iterator it_folder = map_folder2interpolatedgraph.begin(); it_folder != map_folder2interpolatedgraph.end(); ++it_folder) {
      TGraph* g_interpolated = it_folder->second;

      for (int i_point = 0; i_point < g_interpolated->GetN(); ++i_point) {
        double x, y = 0;
        g_interpolated->GetPoint(i_point, x, y);
        if (isFirst) {
          xmin.push_back(x);
          ymin.push_back(y);
        } else {
          if (y < ymin[i_point]) ymin[i_point] = y;
        }

        for (size_t inuis = 0; inuis < nuis.size(); inuis++) {
          string thisNuis = nuis[inuis];
          TGraph* gnuis_interpolated = map_folder2nuis2interpolatedgraph[it_folder->first][thisNuis];
          double xnuis, ynuis = 0;
          gnuis_interpolated->GetPoint(i_point, xnuis, ynuis);
          if (isFirst) {
            xminnuis[thisNuis].push_back(xnuis);
            yminnuis[thisNuis].push_back(ynuis);
          } else {
            if (y < ymin[i_point]) yminnuis[thisNuis][i_point] = ynuis;
          }
        }
      }
      isFirst = false;
    }

    int nrPoints_min = xmin.size();
    TGraph* g_min = new TGraph(nrPoints_min, getAry(xmin), getAry(ymin));

    map_overlay2bestinterpolatedgraph[thisOverlay] = g_min;

    for (size_t inuis = 0; inuis < nuis.size(); inuis++) {
      string thisNuis = nuis[inuis];

      int nrPointsNuis_min = xmin.size();
      TGraph* gnuis_min = new TGraph(nrPointsNuis_min, getAry(xminnuis[thisNuis]), getAry(yminnuis[thisNuis]));

      map_overlay2nuis2bestinterpolatedgraph[thisOverlay][thisNuis] = gnuis_min;
    }
  }

  // Final 1D plotting
  TCanvas* c1 = new TCanvas("c1", "c1", 1024, 768);
  TPad *pad1 = new TPad("pad1", "pad1", 0.0 , 0.0 , 1.0 , 1.0 , 0);
  pad1->Draw();
  pad1->cd();

  TMultiGraph* g_all = new TMultiGraph();
  TLegend* leg = new TLegend(0.0, 0.0, 0.0, 0.0, "", "NDC");
  TLegend* legnuis = new TLegend(0.0, 0.0, 0.0, 0.0, "", "NDC");
  map<string, TGraph* > map_overlay2bestfitgraph;
  map<string, TGraph* > map_overlay21sigmagraph;
  map<string, TGraph* > map_overlay22sigmagraph;
  map<string, TGraph* > map_overlay23sigmagraph;
  map<string, TGraph* > map_overlay24sigmagraph;
  map<string, TGraph* > map_overlay25sigmagraph;
  double intervalBoundary_lo = TMath::Infinity();
  double intervalBoundary_hi = -TMath::Infinity();
  double leftmost_minimum = TMath::Infinity();
  double rightmost_minimum = -TMath::Infinity();

  int i_overlay = 0;
  for (auto it_overlay : map_overlay2folder2poi2nll) {
    TString thisOverlay = it_overlay.first;

  // for(map<string, map<string, TGraph* > >::iterator it_overlay = map_overlay2folder2interpolatedgraph.begin(); it_overlay != map_overlay2folder2interpolatedgraph.end(); ++it_overlay) {

    map<string, TGraph* > map_folder2interpolatedgraph = map_overlay2folder2interpolatedgraph[thisOverlay.Data()];
    TGraph* g_min = map_overlay2bestinterpolatedgraph[thisOverlay.Data()];

    Color_t mycolor = TColor::GetColor(color[i_overlay].c_str());
    g_min->SetLineColor(mycolor);
    g_min->SetLineStyle(style[i_overlay]);
    g_min->SetLineWidth(3);

    for (map<string, TGraph* >::iterator it_folder = map_folder2interpolatedgraph.begin(); it_folder != map_folder2interpolatedgraph.end(); ++it_folder) {
      string thisFolder = it_folder->first;

      TGraph* g = map_overlay2folder2graph[thisOverlay.Data()][thisFolder];
      TGraph* g_interpolated = map_overlay2folder2interpolatedgraph[thisOverlay.Data()][thisFolder];

      g->SetMarkerSize(0.8);
      g->SetMarkerStyle(20);

      g->SetMarkerColor(mycolor);
      g_interpolated->SetLineColor(mycolor);
      g_interpolated->SetLineStyle(style[i_overlay]);
      g_interpolated->SetLineWidth(1);
    }

    g_all->Add(g_min, "AL");
    string legendText = legend[i_overlay];

    leg->AddEntry(g_min, legendText.c_str(), "L");

    for (size_t inuis = 0; inuis < nuis.size(); inuis++) {
      string thisNuis = nuis[inuis];
      TGraph* gnuis_min = map_overlay2nuis2bestinterpolatedgraph[thisOverlay.Data()][thisNuis];
      gnuis_min->SetLineColor(2+inuis);
      g_all->Add(gnuis_min, "L");

      for (map<string, TGraph* >::iterator it_folder = map_folder2interpolatedgraph.begin(); it_folder != map_folder2interpolatedgraph.end(); ++it_folder) {
        string thisFolder = it_folder->first;
        TGraph* g = map_overlay2folder2nuis2graph[thisOverlay.Data()][thisFolder][thisNuis];
        g->SetMarkerSize(0.8);
        g->SetMarkerStyle(20);
        g->SetMarkerColor(gnuis_min->GetLineColor());
      }

      // TODO: add nuisances to legend
      // legnuis->AddEntry(gnuis_min, labelMap[thisNuis].c_str(), "L");
    }

    // Find best fit value, 1, 2, and 3 sigma values as intersections with straight line, assuming chi2
    vector<double> x_bestfit, y_bestfit;
    vector<double> x_1sigma, y_1sigma;
    vector<double> x_2sigma, y_2sigma;
    vector<double> x_3sigma, y_3sigma;
    vector<double> x_4sigma, y_4sigma;
    vector<double> x_5sigma, y_5sigma;

    vector<double> x_1sigma_wrt_0, y_1sigma_wrt_0;
    vector<double> x_2sigma_wrt_0, y_2sigma_wrt_0;

    vector<double> x_1sigma_wrt_1, y_1sigma_wrt_1;
    vector<double> x_2sigma_wrt_1, y_2sigma_wrt_1;

    for (int i_point = 0; i_point < g_min->GetN(); ++i_point) {
      double x = x_lo + i_point * (x_hi - x_lo) / g_min->GetN();
      x_bestfit.push_back(x);
      y_bestfit.push_back(0.0001);
      x_1sigma.push_back(x);
      y_1sigma.push_back(1.0);
      x_2sigma.push_back(x);
      y_2sigma.push_back(4.0);
      x_3sigma.push_back(x);
      y_3sigma.push_back(9.0);
      x_4sigma.push_back(x);
      y_4sigma.push_back(16.0);
      x_5sigma.push_back(x);
      y_5sigma.push_back(25.0);

      x_1sigma_wrt_0.push_back(x);
      x_2sigma_wrt_0.push_back(x);
      x_1sigma_wrt_1.push_back(x);
      x_2sigma_wrt_1.push_back(x);

      double y_at_1_plus_1sigma = g_min->Eval(1.0, 0, "S") + TMath::ChisquareQuantile( 0.68, 1 );
      double y_at_1_plus_2sigma = g_min->Eval(1.0, 0, "S") + TMath::ChisquareQuantile( 0.95, 1 ) ;

      double y_at_0_plus_1sigma = g_min->Eval(0.0, 0, "S") + TMath::ChisquareQuantile( 0.68, 1 );
      double y_at_0_plus_2sigma = g_min->Eval(0.0, 0, "S") + TMath::ChisquareQuantile( 0.95, 1 ) ;

      y_1sigma_wrt_0.push_back(y_at_0_plus_1sigma);
      y_2sigma_wrt_0.push_back(y_at_0_plus_2sigma);
      y_1sigma_wrt_1.push_back(y_at_1_plus_1sigma);
      y_2sigma_wrt_1.push_back(y_at_1_plus_2sigma);
    }

    int nrPoints = x_bestfit.size();
    TGraph* g_bestfit = new TGraph(nrPoints, getAry(x_bestfit), getAry(y_bestfit));
    TGraph* g_1sigma = new TGraph(nrPoints, getAry(x_1sigma), getAry(y_1sigma));
    TGraph* g_2sigma = new TGraph(nrPoints, getAry(x_2sigma), getAry(y_2sigma));
    TGraph* g_3sigma = new TGraph(nrPoints, getAry(x_3sigma), getAry(y_3sigma));
    TGraph* g_4sigma = new TGraph(nrPoints, getAry(x_4sigma), getAry(y_4sigma));
    TGraph* g_5sigma = new TGraph(nrPoints, getAry(x_5sigma), getAry(y_5sigma));

    TGraph* g_1sigma_wrt_0 = new TGraph(nrPoints, getAry(x_1sigma_wrt_0), getAry(y_1sigma_wrt_0));
    TGraph* g_2sigma_wrt_0 = new TGraph(nrPoints, getAry(x_2sigma_wrt_0), getAry(y_2sigma_wrt_0));

    TGraph* g_1sigma_wrt_1 = new TGraph(nrPoints, getAry(x_1sigma_wrt_1), getAry(y_1sigma_wrt_1));
    TGraph* g_2sigma_wrt_1 = new TGraph(nrPoints, getAry(x_2sigma_wrt_1), getAry(y_2sigma_wrt_1));

    g_1sigma->SetLineColor(color_1sigma);
    g_2sigma->SetLineColor(color_2sigma);
    g_3sigma->SetLineColor(color_3sigma);
    g_4sigma->SetLineColor(color_4sigma);
    g_5sigma->SetLineColor(color_5sigma);

    g_1sigma->SetLineStyle(style[i_overlay]);
    g_2sigma->SetLineStyle(style[i_overlay]);
    g_3sigma->SetLineStyle(style[i_overlay]);
    g_4sigma->SetLineStyle(style[i_overlay]);
    g_5sigma->SetLineStyle(style[i_overlay]);

    if (i_overlay == 0) {
      g_all->Add(g_1sigma, "L");
      g_all->Add(g_2sigma, "L");
      g_all->Add(g_3sigma, "L");
      g_all->Add(g_4sigma, "L");
      g_all->Add(g_5sigma, "L");
    }

    vector<double> x_bg, y_bg;
    vector<double> x_sm, y_sm;

    for (double scale_y = 0.0; scale_y < 49; scale_y += 0.1) {
      x_bg.push_back(0.0);
      x_sm.push_back(1.0);
      y_bg.push_back(scale_y);
      y_sm.push_back(scale_y);
    }

    int nrBgPoints = x_bg.size();
    TGraph* g_bg = new TGraph(nrBgPoints, getAry(x_bg), getAry(y_bg));
    TGraph* g_sm = new TGraph(nrBgPoints, getAry(x_sm), getAry(y_sm));

    TGraph* g_bestfit_numeric = findIntersection(*g_min, *g_bestfit);
    TGraph* g_1sigma_numeric = findIntersection(*g_min, *g_1sigma);
    TGraph* g_2sigma_numeric = findIntersection(*g_min, *g_2sigma);
    TGraph* g_3sigma_numeric = findIntersection(*g_min, *g_3sigma);

    TGraph* g_bg_numeric = findIntersection(*g_min, *g_bg);
    TGraph* g_sm_numeric = findIntersection(*g_min, *g_sm);

    TGraph* g_1sigma_wrt_0_numeric = findIntersection(*g_min, *g_1sigma_wrt_0);
    TGraph* g_2sigma_wrt_0_numeric = findIntersection(*g_min, *g_2sigma_wrt_0);

    TGraph* g_1sigma_wrt_1_numeric = findIntersection(*g_min, *g_1sigma_wrt_1);
    TGraph* g_2sigma_wrt_1_numeric = findIntersection(*g_min, *g_2sigma_wrt_1);

    g_bestfit_numeric->SetMarkerColor(color_bestfit);
    g_bestfit_numeric->SetMarkerStyle(22);

    cout << "Significance " << legend[i_overlay] << ":" << pow(g_min->Eval(0.0, 0, "S"), 0.5) << endl;

    for (int i_point = 0; i_point < g_1sigma_numeric->GetN(); ++i_point) {
      vector<double> x_vertical, y_vertical;
      double x, y = 0;
      g_1sigma_numeric->GetPoint(i_point, x, y);
      x_vertical.push_back(x);
      x_vertical.push_back(x);
      y_vertical.push_back(0.0);
      y_vertical.push_back(1.0);
      TGraph* g_vertical = new TGraph(2, getAry(x_vertical), getAry(y_vertical));
      g_vertical->SetLineColor(color_1sigma);
      g_vertical->SetLineStyle(style[i_overlay]);
      // g_all->Add(g_vertical, "L");
    }

    for (int i_point = 0; i_point < g_2sigma_numeric->GetN(); ++i_point) {
      vector<double> x_vertical, y_vertical;
      double x, y = 0;
      g_2sigma_numeric->GetPoint(i_point, x, y);
      x_vertical.push_back(x);
      x_vertical.push_back(x);
      y_vertical.push_back(0.0);
      y_vertical.push_back(4.0);
      TGraph* g_vertical = new TGraph(2, getAry(x_vertical), getAry(y_vertical));
      g_vertical->SetLineColor(color_2sigma);
      g_vertical->SetLineStyle(style[i_overlay]);
      // g_all->Add(g_vertical, "L");
    }

    for (int i_point = 0; i_point < g_3sigma_numeric->GetN(); ++i_point) {
      vector<double> x_vertical, y_vertical;
      double x, y = 0;
      g_3sigma_numeric->GetPoint(i_point, x, y);
      x_vertical.push_back(x);
      x_vertical.push_back(x);
      y_vertical.push_back(0.0);
      y_vertical.push_back(9.0);
      TGraph* g_vertical = new TGraph(2, getAry(x_vertical), getAry(y_vertical));
      g_vertical->SetLineColor(color_3sigma);
      g_vertical->SetLineStyle(style[i_overlay]);
      // g_all->Add(g_vertical, "L");
    }

    for (int i_point = 0; i_point < g_1sigma_numeric->GetN(); ++i_point) {
      double x, y = 0;
      g_1sigma_numeric->GetPoint(i_point, x, y);
      if (x > intervalBoundary_hi) intervalBoundary_hi = x;
      if (x < intervalBoundary_lo) intervalBoundary_lo = x;
    }

    for (int i_point = 1; i_point < g_min->GetN(); ++i_point) {
      double x, y = 0;
      double x_previous, y_previous = 0;
      g_min->GetPoint(i_point, x, y);
      g_min->GetPoint(i_point-1, x_previous, y_previous);
      if (y < y_previous) leftmost_minimum = x;
      else if (y > y_previous) break;
    }

    for (int i_point = g_min->GetN()-2; i_point >= 0 ; --i_point) {
      double x, y = 0;
      double x_next, y_next = 0;
      g_min->GetPoint(i_point, x, y);
      g_min->GetPoint(i_point+1, x_next, y_next);
      if (y < y_next) rightmost_minimum = x;
      else if (y > y_next) break;
    }

    map_overlay2bestfitgraph[thisOverlay.Data()] = g_bestfit_numeric;
    map_overlay21sigmagraph[thisOverlay.Data()] = g_1sigma_numeric;
    map_overlay22sigmagraph[thisOverlay.Data()] = g_2sigma_numeric;
    map_overlay23sigmagraph[thisOverlay.Data()] = g_3sigma_numeric;

    delete g_bg_numeric;
    delete g_sm_numeric;

    i_overlay++;
  }

  // Prepare boundary histogram
  TH2D* hist_boundaries = new TH2D("hist_boundaries", "hist_boundaries", 2, x_lo, x_hi, 2, y_lo, y_hi);
  TString ylabel = "-2 ln #Lambda(" + axis_label[0] + ")";
  hist_boundaries->SetTitle((";" + axis_label[0] + ";" + ylabel.Data()).c_str());
  hist_boundaries->Draw();
  gPad->Modified();

  // Draw all graphs
  TList* g_all_list = g_all->GetListOfGraphs();
  int nrGraphs = g_all_list->GetEntries();
  for (int itrChan = 0; itrChan < nrGraphs; ++itrChan) {
    TGraph* thisGraph = (TGraph*)g_all_list->At(itrChan);
    Option_t* thisOptions = g_all->GetGraphDrawOption(thisGraph);
    TString thisOptionsPlusSame = thisOptions;
    thisOptionsPlusSame.ReplaceAll("A", "");
    thisOptionsPlusSame += " SAME";
    thisGraph->DrawClone(thisOptionsPlusSame.Data());
  }

  // Cover top part of canvas for printing labels, etc.
  std::vector<double> whiteX;
  std::vector<double> whiteY;

  whiteX.push_back(x_lo); whiteY.push_back(0.725 * y_hi);
  whiteX.push_back(x_hi); whiteY.push_back(0.725 * y_hi);
  whiteX.push_back(x_hi); whiteY.push_back(y_hi);
  whiteX.push_back(x_lo); whiteY.push_back(y_hi);
  whiteX.push_back(x_lo); whiteY.push_back(0.725 * y_hi);
  TGraph* gwhite = new TGraph(5, getAry(whiteX), getAry(whiteY));
  gwhite->SetFillColor(kWhite);
  gwhite->SetLineColor(kWhite);
  gwhite->SetLineWidth(0);
  gwhite->DrawClone("F SAME");

  // Redraw axes
  TGaxis *axis1 = new TGaxis(x_lo, y_lo, x_lo, y_hi, y_lo, y_hi, 0, "-");
  axis1->SetNdivisions(hist_boundaries->GetYaxis()->GetNdivisions());
  axis1->ImportAxisAttributes(hist_boundaries->GetYaxis());
  axis1->SetTitle(hist_boundaries->GetYaxis()->GetTitle());
  axis1->SetLabelFont(42);
  axis1->SetTitleFont(42);
  axis1->SetNoExponent(kTRUE);
  axis1->SetMoreLogLabels();
  axis1->SetLineColor(kBlack);
  hist_boundaries->GetYaxis()->SetLabelSize(0);
  hist_boundaries->GetYaxis()->SetTitleSize(0);
  axis1->Draw();

  TGaxis *axis2 = new TGaxis(x_hi, y_lo, x_hi, y_hi, y_lo, y_hi, 0, "+");
  axis2->SetNdivisions(hist_boundaries->GetYaxis()->GetNdivisions());
  axis2->ImportAxisAttributes(hist_boundaries->GetYaxis());
  axis2->SetTitle("");
  axis2->SetLabelFont(42);
  axis2->SetTitleFont(42);
  axis2->SetNoExponent(kTRUE);
  axis2->SetMoreLogLabels();
  axis2->SetLineColor(kBlack);
  hist_boundaries->GetYaxis()->SetLabelSize(0);
  axis2->Draw();

  TGaxis *axis0 = new TGaxis(x_lo, y_lo, x_hi, y_lo, x_lo, x_hi, hist_boundaries->GetXaxis()->GetNdivisions(), "+");
  axis0->ImportAxisAttributes(hist_boundaries->GetXaxis());
  axis0->SetTextFont(42);
  axis0->SetLabelFont(42);
  axis0->SetTitleFont(42);
  axis0->SetLabelSize(1.0*hist_boundaries->GetXaxis()->GetLabelSize());
  axis0->SetLineColor(kBlack);
  axis0->SetTitleOffset(1.0);
  axis0->SetTitle(hist_boundaries->GetXaxis()->GetTitle());
  // axis0->SetTitle("");
  hist_boundaries->GetXaxis()->SetLabelSize(0);
  // axis0->SetNdivisions(505);
  axis0->Draw();

  TGaxis *axis10 = new TGaxis(x_lo, y_hi, x_hi, y_hi, x_lo, x_hi, hist_boundaries->GetXaxis()->GetNdivisions(), "-");
  axis10->ImportAxisAttributes(hist_boundaries->GetXaxis());
  axis10->SetTextFont(42);
  axis10->SetLabelFont(42);
  axis10->SetTitleFont(42);
  axis10->SetLabelSize(1.0*hist_boundaries->GetXaxis()->GetLabelSize());
  axis10->SetTitle("");
  axis10->SetLineColor(kBlack);
  hist_boundaries->GetXaxis()->SetLabelSize(0);
  axis10->Draw();

  c1->Update();

  // ATLAS label
  double position_label_x = 0.2;
  double position_label_y = 0.885;

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(72);
  tex.SetTextColor(kBlack);
  tex.DrawLatex(position_label_x, position_label_y, "ATLAS");

  // Internal, Preliminary, etc.
  TLatex tex_label;
  tex_label.SetNDC();
  tex_label.SetTextFont(42);
  tex_label.SetTextColor(1);
  tex_label.DrawLatex(position_label_x + 0.1275, position_label_y, label.c_str());

  // Luminosity
  TLatex tex_luminosity;
  tex_luminosity.SetNDC();
  tex_luminosity.SetTextFont(42);
  tex_luminosity.SetTextColor(1);
  tex_luminosity.DrawLatex(position_label_x, position_label_y - 1 * 0.055, luminosity.c_str());

  // Legend
  double position_legend_x = position_label_x + 0.4;
  double position_legend_y = position_label_y + 0.035 - 0.043 * leg->GetListOfPrimitives()->GetSize();

  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetX1(position_legend_x + 0.0);
  leg->SetY2(position_legend_y + 0.043 * leg->GetListOfPrimitives()->GetSize());
  leg->SetX2(position_legend_x + 0.25);
  leg->SetY1(position_legend_y - 0.0);
  leg->SetTextSize(0.034);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("F");

  if (legnuis->GetListOfPrimitives()->GetSize() > 0) {
    legnuis->SetTextFont(42);
    legnuis->SetFillColor(kWhite);
    legnuis->SetX1(0.0);
    legnuis->SetY2(1.0);
    legnuis->SetX2(1.0);
    legnuis->SetY1(0.95);
    legnuis->SetTextSize(0.025);
    legnuis->SetFillStyle(0);
    legnuis->SetBorderSize(0);
    legnuis->SetNColumns(ceil(legnuis->GetListOfPrimitives()->GetSize() / 2.0));
    legnuis->Draw("F");
  }

  // Save the canvas
  c1->Update();
  pad1->Update();

  stringstream saveName;
  saveName << "scan_" << poi;
  save(saveName.str(), "eps", c1);
  save(saveName.str(), "pdf", c1);
  save(saveName.str(), "C", c1);
  save(saveName.str(), "png", c1);
}

// ____________________________________________________________________________|__________
// Plot 2D profile likelihood scan
void plot2D(vector<string> filenames, vector<string> poi, vector<string> nuis,
            vector<pair<string, map<string, map< vector<double>, double > > > > map_overlay2folder2poi2nll,
            vector<pair<string, map<string, map< vector<double>, int > > > > map_overlay2folder2poi2status,
            vector<pair<string, map<string, map< vector<double>, map<string, double> > > > > map_overlay2folder2poi2nuis,
            double x_lo, double x_hi, double y_lo, double y_hi,
            vector<string> color, vector<int> style, vector<string> legend, vector<string> axis_label,
            string label, string luminosity)
{
  LOG(logINFO) << "Plotting 2D likelihood scan.";

  map<string, map<string, TH2D* > > map_overlay2folder2hist;
  map<string, map<string, TGraph2D* > > map_overlay2folder2graph;
  map<string, map<string, TH2D* > > map_overlay2folder2interpolatedhist;

  for (auto it_overlay : map_overlay2folder2poi2nll) {
    int i_overlay = &it_overlay - &map_overlay2folder2poi2nll[0];
    string thisOverlay = it_overlay.first;
    map<string, map< vector<double>, double > > map_folder2poi2nll = it_overlay.second;

    map<string, TH2D* > map_folder2hist;
    map<string, TGraph2D* > map_folder2graph;
    map<string, TH2D* > map_folder2interpolatedhist;

    double xlo =  numeric_limits<double>::infinity();
    double xhi = -numeric_limits<double>::infinity();
    double ylo =  numeric_limits<double>::infinity();
    double yhi = -numeric_limits<double>::infinity();

    double xlo_raw =  numeric_limits<double>::infinity();
    double xhi_raw = -numeric_limits<double>::infinity();
    double ylo_raw =  numeric_limits<double>::infinity();
    double yhi_raw = -numeric_limits<double>::infinity();

    // Order coordinates for every folder and make histogram
    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      map< vector<double>, double > map_poi2nll = it_folder->second;

      set<double> sx, sy;
      vector<double> x, y;

      for(map< vector<double>, double>::iterator it_poi = map_poi2nll.begin(); it_poi != map_poi2nll.end(); ++it_poi) {
        sx.insert(it_poi->first[0]);
        sy.insert(it_poi->first[1]);
      }

      x.assign(sx.begin(), sx.end());
      y.assign(sy.begin(), sy.end());

      if (x[0] < xlo_raw) xlo_raw = x[0];
      if (x[x.size()-1] > xhi_raw) xhi_raw = x[x.size()-1];

      if (y[0] < ylo_raw) ylo_raw = y[0];
      if (y[y.size()-1] > yhi_raw) yhi_raw = y[y.size()-1];

      int nrPoints_in_x = x.size();
      int nrPoints_in_y = y.size();

      deque<double> x_bin_boundaries, y_bin_boundaries;

      for (int i_point = 0; i_point < nrPoints_in_x - 1; i_point++) {
        double point_lo = x[i_point];
        double point_hi = x[i_point+1];
        double edge = point_lo + fabs(point_hi - point_lo) / 2.0;
        x_bin_boundaries.push_back(edge);
      }

      x_bin_boundaries.push_front(x_bin_boundaries[0] - fabs(x_bin_boundaries[1] - x_bin_boundaries[0]) / 1.0);
      x_bin_boundaries.push_front(x_bin_boundaries[0] - fabs(x_bin_boundaries[1] - x_bin_boundaries[0]) / 1.0);
      x_bin_boundaries.push_back(x_bin_boundaries[x_bin_boundaries.size()-1] + fabs(x_bin_boundaries[x_bin_boundaries.size()-1] - x_bin_boundaries[x_bin_boundaries.size()-2]) / 1.0);
      x_bin_boundaries.push_back(x_bin_boundaries[x_bin_boundaries.size()-1] + fabs(x_bin_boundaries[x_bin_boundaries.size()-1] - x_bin_boundaries[x_bin_boundaries.size()-2]) / 1.0);

      for (int i_point = 0; i_point < nrPoints_in_y - 1; i_point++) {
        double point_lo = y[i_point];
        double point_hi = y[i_point+1];
        double edge = point_lo + (point_hi - point_lo) / 2.0;
        y_bin_boundaries.push_back(edge);
      }

      y_bin_boundaries.push_front(y_bin_boundaries[0] - fabs(y_bin_boundaries[1] - y_bin_boundaries[0]) / 1.0);
      y_bin_boundaries.push_front(y_bin_boundaries[0] - fabs(y_bin_boundaries[1] - y_bin_boundaries[0]) / 1.0);
      y_bin_boundaries.push_back(y_bin_boundaries[y_bin_boundaries.size()-1] + fabs(y_bin_boundaries[y_bin_boundaries.size()-1] - y_bin_boundaries[y_bin_boundaries.size()-2]) / 1.0);
      y_bin_boundaries.push_back(y_bin_boundaries[y_bin_boundaries.size()-1] + fabs(y_bin_boundaries[y_bin_boundaries.size()-1] - y_bin_boundaries[y_bin_boundaries.size()-2]) / 1.0);

      if (x_bin_boundaries[0] < xlo) xlo = x_bin_boundaries[0];
      if (x_bin_boundaries[x_bin_boundaries.size()-1] > xhi) xhi = x_bin_boundaries[x_bin_boundaries.size()-1];

      if (y_bin_boundaries[0] < ylo) ylo = y_bin_boundaries[0];
      if (y_bin_boundaries[y_bin_boundaries.size()-1] > yhi) yhi = y_bin_boundaries[y_bin_boundaries.size()-1];

      TString histName = "hist_" + thisOverlay + "_" + thisFolder;
      TH2D *h = new TH2D(histName.Data(), histName.Data(), x_bin_boundaries.size()-1, getAry(x_bin_boundaries), y_bin_boundaries.size()-1, getAry(y_bin_boundaries));

      std::vector<double> x_graph, y_graph, z_graph;

      for(map< vector<double>, double>::iterator it_poi = map_poi2nll.begin(); it_poi != map_poi2nll.end(); ++it_poi) {
        double x_val = it_poi->first[0];
        double y_val = it_poi->first[1];
        double nll_val = it_poi->second;
        int bin = h->FindBin(x_val, y_val);

        h->SetBinContent(bin, nll_val);
        x_graph.push_back(x_val);
        y_graph.push_back(y_val);
        z_graph.push_back(nll_val);
      }

      for (int i = 1 ; i < h->GetNbinsX()+1; i++) {
        double x = h->GetXaxis()->GetBinCenter(i);
        double y_lo = h->GetYaxis()->GetBinCenter(1);
        double y_hi = h->GetYaxis()->GetBinCenter(h->GetNbinsY());
        int bin_lo = h->FindBin(x, y_lo);
        int bin_hi = h->FindBin(x, y_hi);
        double z = h->GetMaximum();
        h->SetBinContent(bin_lo, z);
        h->SetBinContent(bin_hi, z);
      }

      for (int i = 1 ; i < h->GetNbinsY()+1; i++) {
        double y = h->GetYaxis()->GetBinCenter(i);
        double x_lo = h->GetXaxis()->GetBinCenter(1);
        double x_hi = h->GetXaxis()->GetBinCenter(h->GetNbinsX());
        int bin_lo = h->FindBin(x_lo, y);
        int bin_hi = h->FindBin(x_hi, y);
        double z = h->GetMaximum();
        h->SetBinContent(bin_lo, z);
        h->SetBinContent(bin_hi, z);
      }

      int n_graph = x_graph.size();
      TGraph2D *g = new TGraph2D(n_graph, getAry(x_graph), getAry(y_graph), getAry(z_graph));

      map_folder2hist[thisFolder] = h;
      map_folder2graph[thisFolder] = g;
    }

    // Fill empty bins by 2D linear interpolation of closest filled bins
    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TH2D* h = map_folder2hist[thisFolder];
      TGraph2D* g = map_folder2graph[thisFolder];

      // find empty bins
      for (int i = 1 ; i < h->GetNbinsX()+1; i++) {
        for (int j = 1; j < h->GetNbinsY()+1; j++) {
          double x = h->GetXaxis()->GetBinCenter(i);
          double y = h->GetYaxis()->GetBinCenter(j);
          int bin = h->FindBin(x, y);
          double z = h->GetBinContent(bin);

          if (z == 0) {
            bool interpolated = false;

            z = g->Interpolate(x, y);

            if (AlmostEqualUlpsAndAbs(z, 0.0, 0.0001, 4)) z = h->GetMaximum();
            if (fabs(z) < 100.0) z = h->GetMaximum();

            h->SetBinContent(bin, z);
          }
        }
      }
    }

    // Make smooth interpolated graph for every folder in poi range, find minimum nll
    double minNll = TMath::Infinity();

    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TH2D* h = map_folder2hist[thisFolder];
      TGraph2D* g = map_folder2graph[thisFolder];

      int n_interp_bins = 50;

      deque<double> x_bin_boundaries, y_bin_boundaries;
      double xspacing = fabs(xhi - xlo) / n_interp_bins;
      double yspacing = fabs(yhi - ylo) / n_interp_bins;

      for (int i_point = 0; i_point < n_interp_bins + 1; i_point++) {
        double edge = xlo + i_point * xspacing;
        x_bin_boundaries.push_back(edge);
      }

      for (int i_point = 0; i_point < n_interp_bins + 1; i_point++) {
        double edge = ylo + i_point * yspacing;
        y_bin_boundaries.push_back(edge);
      }

      TH2D *h_interpolated = new TH2D(h->GetName(), "", x_bin_boundaries.size()-1, getAry(x_bin_boundaries), y_bin_boundaries.size()-1, getAry(y_bin_boundaries));

      for (int i = 1 ; i < h_interpolated->GetNbinsX()+1; i++) {
        for (int j = 1; j < h_interpolated->GetNbinsY()+1; j++) {
          double x = h_interpolated->GetXaxis()->GetBinCenter(i);
          double y = h_interpolated->GetYaxis()->GetBinCenter(j);

          double z = g->Interpolate(x, y);

          if (AlmostEqualUlpsAndAbs(z, 0.0, 0.0001, 4)) z = h->GetMaximum();
          if (fabs(z) < 100.0) z = h->GetMaximum();

          int bin = h_interpolated->FindBin(x, y);
          h_interpolated->SetBinContent(bin, z);
        }
      }

      // h_interpolated = IncreaseResolutionAndSmooth(h, 300, xlo, xhi, 300, ylo, yhi);
      h_interpolated = h;

      for (int i = 1 ; i < h_interpolated->GetNbinsX()+1; i++) {
        for (int j = 1; j < h_interpolated->GetNbinsY()+1; j++) {
          double x = h_interpolated->GetXaxis()->GetBinCenter(i);
          double y = h_interpolated->GetYaxis()->GetBinCenter(j);
          int bin = h_interpolated->FindBin(x, y);
          double z = h_interpolated->GetBinContent(bin);

          if (z < minNll) {
            minNll = z;
          }
        }
      }

      map_folder2interpolatedhist[thisFolder] = h_interpolated;
    }

    // Subtract minimum nll in every point of coarse and smooth graph
    for (map<string, map< vector<double>, double > >::iterator it_folder = map_folder2poi2nll.begin(); it_folder != map_folder2poi2nll.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TH2D* h_interpolated = map_folder2interpolatedhist[thisFolder];

      for (int i = 1 ; i < h_interpolated->GetNbinsX()+1; i++) {
        for (int j = 1; j < h_interpolated->GetNbinsY()+1; j++) {
          double x = h_interpolated->GetXaxis()->GetBinCenter(i);
          double y = h_interpolated->GetYaxis()->GetBinCenter(j);
          int bin = h_interpolated->FindBin(x, y);
          double z = h_interpolated->GetBinContent(bin);
          z -= minNll;

          h_interpolated->SetBinContent(bin, z);
        }
      }

      // h_interpolated = IncreaseResolutionAndSmooth(h_interpolated, 50, xlo, xhi, 50, ylo, yhi);
      map_folder2interpolatedhist[thisFolder] = h_interpolated;
    }

    map_overlay2folder2interpolatedhist[thisOverlay] = map_folder2interpolatedhist;
  }

  map<string, TH2D* > map_overlay2bestinterpolatedhist;

  for(map<string, map<string, TH2D* > >::iterator it_overlay = map_overlay2folder2interpolatedhist.begin(); it_overlay != map_overlay2folder2interpolatedhist.end(); ++it_overlay) {
    string thisOverlay = it_overlay->first;
    map<string, TH2D* > map_folder2interpolatedhist = it_overlay->second;

    bool isFirst = true;
    TH2D* h_min;

    for (map<string, TH2D* >::iterator it_folder = map_folder2interpolatedhist.begin(); it_folder != map_folder2interpolatedhist.end(); ++it_folder) {
      TH2D* h_interpolated = it_folder->second;
      if (isFirst) {
        h_min = (TH2D*)h_interpolated->Clone();
      } else {
        for (int i = 1 ; i < h_min->GetNbinsX()+1; i++) {
          for (int j = 1; j < h_min->GetNbinsY()+1; j++) {
            double x = h_min->GetXaxis()->GetBinCenter(i);
            double y = h_min->GetYaxis()->GetBinCenter(j);
            int bin = h_min->FindBin(x, y);
            double z_min = h_min->GetBinContent(bin);
            double z = h_interpolated->GetBinContent(bin);

            if (z < z_min) {
              z_min = z;
            }

            h_min->SetBinContent(bin, z_min);
          }
        }
      }
    }

    map_overlay2bestinterpolatedhist[thisOverlay] = h_min;
  }

  // Find all contours
  map<string, map<int, deque<TGraph*> > > map_overlay2CL2bestgraph;
  map<string, TGraph* > map_overlay2bestfit;
  map<string, map<string, map<int, deque<TGraph*> > > > map_overlay2folder2CL2graph;

  for(map<string, map<string, TH2D* > >::iterator it_overlay = map_overlay2folder2interpolatedhist.begin(); it_overlay != map_overlay2folder2interpolatedhist.end(); ++it_overlay) {
    string thisOverlay = it_overlay->first;
    map<string, TH2D* > map_folder2interpolatedhist = it_overlay->second;

    TH2D* h_min = map_overlay2bestinterpolatedhist[thisOverlay];
    map<int, deque<TGraph*> > map_CL2bestgraph = GetContours(h_min);
    map_overlay2CL2bestgraph[thisOverlay] = map_CL2bestgraph;
    TGraph* bestfitPoint = GetBestFit(h_min);
    map_overlay2bestfit[thisOverlay] = bestfitPoint;

    for (map<string, TH2D* >::iterator it_folder = map_folder2interpolatedhist.begin(); it_folder != map_folder2interpolatedhist.end(); ++it_folder) {
      string thisFolder = it_folder->first;
      TH2D* h_interpolated = it_folder->second;
      map<int, deque<TGraph*> > map_CL2graph = GetContours(h_interpolated);
      map_overlay2folder2CL2graph[thisOverlay][thisFolder] = map_CL2graph;
    }
  }

  // Final 2D plotting
  TCanvas* c1 = new TCanvas("c1", "c1", 1024, 768);
  TPad *pad1 = new TPad("pad1", "pad1", 0.0 , 0.0 , 1.0 , 1.0 , 0);
  pad1->Draw();
  pad1->cd();

  TMultiGraph* g_all = new TMultiGraph();
  TMultiGraph* g_global = new TMultiGraph();
  TLegend* leg = new TLegend(0.0, 0.0, 0.0, 0.0, "", "NDC");

  // Prepare boundary histogram
  TH2D* hist_boundaries = new TH2D("hist_boundaries", "hist_boundaries", 2, x_lo, x_hi, 2, y_lo, y_hi);
  hist_boundaries->SetTitle((";" + axis_label[0] + ";" + axis_label[1]).c_str());
  hist_boundaries->Draw();
  gPad->Modified();

  // Draw all graphs
  for (map<string, map<string, map<int, deque<TGraph*> > > >::iterator it_overlay = map_overlay2folder2CL2graph.begin(); it_overlay != map_overlay2folder2CL2graph.end(); ++it_overlay) {
    TString thisOverlay = it_overlay->first;
    int i_overlay = std::distance(map_overlay2folder2CL2graph.begin(), it_overlay);

    string legendText = legend[i_overlay];
    Color_t mycolor = TColor::GetColor(color[i_overlay].c_str());

    map<int, deque<TGraph*> > map_CL2bestgraph = map_overlay2CL2bestgraph[thisOverlay.Data()];

    for(map<int, deque<TGraph*> >::iterator it_cl = map_CL2bestgraph.begin(); it_cl != map_CL2bestgraph.end(); ++it_cl) {
      int thisCL = it_cl->first;
      deque<TGraph*> bestgraph = it_cl->second;
      TGraph* bestfitPoint = map_overlay2bestfit[thisOverlay.Data()];
      bestfitPoint->SetMarkerStyle(34);
      bestfitPoint->SetMarkerSize(2.0);
      bestfitPoint->SetMarkerColor(mycolor);

      g_all->Add(bestfitPoint, "P");
      g_global->Add(bestfitPoint, "P");

      for (int i_graph = 0; i_graph < bestgraph.size(); i_graph++) {
        TGraph* g = bestgraph[i_graph];

        g->SetLineWidth(3);
        g->SetLineColor(mycolor);
        g->SetLineStyle(style[i_overlay]);

        if (thisCL == 1) {
          g->SetLineStyle(kSolid);
          g_all->Add(g, "l");
          leg->AddEntry(g, legendText.c_str(), "L");
        }

        if (thisCL == 2) {
          g->SetLineStyle(kDashed);
          g_all->Add(g, "l");
        }
      }
    }
  }

  g_all->Draw("L");

  TList* g_global_list = g_global->GetListOfGraphs();
  int nrGraphs = g_global_list->GetEntries();

  for (int itrChan = 0; itrChan < nrGraphs; ++itrChan) {
    TGraph* thisGraph = (TGraph*)g_global_list->At(itrChan);
    Option_t* thisOptions = g_global->GetGraphDrawOption(thisGraph);
    TString thisOptionsPlusSame = thisOptions;
    thisOptionsPlusSame.ReplaceAll("A", "");
    thisOptionsPlusSame += " SAME";
    thisGraph->DrawClone(thisOptionsPlusSame.Data());
  }

  // Cover top part of canvas for printing labels, etc.
  std::vector<double> whiteX;
  std::vector<double> whiteY;

  whiteX.push_back(x_lo); whiteY.push_back(0.725 * y_hi);
  whiteX.push_back(x_hi); whiteY.push_back(0.725 * y_hi);
  whiteX.push_back(x_hi); whiteY.push_back(y_hi);
  whiteX.push_back(x_lo); whiteY.push_back(y_hi);
  whiteX.push_back(x_lo); whiteY.push_back(0.725 * y_hi);
  TGraph* gwhite = new TGraph(5, getAry(whiteX), getAry(whiteY));
  gwhite->SetFillColor(kWhite);
  gwhite->SetLineColor(kWhite);
  gwhite->SetLineWidth(0);
  gwhite->DrawClone("F SAME");

  // Redraw axes
  TGaxis *axis1 = new TGaxis(x_lo, y_lo, x_lo, y_hi, y_lo, y_hi, 0, "-");
  axis1->SetNdivisions(hist_boundaries->GetYaxis()->GetNdivisions());
  axis1->ImportAxisAttributes(hist_boundaries->GetYaxis());
  axis1->SetTitle(hist_boundaries->GetYaxis()->GetTitle());
  axis1->SetLabelFont(42);
  axis1->SetTitleFont(42);
  axis1->SetNoExponent(kTRUE);
  axis1->SetMoreLogLabels();
  axis1->SetLineColor(kBlack);
  hist_boundaries->GetYaxis()->SetLabelSize(0);
  hist_boundaries->GetYaxis()->SetTitleSize(0);
  axis1->Draw();

  TGaxis *axis2 = new TGaxis(x_hi, y_lo, x_hi, y_hi, y_lo, y_hi, 0, "+");
  axis2->SetNdivisions(hist_boundaries->GetYaxis()->GetNdivisions());
  axis2->ImportAxisAttributes(hist_boundaries->GetYaxis());
  axis2->SetTitle("");
  axis2->SetLabelFont(42);
  axis2->SetTitleFont(42);
  axis2->SetNoExponent(kTRUE);
  axis2->SetMoreLogLabels();
  axis2->SetLineColor(kBlack);
  hist_boundaries->GetYaxis()->SetLabelSize(0);
  axis2->Draw();

  TGaxis *axis0 = new TGaxis(x_lo, y_lo, x_hi, y_lo, x_lo, x_hi, hist_boundaries->GetXaxis()->GetNdivisions(), "+");
  axis0->ImportAxisAttributes(hist_boundaries->GetXaxis());
  axis0->SetTextFont(42);
  axis0->SetLabelFont(42);
  axis0->SetTitleFont(42);
  axis0->SetLabelSize(1.0*hist_boundaries->GetXaxis()->GetLabelSize());
  axis0->SetLineColor(kBlack);
  axis0->SetTitleOffset(1.0);
  axis0->SetTitle(hist_boundaries->GetXaxis()->GetTitle());
  // axis0->SetTitle("");
  hist_boundaries->GetXaxis()->SetLabelSize(0);
  // axis0->SetNdivisions(505);
  axis0->Draw();

  TGaxis *axis10 = new TGaxis(x_lo, y_hi, x_hi, y_hi, x_lo, x_hi, hist_boundaries->GetXaxis()->GetNdivisions(), "-");
  axis10->ImportAxisAttributes(hist_boundaries->GetXaxis());
  axis10->SetTextFont(42);
  axis10->SetLabelFont(42);
  axis10->SetTitleFont(42);
  axis10->SetLabelSize(1.0*hist_boundaries->GetXaxis()->GetLabelSize());
  axis10->SetTitle("");
  axis10->SetLineColor(kBlack);
  hist_boundaries->GetXaxis()->SetLabelSize(0);
  axis10->Draw();

  c1->Update();

  // ATLAS label
  double position_label_x = 0.2;
  double position_label_y = 0.885;

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(72);
  tex.SetTextColor(kBlack);
  tex.DrawLatex(position_label_x, position_label_y, "ATLAS");

  // Internal, Preliminary, etc.
  TLatex tex_label;
  tex_label.SetNDC();
  tex_label.SetTextFont(42);
  tex_label.SetTextColor(1);
  tex_label.DrawLatex(position_label_x + 0.1275, position_label_y, label.c_str());

  // Luminosity
  TLatex tex_luminosity;
  tex_luminosity.SetNDC();
  tex_luminosity.SetTextFont(42);
  tex_luminosity.SetTextColor(1);
  tex_luminosity.DrawLatex(position_label_x, position_label_y - 1 * 0.055, luminosity.c_str());

  // Legend
  double position_legend_x = position_label_x + 0.4;
  double position_legend_y = position_label_y + 0.035 - 0.043 * leg->GetListOfPrimitives()->GetSize();

  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetX1(position_legend_x + 0.0);
  leg->SetY2(position_legend_y + 0.043 * leg->GetListOfPrimitives()->GetSize());
  leg->SetX2(position_legend_x + 0.25);
  leg->SetY1(position_legend_y - 0.0);
  leg->SetTextSize(0.034);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("F");

  // General legend

  TLegend* leg_general = new TLegend(0.0, 0.0, 0.0, 0.0, "", "NDC");
  leg_general->SetNColumns(2);
  leg_general->SetTextFont(42);
  leg_general->SetFillColor(kWhite);
  leg_general->SetX1(0.195);
  leg_general->SetY2(0.275);
  leg_general->SetX2(0.5);
  leg_general->SetY1(0.175);
  leg_general->SetTextSize(0.034);
  leg_general->SetFillStyle(0);
  leg_general->SetBorderSize(0);

  const Int_t n_leg = 1;
  Double_t x0_leg[n_leg], y0_leg[n_leg];
  Double_t x1_leg[n_leg], y1_leg[n_leg];

  x0_leg[0] = 1;
  y0_leg[0] = 1;
  x1_leg[0] = 1;
  y1_leg[0] = 1;

  TGraph* g0_leg = new TGraph(n_leg, x0_leg, y0_leg);
  TGraph* g1_leg = new TGraph(n_leg, x1_leg, y1_leg);

  g0_leg->SetMarkerStyle(29);
  g0_leg->SetMarkerSize(3);
  g0_leg->SetLineStyle(kSolid);
  g0_leg->SetLineWidth(3);
  g0_leg->SetMarkerSize(2.5);

  g1_leg->SetMarkerStyle(34);
  g1_leg->SetMarkerSize(3);
  g1_leg->SetLineStyle(kDashed);
  g1_leg->SetLineWidth(3);
  g1_leg->SetMarkerSize(2.0);

  leg_general->AddEntry(g0_leg, "SM", "P");
  leg_general->AddEntry(g0_leg, "68% CL", "L");
  leg_general->AddEntry(g1_leg, "Best fit", "P");
  leg_general->AddEntry(g1_leg, "95% CL", "L");

  g0_leg->Draw("PSAME");
  leg_general->Draw("F");

  // Save the canvas
  c1->Update();
  pad1->Update();

  stringstream saveName;
  saveName << "scan_" << poi[0] << "_" << poi[1];
  save(saveName.str(), "eps", c1);
  save(saveName.str(), "pdf", c1);
  save(saveName.str(), "C", c1);
  save(saveName.str(), "png", c1);
}


// _____________________________________________________________________________
TH2D* IncreaseResolutionAndSmooth(TH2D* h, int xbins, double xlo, double xhi,
                                  int ybins, double ylo, double yhi)
{
  deque<double> x_bin_boundaries, y_bin_boundaries;
  double xspacing = fabs(xhi - xlo) / xbins;
  double yspacing = fabs(yhi - ylo) / ybins;

  for (int i_point = 0; i_point < xbins + 1; i_point++) {
    double edge = xlo + i_point * xspacing;
    x_bin_boundaries.push_back(edge);
  }

  for (int i_point = 0; i_point < ybins + 1; i_point++) {
    double edge = ylo + i_point * yspacing;
    y_bin_boundaries.push_back(edge);
  }

  RooWorkspace* ws = new RooWorkspace();

  stringstream lookup, rawx, rawy;
  lookup << "{lookup_x[" << xlo << "," << xhi << "],lookup_y[" << ylo << "," << yhi << "]}";
  rawx << "expr::func_x('@0',x[" << xlo << "," << xhi << "])";
  rawx << "expr::func_y('@0',y[" << ylo << "," << yhi << "])";

  ws->factory(lookup.str().c_str()) ;
  RooRealVar* lookup_x = ws->var("lookup_x");
  RooRealVar* lookup_y = ws->var("lookup_y");
  RooDataHist dh("dh", "dh", RooArgList(*lookup_x, *lookup_y), h);
  ws->factory(rawx.str().c_str());
  ws->factory(rawy.str().c_str());
  ws->import(dh);
  ws->factory("HistFunc::z_lookup({x,y},{lookup_x,lookup_y},dh,2)");
  RooHistFunc* hf = (RooHistFunc*)ws->obj("z_lookup");

  TH2D *h_interpolated = new TH2D(h->GetName(), "", x_bin_boundaries.size()-1, getAry(x_bin_boundaries), y_bin_boundaries.size()-1, getAry(y_bin_boundaries));

  for (int i = 1 ; i < h_interpolated->GetNbinsX()+1; i++) {
    for (int j = 1; j < h_interpolated->GetNbinsY()+1; j++) {
      double x = h_interpolated->GetXaxis()->GetBinCenter(i);
      double y = h_interpolated->GetYaxis()->GetBinCenter(j);

      ws->var("x")->setVal(x);
      ws->var("y")->setVal(y);

      double z = hf->getVal();

      int bin = h_interpolated->FindBin(x, y);
      h_interpolated->SetBinContent(bin, z);
    }
  }

  return h_interpolated;
}


// _____________________________________________________________________________
map<int, deque<TGraph*> > GetContours(TH2D* h)
{
  map<int, deque<TGraph*> > map_CL2graph;

  double def1s = TMath::ChisquareQuantile( 0.68, 2 );
  double def2s = TMath::ChisquareQuantile( 0.95, 2 );
  double def3s = TMath::ChisquareQuantile( 0.997, 2 );

  // find and draw n sigma contours on top
  double contours[3];
  contours[0] = def1s;
  contours[1] = def2s;
  contours[2] = def3s;

  stringstream ss;
  ss << "tmpc_" << h->GetName();

  TCanvas* tmpc = new TCanvas(ss.str().c_str(), "Contour List", 0, 0, 600, 600);
  TH2D *h3 = new TH2D(*static_cast<TH2D*>(h));

  h3->SetContour(3, contours);

  tmpc->cd();
  h3->Draw("CONT LIST");
  tmpc->Update();

  // Get Contours
  TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  TList* contLevel = NULL;
  TGraph* curv     = NULL;
  TGraph* gc       = NULL;

  Int_t nGraphs    = 0;
  Int_t TotalConts = 0;

  if (conts == NULL) {
    cout << "WARNING::No Contours Were Extracted!" << endl;
    TotalConts = 0;
    return map_CL2graph;
  } else {
    TotalConts = conts->GetSize();
  }

  for(int i = 0; i < TotalConts; i++){
    contLevel = (TList*)conts->At(i);
    nGraphs += contLevel->GetSize();
  }

  nGraphs = 0;
  double x0, y0, z0;

  for(int i = 0; i < TotalConts; i++){
    contLevel = (TList*)conts->At(i);
    if (i == 0) z0 = contours[2];
    if (i == 1) z0 = contours[1];
    if (i == 2) z0 = contours[0];

    // Get first graph from list on curves on this level
    curv = (TGraph*)contLevel->First();
    for(int icont = 0; icont < contLevel->GetSize(); icont++) {
      curv->GetPoint(0, x0, y0);

      gc = (TGraph*)curv->Clone();

      if (z0 == def1s) {
        map_CL2graph[3].push_back(gc);
      } else if (z0 == def2s) {
        map_CL2graph[2].push_back(gc);
      } else if (z0 == def3s) {
        map_CL2graph[1].push_back(gc);
      }

      nGraphs++;

      curv = (TGraph*)contLevel->After(curv); // Get Next graph
    }
  }

  return map_CL2graph;
}


// _____________________________________________________________________________
TGraph* GetBestFit(TH2D* h)
{
  TGraph* bestfitPoint;

  double minNll = TMath::Infinity();
  const Int_t n = 1;
  Double_t x0[n], y0[n];

  for (int i = 1 ; i < h->GetNbinsX()+1; i++) {
    for (int j = 1; j < h->GetNbinsY()+1; j++) {
      double x = h->GetXaxis()->GetBinCenter(i);
      double y = h->GetYaxis()->GetBinCenter(j);
      int bin = h->FindBin(x, y);
      double z = h->GetBinContent(bin);

      h->SetBinContent(bin, z);

      if (z < minNll) {
        minNll = z;
        x0[0] = x;
        y0[0] = y;
      }
    }
  }

  bestfitPoint = new TGraph(n, x0, y0);

  return bestfitPoint;
}
