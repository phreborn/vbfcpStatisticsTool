// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2018-05-18
// Description : Plot log10(S/B)

#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLine.h"
#include "TSystem.h"
#include "TTime.h"

#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooMultiVarGaussian.h"
#include "RooNLLVar.h"
#include "RooRandom.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRegexp.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "Math/MinimizerOptions.h"
#include "RooBifurGauss.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooPoisson.h"

#include "ExtendedMinimizer.hxx"
#include "ExtendedModel.hxx"
#include "RooMultiVarGaussianHighPrecision.h"

#include "TextTable.hxx"
#include "log.hxx"
#include "utils.hxx"

#include "atlasrootstyle/AtlasLabels.h"
#include "atlasrootstyle/AtlasStyle.h"
#include "atlasrootstyle/AtlasUtils.h"

#include <chrono>
#include <iomanip>
#include <list>
#include <math.h>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::reverse
#include <vector>       // std::vector

#if defined __GNUC__ || defined __APPLE__
#include <Eigen/Dense>
#else
#include <eigen3/Eigen/Dense>
#endif

#pragma warning(push)
#pragma warning(disable : 4819)
#include "boost/program_options.hpp"
#include "boost/program_options/cmdline.hpp"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/variables_map.hpp"
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/format.hpp>
#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>
#include <boost/range/adaptor/reversed.hpp>
#pragma warning(pop)
#pragma warning(disable : 4503)

#include "yaml-cpp/yaml.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace Eigen;
using boost::multiprecision::cpp_dec_float_50;

// _____________________________________________________________________________
// Declarations of functions used in this file
RooCurve *createBand(RooCurve *cenCurve, RooAbsPdf *pdf, RooRealVar *x, RooFitResult *res, RooArgSet *params, double Z = 1, bool isHist = false, bool mergeBins = false, double edgeLo = -999, double edgeHi = -999);
void releaseTheGhost(RooDataSet *obsdata, RooRealVar *x, RooRealVar *w, double ghostwt);
RooCurve *convertAsimovToCurve(RooDataSet *asimov, TString observableName, bool isHist, bool mergeBins = false, double epsilon = 1e-6);
void copyDataSet(RooDataSet *datadebit, RooRealVar *obsdebit, RooDataSet *datacredit, RooRealVar *obscredit, RooRealVar *wt, double scale = 1, bool verbose = false);
vector<double> makeHist(TString name, RooCurve *curve);
vector<vector<double>> makeHistBand(TString name, RooCurve *curve, RooCurve *band);
double calcPoissonCLLower(double q, double obs);
double calcPoissonCLUpper(double q, double obs);
RooAbsPdf *createHessePdf(RooFitResult *res, const RooArgSet &params);

typedef boost::multiprecision::cpp_dec_float<50> mp_backend;
typedef boost::multiprecision::number<mp_backend, boost::multiprecision::et_off> SuperFloat;
typedef Eigen::Matrix<cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixR;

// _____________________________________________________________________________
class bin {
public:
  double d;
  double b;
  double b0;
  double s;
  double s1;
  double sb;
  double log_sb;
  double b_hi;
  double b_lo;
  string channel;
};

vector<bin> getBinnedValues(RooAbsPdf *pdfi, RooDataSet *datai, RooRealVar *x, RooRealVar *poi, bool bandFromToys, RooFitResult *fr, RooFitResult *fr_bkg, double muhat, RooArgSet *nuis, double edgeLo = -999, double edgeHi = -999);

// _____________________________________________________________________________
// Main routine
int main(int argc, char **argv) {
  TTime thistime = gSystem->Now();

  // Load custom classes
  loadCustom();

  // ATLAS style
  SetAtlasStyle();

  // Model information
  string inFileName = "path/to/workspace.root";
  string wsName = "combined";
  string modelConfigName = "ModelConfig";
  string dataName = "combData";
  string snapshot = "ucmles";
  string folder = "test";

  // Parameter settings
  string poiName = "mu";

  // Misc settings
  string loglevel = "INFO";
  int fixCache = 1;
  int fixMulti = 1;
  int binnedLikelihood = 1;
  string fitresult = "path/to/fitresult.root";
  string fitresult_bkg = "path/to/fitresult.root";
  string fitresultName = "fitresult_minimizer_combData";
  bool bandFromToys = false;
  bool showOrigin = false;
  string bandFileName = "";
  int nTargetBins = 15;

  Color_t myOrange = TColor::GetColor("#E08214");
  Color_t myRed = TColor::GetColor("#991616");

  // Bookkeeping

  using namespace boost;
  namespace po = boost::program_options;
  po::options_description desc("Program options");
  desc.add_options()
  ("help           , h", "Print this help message")
  ("input", po::value<string>(&inFileName)->default_value(inFileName), "File to run over.")
  ("result", po::value<string>(&fitresult)->default_value(fitresult), "Path to mu=muhat fitresult.")
  ("result_bkg", po::value<string>(&fitresult_bkg)->default_value(fitresult_bkg), "Path to mu=0 fitresult.")
  ("band", po::value<string>(&bandFileName)->default_value(bandFileName), "Path to band.")
  ("resultName", po::value<string>(&fitresultName)->default_value(fitresultName), "Name of the fitresult.")
  ("poi", po::value<string>(&poiName)->default_value(poiName), "POIs to measure.")
  ("snapshot", po::value<string>(&snapshot)->default_value(snapshot), "Initial snapshot.")
  ("folder", po::value<string>(&folder)->default_value(folder), "Output folder.")
  ("workspace", po::value<string>(&wsName)->default_value(wsName), "WS to grab.")
  ("modelconfig", po::value<string>(&modelConfigName)->default_value(modelConfigName), "MC to load.")
  ("data", po::value<string>(&dataName)->default_value(dataName), "Data to use.")
  ("loglevel", po::value<string>(&loglevel)->default_value(loglevel), "Control verbosity.")
  ("nTargetBins", po::value<int>(&nTargetBins)->default_value(nTargetBins), "Number of bins in final plot.")
  ("bandFromToys", po::bool_switch(&bandFromToys), "Get uncertainty band from toys.")
  ("showOrigin", po::bool_switch(&showOrigin), "Show third panel with origin.")
  ;

  po::variables_map vm0;

  try {
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm0);
    po::notify(vm0);
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

  // Parse options
  vector<string> parsed = parseString(poiName, ",");

  // Load the model
  ExtendedModel *model = new ExtendedModel("model", inFileName, wsName, modelConfigName, dataName, snapshot, binnedLikelihood, "pdf_", fixCache, fixMulti);

  RooWorkspace *ws = model->GetWorkspace();
  ModelConfig *mc = model->GetModelConfig();
  RooAbsPdf *pdf = model->GetPdf();
  RooAbsData *data = model->GetData();
  RooArgSet *nuis = model->GetNuisanceParameters();
  RooArgSet *globs = model->GetGlobalObservables();
  RooArgSet *pois = model->GetParametersOfInterest();
  RooArgSet *obs = model->GetObservables();
  RooRealVar *poi = ws->var(poiName.c_str());

  TFile *resultFile = new TFile(fitresult.c_str());
  RooFitResult *fr = (RooFitResult *)(resultFile->Get(fitresultName.c_str()));

  double muhat = ((RooRealVar *)fr->floatParsFinal().find(poiName.c_str()))->getVal();

  TFile *resultFile_bkg = new TFile(fitresult_bkg.c_str());
  RooFitResult *fr_bkg = (RooFitResult *)(resultFile_bkg->Get(fitresultName.c_str()));

  std::system(("mkdir -vp root-files/" + folder).c_str());

  vector<RooAbsArg *> binned;
  RooFIter iter = ws->components().fwdIterator();
  RooAbsArg *arg;
  while ((arg = iter.next())) {
    if (arg->IsA() == RooRealSumPdf::Class()) {
      binned.push_back(arg);
    }
  }

  vector<bin> bins;

  RooAbsCategoryLValue *m_cat = (RooAbsCategoryLValue *)(&((RooSimultaneous *)pdf)->indexCat());
  int numChannels = m_cat->numBins(0);
  TList *m_dataList = data->split(*m_cat, true);

  // harvest bin origins
  vector<string> channels;

  for (int i = 0; i < numChannels; i++) {
    m_cat->setBin(i);
    TString channelname = m_cat->getLabel();
    cout << channelname.Data() << endl;

    if (channelname.Contains("channel_HZZ")) {
      channels.push_back("H4l");
    } else if (channelname.Contains("channel_HGam")) {
      channels.push_back("Hyy");
    } else if (channelname.Contains("channel_ttHbbdil") || channelname.Contains("channel_ttHbbljets")) {
      channels.push_back("Hbb");
    } else if (channelname.Contains("channel_ttHMLnotau") || channelname.Contains("channel_ttHMLtau")) {
      channels.push_back("HML");
    } else {
      channels.push_back("other");
    }
  }

  if (bandFileName == "") {
    for (int i = 0; i < numChannels; i++) {
      m_cat->setBin(i);
      TString channelname = m_cat->getLabel();

      string channel = "";

      if (channelname.Contains("channel_HZZ")) {
        channel = "H4l";
        channels.push_back("H4l");
      } else if (channelname.Contains("channel_HGam")) {
        channel = "Hyy";
        channels.push_back("Hyy");
      } else if (channelname.Contains("channel_ttHbbdil") || channelname.Contains("channel_ttHbbljets")) {
        channel = "Hbb";
        channels.push_back("Hbb");
      } else if (channelname.Contains("channel_ttHMLnotau") || channelname.Contains("channel_ttHMLtau")) {
        channel = "HML";
        channels.push_back("HML");
      } else {
        channel = "other";
        channels.push_back("other");
      }

      RooAbsPdf *pdfi = ((RooSimultaneous *)pdf)->getPdf(m_cat->getLabel());

      bool isBinned = false;
      for (auto p : binned) {
        if (pdfi->dependsOn(*p)) {
          isBinned = true;
          break;
        }
      }

      RooDataSet *datai = (RooDataSet *)(m_dataList->At(i));
      RooRealVar *x = (RooRealVar *)pdfi->getObservables(datai)->first();

      if (!isBinned) {
        LOG(logINFO) << "Category " << channelname << " is not binned and requires special treatment";

        vector<string> diphotonCats = {"LHCP_ttH_c33", "LHCP_ttH_c32", "LHCP_ttH_c31", "LHCP_ttH_c30", "LHCP_ttH_c29", "LHCP_ttH_c28", "LHCP_ttH_c27"};
        bool isDiphoton = false;
        for (auto c : diphotonCats) {
          if (channelname.Contains(c)) {
            isDiphoton = true;
            break;
          }
        }

        if (isDiphoton) {
          int maxtrys = 50;
          double tolerance = 0.000000001;
          double massRangeLo = 105.;
          double massRangeHi = 160.;
          double CLbound = 90;
          string sigModPrefix = "pdf_ttH_";
          string sbModPrefix = "modelSB_";

          double edgeLo = -999;
          double edgeHi = -999;

          TString thisChannelName = channelname.ReplaceAll("combCat_", "");

          string meanName = string("muCB_") + thisChannelName.Data();
          string loName = string("alphaCBLo_") + thisChannelName.Data();
          string hiName = string("alphaCBHi_") + thisChannelName.Data();

          double mean = ((RooFormulaVar *)ws->obj(meanName.c_str()))->getVal();
          double alphaLo = ((RooFormulaVar *)ws->obj(loName.c_str()))->getVal();
          double alphaHi = ((RooFormulaVar *)ws->obj(hiName.c_str()))->getVal();

          x->setRange("PartRange", massRangeLo, massRangeHi);
          x->Print("v");

          string sigmodName = sigModPrefix + thisChannelName.Data();
          RooAbsPdf *pdf_yy = (RooAbsPdf *)ws->obj(sigmodName.c_str());

          RooAbsReal *integ = pdf_yy->createIntegral(RooArgSet(*x), RooArgSet(*x), "PartRange");

          LOG(logINFO) << "Searching upper edge ...";
          double xmin = mean;
          double xmax = massRangeHi;
          int itry = 0;

          LOG(logDEBUG) << "Integral: " << integ->getVal();

          while (true) {
            edgeHi = (xmax + xmin) * 0.5;
            x->setRange("PartRange", mean, edgeHi);
            double quant = integ->getVal();
            double diff = (quant - 0.5 * 0.01 * CLbound);
            LOG(logINFO) << "edgeHi = " << edgeHi << " quantile = " << quant << " (" << diff << "%)";

            if (fabs(diff) < tolerance) {
              break;
            }

            if (diff < 0) {
              xmin = edgeHi;
            } else {
              xmax = edgeHi;
            }

            itry++;

            if (itry > maxtrys) {
              LOG(logWARNING) << "Did not converge in " << maxtrys << " attempts. Keep last high edge.";
              break;
            }
          }

          LOG(logINFO) << "Searching lower edge ...";
          xmin = massRangeLo;
          xmax = mean;
          itry = 0;

          while (true) {
            edgeLo = (xmax + xmin) * 0.5;
            x->setRange("PartRange", edgeLo, mean);
            double quant = integ->getVal();
            double diff = (quant - 0.5 * 0.01 * CLbound);
            LOG(logINFO) << "edgeLo = " << edgeLo << " quantile = " << quant << " (" << diff << "%)";

            if (fabs(diff) < tolerance) {
              break;
            }

            if (diff > 0) {
              xmin = edgeLo;
            } else {
              xmax = edgeLo;
            }

            itry++;

            if (itry > maxtrys) {
              LOG(logWARNING) << "Did not converge in " << maxtrys << " attempts. Keep last low edge.";
              break;
            }
          }

          LOG(logINFO) << "Found range [ " << edgeLo << " , " << edgeHi << " ] " << thisChannelName.Data();
          LOG(logINFO) << "sigma_" << CLbound << ": " << (edgeHi - edgeLo) / 2.;

// if (thisChannelName.Contains("LHCP_ttH_c33")) {
//   edgeLo =  122.9015;
//   edgeHi = 128.5467;
// }
// if (thisChannelName.Contains("LHCP_ttH_c32")) {
//   edgeLo =  122.3759;
//   edgeHi = 128.9636;
// }
// if (thisChannelName.Contains("LHCP_ttH_c31")) {
//   edgeLo =  122.1479;
//   edgeHi = 129.2057;
// }
// if (thisChannelName.Contains("LHCP_ttH_c30")) {
//   edgeLo =  122.9197;
//   edgeHi = 128.5586;
// }
// if (thisChannelName.Contains("LHCP_ttH_c29")) {
//   edgeLo =  122.5933;
//   edgeHi = 128.8744;
// }
// if (thisChannelName.Contains("LHCP_ttH_c28")) {
//   edgeLo =  122.5123;
//   edgeHi = 128.9521;
// }
// if (thisChannelName.Contains("LHCP_ttH_c27")) {
//   edgeLo =  122.5274;
//   edgeHi = 128.9206;
// }

          vector<bin> thisBins = getBinnedValues(pdfi, datai, x, poi, bandFromToys, fr, fr_bkg, muhat, nuis, edgeLo, edgeHi);
          for (auto thisBin : thisBins) {
            thisBin.channel = channel;
            bins.push_back(thisBin);
          }

        } else {
          // H4l
          vector<bin> thisBins = getBinnedValues(pdfi, datai, x, poi, bandFromToys, fr, fr_bkg, muhat, nuis);
          for (auto thisBin : thisBins) {
            thisBin.channel = channel;
            bins.push_back(thisBin);
          }
        }
      } else {
        // ML and bb
        vector<bin> thisBins = getBinnedValues(pdfi, datai, x, poi, bandFromToys, fr, fr_bkg, muhat, nuis);
        for (auto thisBin : thisBins) {
          thisBin.channel = channel;
          bins.push_back(thisBin);
        }
      }

      // write bins to file
      TFile fout("debug_bands_playing.root", "recreate");

      TTree *bandTree = new TTree("band", "band");
      bandTree->SetDirectory(0);

      bin this_bin;
      this_bin.d = -999.;
      this_bin.b = -999.;
      this_bin.b0 = -999.;
      this_bin.s = -999.;
      this_bin.s1 = -999.;
      this_bin.sb = -999.;
      this_bin.log_sb = -999.;
      this_bin.b_hi = -999.;
      this_bin.b_lo = -999.;
      this_bin.channel = "";

      bandTree->Branch("d", &this_bin.d);
      bandTree->Branch("b", &this_bin.b);
      bandTree->Branch("b0", &this_bin.b0);
      bandTree->Branch("s", &this_bin.s);
      bandTree->Branch("s1", &this_bin.s1);
      bandTree->Branch("sb", &this_bin.sb);
      bandTree->Branch("log_sb", &this_bin.log_sb);
      bandTree->Branch("b_hi", &this_bin.b_hi);
      bandTree->Branch("b_lo", &this_bin.b_lo);
      bandTree->Branch("channel", &this_bin.channel);

      for (auto b : bins) {
        this_bin.d = b.d;
        this_bin.b = b.b;
        this_bin.b0 = b.b0;
        this_bin.s = b.s;
        this_bin.s1 = b.s1;
        this_bin.sb = b.sb;
        this_bin.log_sb = b.log_sb;
        this_bin.b_hi = b.b_hi;
        this_bin.b_lo = b.b_lo;
        this_bin.channel = b.channel;

        bandTree->Fill();
      }

      bandTree->ResetBranchAddresses();
      bandTree->Write("", TObject::kOverwrite);
      fout.Close();
    }
  }

  else {
    TFile fileIn(bandFileName.c_str());
    TTree *theTree = nullptr;
    TTreeReader theReader("band", &fileIn);

    TTreeReaderValue<Double_t> d(theReader, "d");
    TTreeReaderValue<Double_t> b(theReader, "b");
    TTreeReaderValue<Double_t> b0(theReader, "b0");
    TTreeReaderValue<Double_t> s(theReader, "s");
    TTreeReaderValue<Double_t> s1(theReader, "s1");
    TTreeReaderValue<Double_t> sb(theReader, "sb");
    TTreeReaderValue<Double_t> log_sb(theReader, "log_sb");
    TTreeReaderValue<Double_t> b_hi(theReader, "b_hi");
    TTreeReaderValue<Double_t> b_lo(theReader, "b_lo");
    TTreeReaderValue<string> channel(theReader, "channel");

    while (theReader.Next()) {
      bin this_bin;
      this_bin.d = *d;
      this_bin.b = *b;
      this_bin.b0 = *b0;
      this_bin.s = *s;
      this_bin.s1 = *s1;
      this_bin.sb = *sb;
      this_bin.log_sb = *log_sb;
      this_bin.b_hi = *b_hi;
      this_bin.b_lo = *b_lo;
      this_bin.channel = *channel;

      // if (this_bin.channel != "HML") continue;

      bins.push_back(this_bin);
    }
  }

  // Get the total signal
  double nSigTot = 0.;
  double nTot = 0.;
  for (int i = 0; i < bins.size(); ++i) {
    double s = bins[i].s + bins[i].s1;
    nSigTot += s;
    nTot += bins[i].b + s;
  }

  LOG(logDEBUG) << "Total signal: " << nSigTot;

  vector<string> unique_channels = channels;

  sort(unique_channels.begin(), unique_channels.end());
  unique_channels.erase(unique(unique_channels.begin(), unique_channels.end()), unique_channels.end());

  // sort by log10(s/b)
  std::sort(bins.begin(), bins.end(), [](const bin &a, const bin &b) { return a.log_sb < b.log_sb; });

  // binning
  int nbins = bins.size();

  vector<double> bin_boundaries;
  bin_boundaries.push_back(bins[0].log_sb - (bins[1].log_sb - bins[0].log_sb) / 2.);

  for (int i = 0; i < nbins - 1; ++i) {
    double this_bin_center = bins[i].log_sb;
    double next_bin_center = bins[i + 1].log_sb;
    double edge = this_bin_center + (next_bin_center - this_bin_center) / 2.;
    bin_boundaries.push_back(edge);
  }

  bin_boundaries.push_back(bins[nbins - 1].log_sb + (bins[nbins - 1].log_sb - bins[nbins - 2].log_sb) / 2.);

  // fill histograms
  TH1D *h_d = new TH1D("h_d", "h_d", nbins, &bin_boundaries[0]);

  TH1D *h_b = new TH1D("h_b", "h_b", nbins, &bin_boundaries[0]);
  TH1D *h_b_hi = new TH1D("h_b_hi", "h_b_hi", nbins, &bin_boundaries[0]);
  TH1D *h_b_lo = new TH1D("h_b_lo", "h_b_lo", nbins, &bin_boundaries[0]);

  TH1D *h_b0 = new TH1D("h_b0", "h_b0", nbins, &bin_boundaries[0]);
  TH1D *h_s = new TH1D("h_s", "h_s", nbins, &bin_boundaries[0]);
  TH1D *h_s1 = new TH1D("h_s1", "h_s1", nbins, &bin_boundaries[0]);

  TH1D *h_tot = new TH1D("h_tot", "h_tot", nbins, &bin_boundaries[0]);
  TH1D *h_tot_sig = new TH1D("h_tot_sig", "h_tot_sig", nbins, &bin_boundaries[0]);
  map<string, TH1D *> component_histograms;
  map<string, TH1D *> component_histograms_significance;

  for (auto c : unique_channels) {
    string name = "h_" + c;
    string name_sig = "h_sig_" + c;
    component_histograms[c] = new TH1D(name.c_str(), name.c_str(), nbins, &bin_boundaries[0]);
    component_histograms_significance[c] = new TH1D(name_sig.c_str(), name_sig.c_str(), nbins, &bin_boundaries[0]);
  }

  for (int i = 0; i < bins.size(); ++i) {
    h_d->SetBinContent(i + 1, bins[i].d);

    h_b->SetBinContent(i + 1, bins[i].b);
    h_b_hi->SetBinContent(i + 1, bins[i].b);
    h_b_lo->SetBinContent(i + 1, bins[i].b);

    h_b0->SetBinContent(i + 1, bins[i].b0);
    h_b_hi->SetBinError(i + 1, bins[i].b_hi);
    h_b_lo->SetBinError(i + 1, bins[i].b_lo);

    h_s->SetBinContent(i + 1, bins[i].s);
    h_s1->SetBinContent(i + 1, bins[i].s1);

    LOG(logDEBUG) << bins[i].channel << " " << bins[i].s1 + bins[i].s + bins[i].b;

    double s = bins[i].s1 + bins[i].s;
    double b = bins[i].b;
    double sb = (h_b_hi->GetBinError(i + 1) + h_b_lo->GetBinError(i + 1)) / 2.;

    double z2 = (2 * ((s + b) * std::log(((s + b) * (b + sb * sb)) / (b * b + (s + b) * sb * sb)) - (b * b) / (sb * sb) * std::log(1 + (sb * sb * s) / (b * (b + sb * sb)))));

    component_histograms[bins[i].channel]->SetBinContent(i + 1, bins[i].s1 + bins[i].s + bins[i].b);
    component_histograms_significance[bins[i].channel]->SetBinContent(i + 1, z2);
    h_tot->SetBinContent(i + 1, bins[i].s1 + bins[i].s + bins[i].b);
    h_tot_sig->SetBinContent(i + 1, z2);
  }

  // rebin
  int nn = 0;

  for (auto b : bins) {
    LOG(logDEBUG) << nn - 1 - nbins;
    LOG(logDEBUG) << "bin " << nn << " log s/b=" << b.log_sb << " channel: " << b.channel << " b=" << b.b << " b_hi=" << b.b_hi << " b_lo=" << b.b_lo  << " s=" << b.s1 + b.s << " s+b=" << b.b + b.s1 + b.s << " z=" << h_tot_sig->GetBinContent(nn + 1) << " data=" << b.d;
    nn++;
  }

  double bin_width = std::log10(nTot) / nTargetBins;

  LOG(logDEBUG) << "nTot=" << nTot;
  LOG(logDEBUG) << "log10(nTot)=" << std::log10(nTot);
  LOG(logDEBUG) << "nTargetBins=" << nTargetBins;
  LOG(logDEBUG) << "bin_width=" << bin_width;

  vector<double> bin_boundaries_rebin;
  bin_boundaries_rebin.push_back(bins[nbins - 1].log_sb + (bins[nbins - 1].log_sb - bins[nbins - 2].log_sb) / 2.);

  Double_t xmin = 1.0;
  Double_t xmax = nTot;
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin) / (nTargetBins-1);
  vector<Double_t> xbins;
  xbins.push_back(xmin);

  for (Int_t i=1;i<=(nTargetBins-1);i++) {
    xbins.push_back(xmin + TMath::Power(10,logxmin+i*binwidth));
  }

  double this_yield = 0.0;
  int n = 0;
  int bcount = 0;

  for (auto b : boost::adaptors::reverse(bins)) {
    this_yield += b.b + b.s1 + b.s;

    LOG(logDEBUG) << "this_yield=" << this_yield << " target=" << xbins[n];
    bcount++;

    double edge_candidate = bin_boundaries[bin_boundaries.size() - bcount -1];

    double this_bin_wdith = fabs(edge_candidate - bin_boundaries_rebin.back());
    double threshold = fabs(bin_boundaries.front() - bin_boundaries.back()) / (1.5 * (nTargetBins-1));

    if (this_yield > xbins[n] && this_bin_wdith > threshold) {
      bin_boundaries_rebin.push_back(edge_candidate);
      n++;
    }
  }

// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 128]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 121]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 105]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 89]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 66]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 52]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 39]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 33]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 26]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 20]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 16]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 13]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 8]);
// bin_boundaries_rebin.push_back(bin_boundaries[bin_boundaries.size() - 4]);

  bin_boundaries_rebin.push_back(bins[0].log_sb - (bins[1].log_sb - bins[0].log_sb) / 2.);
  std::reverse(bin_boundaries_rebin.begin(),bin_boundaries_rebin.end());

  // bin_boundaries_rebin = bin_boundaries;

  h_d = (TH1D *)h_d->Rebin(bin_boundaries_rebin.size() - 1, "h_d_rebin", &bin_boundaries_rebin[0]);
  h_b0 = (TH1D *)h_b0->Rebin(bin_boundaries_rebin.size() - 1, "h_b0_rebin", &bin_boundaries_rebin[0]);
  h_b = (TH1D *)h_b->Rebin(bin_boundaries_rebin.size() - 1, "h_b_rebin", &bin_boundaries_rebin[0]);
  h_b_hi = (TH1D *)h_b_hi->Rebin(bin_boundaries_rebin.size() - 1, "h_b_hi_rebin", &bin_boundaries_rebin[0]);
  h_b_lo = (TH1D *)h_b_lo->Rebin(bin_boundaries_rebin.size() - 1, "h_b_lo_rebin", &bin_boundaries_rebin[0]);
  h_s = (TH1D *)h_s->Rebin(bin_boundaries_rebin.size() - 1, "h_s_rebin", &bin_boundaries_rebin[0]);
  h_s1 = (TH1D *)h_s1->Rebin(bin_boundaries_rebin.size() - 1, "h_s1_rebin", &bin_boundaries_rebin[0]);

  h_tot = (TH1D *)h_tot->Rebin(bin_boundaries_rebin.size() - 1, "h_tot_rebin", &bin_boundaries_rebin[0]);
  h_tot_sig = (TH1D *)h_tot_sig->Rebin(bin_boundaries_rebin.size() - 1, "h_tot_sig_rebin", &bin_boundaries_rebin[0]);

  for (auto c : unique_channels) {
    string name = "h_" + c + "_rebin";
    string name_sig = "h_sig_" + c + "_rebin";
    component_histograms[c] = (TH1D *)component_histograms[c]->Rebin(bin_boundaries_rebin.size() - 1, name.c_str(), &bin_boundaries_rebin[0]);
    component_histograms_significance[c] = (TH1D *)component_histograms_significance[c]->Rebin(bin_boundaries_rebin.size() - 1, name_sig.c_str(), &bin_boundaries_rebin[0]);
  }

  for (auto c : unique_channels) {
    component_histograms[c]->Divide(h_tot);
    component_histograms_significance[c]->Divide(h_tot_sig);
    // component_histograms[c]->Multiply(h_b);
  }

  LOG(logDEBUG) << "Significnace";

  for (int i = 0; i < h_tot->GetNbinsX(); ++i) {
    LOG(logDEBUG) << h_tot_sig->GetBinContent(i+1);

    double old_fraction = 0.0;
    double old_yield = 0.0;

    double old_fraction_sig = 0.0;
    double old_yield_sig = 0.0;

    bool isFirst = true;

    for (auto c : unique_channels) {
      double b_val = h_b->GetBinContent(i + 1);

      double frac = component_histograms[c]->GetBinContent(i + 1);
      double frac_sig = component_histograms_significance[c]->GetBinContent(i + 1);

      double k = 1.0;
      double cc = 0.0;

      double expr = k * (old_fraction + frac) + cc;
      double old_expr = k * (old_fraction) + cc;

      double expr_sig = k * (old_fraction_sig + frac_sig) + cc;
      double old_expr_sig = k * (old_fraction_sig) + cc;

      double l_val = expr - old_expr;
      double l_val_sig = expr_sig - old_expr_sig;

      component_histograms[c]->SetBinContent(i + 1, l_val);
      component_histograms_significance[c]->SetBinContent(i + 1, l_val_sig);

      old_fraction += frac;
      old_fraction_sig += frac_sig;

      old_yield += l_val;
      old_yield_sig += l_val_sig;
    }
  }

  TGraphAsymmErrors *backround_band = new TGraphAsymmErrors(h_b);
  for (int i = 0; i < backround_band->GetN(); ++i) {
    double x = (bin_boundaries_rebin[i] + bin_boundaries_rebin[i + 1]) / 2;
    double y = h_b->GetBinContent(i + 1);

    backround_band->SetPoint(i, x, y);

    backround_band->SetPointEXlow(i, fabs(x - bin_boundaries_rebin[i]));
    backround_band->SetPointEXhigh(i, fabs(x - bin_boundaries_rebin[i + 1]));

    backround_band->SetPointEYhigh(i, h_b_hi->GetBinError(i + 1));
    backround_band->SetPointEYlow(i, h_b_lo->GetBinError(i + 1));
  }

  h_s1->SetLineColor(myOrange);
  h_s1->SetLineWidth(0);
  h_s1->SetFillColor(myOrange);
  h_s->SetLineColor(myRed);
  h_s->SetLineWidth(0);
  h_s->SetFillColor(myRed);
  backround_band->SetLineColor(kWhite);
  backround_band->SetFillColor(kBlack);
  backround_band->SetFillStyle(3245);
  h_b0->SetLineColor(kBlack);
  h_b0->SetLineStyle(2);

  // data
  TGraphAsymmErrors *data_graph = new TGraphAsymmErrors(h_d);
  for (int p = 0; p < data_graph->GetN(); p++) {
    double x, y;
    data_graph->GetPoint(p, x, y);
    data_graph->SetPointError(p, 0, 0, fabs(calcPoissonCLLower(0.68, y) - y), fabs(calcPoissonCLUpper(0.68, y) - y));

    if (y < 1e-9) {
      data_graph->SetPoint(p, x, 1e-9);
    }
  }

  // stack
  THStack *h_stack = new THStack("stack", "stack");
  h_stack->Add(h_b);
  h_stack->Add(h_s1);
  h_stack->Add(h_s);

  // Color_t color_H4l = TColor::GetColor("#272943");
  // Color_t color_Hyy = TColor::GetColor("#D4FC76");
  // Color_t color_Hbb = TColor::GetColor("#52CCD0");
  // Color_t color_HML = TColor::GetColor("#31F1B3");

  Color_t color_Hbb = TColor::GetColor("#eff3ff");
  Color_t color_HML = TColor::GetColor("#bdd7e7");
  Color_t color_Hyy = TColor::GetColor("#6baed6");
  Color_t color_H4l = TColor::GetColor("#2171b5");

  // Color_t color_Hbb = TColor::GetColor("#edf8e9");
  // Color_t color_HML = TColor::GetColor("#bae4b3");
  // Color_t color_Hyy = TColor::GetColor("#74c476");
  // Color_t color_H4l = TColor::GetColor("#238b45");

  component_histograms["H4l"]->SetLineColor(color_H4l);
  component_histograms["H4l"]->SetLineWidth(0.0);
  component_histograms["H4l"]->SetFillColor(color_H4l);

  component_histograms["Hyy"]->SetLineColor(color_Hyy);
  component_histograms["Hyy"]->SetLineWidth(0.0);
  component_histograms["Hyy"]->SetFillColor(color_Hyy);

  component_histograms["Hbb"]->SetLineColor(color_Hbb);
  component_histograms["Hbb"]->SetLineWidth(0.0);
  component_histograms["Hbb"]->SetFillColor(color_Hbb);

  component_histograms["HML"]->SetLineColor(color_HML);
  component_histograms["HML"]->SetLineWidth(0.0);
  component_histograms["HML"]->SetFillColor(color_HML);

  component_histograms_significance["H4l"]->SetLineColor(color_H4l);
  component_histograms_significance["H4l"]->SetLineWidth(0.0);
  component_histograms_significance["H4l"]->SetFillColor(color_H4l);

  component_histograms_significance["Hyy"]->SetLineColor(color_Hyy);
  component_histograms_significance["Hyy"]->SetLineWidth(0.0);
  component_histograms_significance["Hyy"]->SetFillColor(color_Hyy);

  component_histograms_significance["Hbb"]->SetLineColor(color_Hbb);
  component_histograms_significance["Hbb"]->SetLineWidth(0.0);
  component_histograms_significance["Hbb"]->SetFillColor(color_Hbb);

  component_histograms_significance["HML"]->SetLineColor(color_HML);
  component_histograms_significance["HML"]->SetLineWidth(0.0);
  component_histograms_significance["HML"]->SetFillColor(color_HML);

  THStack *h_stack_components = new THStack("component_stack", "component_stack");
  // h_stack_components->Add(h_tot);

  THStack *h_stack_components_sig = new THStack("component_stack_sig", "component_stack_sig");
  // h_stack_components_sig->Add(h_tot_sig);

  // for (auto c : unique_channels) {
  //   string name = "h_" + c;
  //   h_stack_components->Add(component_histograms[c]);
  //   h_stack_components_sig->Add(component_histograms_significance[c]);
  // }

  h_stack_components->Add(component_histograms["Hbb"]);
  h_stack_components->Add(component_histograms["HML"]);
  h_stack_components->Add(component_histograms["Hyy"]);
  h_stack_components->Add(component_histograms["H4l"]);

  // histograms for bottom panel
  TH1D *h_d2 = new TH1D("h_d2", "h_d2", bin_boundaries_rebin.size() - 1, &bin_boundaries_rebin[0]);
  TH1D *h_s2 = new TH1D("h_s2", "h_s2", bin_boundaries_rebin.size() - 1, &bin_boundaries_rebin[0]);
  TH1D *h_s12 = new TH1D("h_s12", "h_s12", bin_boundaries_rebin.size() - 1, &bin_boundaries_rebin[0]);
  TH1D *h_b02 = new TH1D("h_b02", "h_b02", bin_boundaries_rebin.size() - 1, &bin_boundaries_rebin[0]);
  TH1D *h_b2 = new TH1D("h_b2", "h_b2", bin_boundaries_rebin.size() - 1, &bin_boundaries_rebin[0]);

  vector<double> values_pad2;

  for (int i = 0; i < backround_band->GetN(); ++i) {

    double e = (h_b_hi->GetBinError(i + 1) + h_b_lo->GetBinError(i + 1)) / 2.;

    // double d2 = (h_d->GetBinContent(i + 1) - h_b->GetBinContent(i + 1));
    // double s2 = (h_s->GetBinContent(i + 1) + h_s1->GetBinContent(i + 1));
    // double s12 = h_s1->GetBinContent(i + 1);
    // double b02 = (h_b0->GetBinContent(i + 1) - h_b->GetBinContent(i + 1));
    // double b2 = (h_b->GetBinContent(i + 1) - h_b->GetBinContent(i + 1));

    double d2 = (h_d->GetBinContent(i + 1)) / h_b->GetBinContent(i+1);
    double s2 = (h_s->GetBinContent(i + 1) + h_s1->GetBinContent(i + 1) + h_b->GetBinContent(i+1)) / h_b->GetBinContent(i+1);
    double s12 = (h_s1->GetBinContent(i + 1) + h_b->GetBinContent(i+1)) / h_b->GetBinContent(i+1);
    double b02 = (h_b0->GetBinContent(i + 1)) / h_b->GetBinContent(i+1);
    double b2 = (h_b->GetBinContent(i + 1)) / h_b->GetBinContent(i+1);

    // double d2 = (h_d->GetBinContent(i + 1) - h_b->GetBinContent(i + 1)) / e;
    // double s2 = (h_s->GetBinContent(i + 1) + h_s1->GetBinContent(i + 1)) / e;
    // double s12 = h_s1->GetBinContent(i + 1) / e;
    // double b02 = (h_b0->GetBinContent(i + 1) - h_b->GetBinContent(i + 1)) / e;
    // double b2 = 0;

    // double d2 = (h_d->GetBinContent(i + 1) - h_b->GetBinContent(i + 1)) / sqrt(h_b->GetBinContent(i+1) + e * e);
    // double s2 = (h_s->GetBinContent(i + 1) + h_s1->GetBinContent(i + 1)) / sqrt(h_b->GetBinContent(i+1) + e * e);
    // double s12 = h_s1->GetBinContent(i + 1) / sqrt(h_b->GetBinContent(i+1) + e * e);
    // double b02 = (h_b0->GetBinContent(i + 1) - h_b->GetBinContent(i + 1)) / sqrt(h_b->GetBinContent(i+1) + e * e);
    // double b2 = 0;

    // double d2 = (h_d->GetBinContent(i + 1) - h_b->GetBinContent(i + 1)) / sqrt(h_b->GetBinContent(i + 1));
    // double s2 = (h_s->GetBinContent(i + 1) + h_s1->GetBinContent(i + 1)) / sqrt(h_b->GetBinContent(i + 1));
    // double s12 = h_s1->GetBinContent(i + 1) / sqrt(h_b->GetBinContent(i + 1));
    // double b02 = (h_b0->GetBinContent(i + 1) - h_b->GetBinContent(i + 1)) / sqrt(h_b->GetBinContent(i + 1));
    // double b2 = 0;

    h_d2->SetBinContent(i + 1, d2);
    h_s2->SetBinContent(i + 1, s2);
    h_s12->SetBinContent(i + 1, s12);
    h_b02->SetBinContent(i + 1, b02);
    h_b2->SetBinContent(i + 1, b2);

    values_pad2.push_back(d2);
    values_pad2.push_back(s2);
    values_pad2.push_back(s12);
    values_pad2.push_back(b02);
    values_pad2.push_back(b2);

  }

  // data for bottom panel
  TGraphAsymmErrors *tmp_data_graph2 = new TGraphAsymmErrors(h_d);
  TGraphAsymmErrors *data_graph2 = new TGraphAsymmErrors(h_d2);
  TGraphAsymmErrors *b2_graph2= new TGraphAsymmErrors(h_b2);

  for (int p = 0; p < data_graph2->GetN(); p++) {
    double x, y, bx, by, xtmp, ytmp;
    tmp_data_graph2->GetPoint(p, xtmp, ytmp);
    data_graph2->GetPoint(p, x, y);
    b2_graph2->GetPoint(p, bx, by);

    double lo = fabs(calcPoissonCLLower(0.68, ytmp) - ytmp);
    double hi = fabs(calcPoissonCLUpper(0.68, ytmp) - ytmp);

    double e_lo = h_b_lo->GetBinError(p + 1);
    double e_hi = h_b_hi->GetBinError(p + 1);
    double e = (e_lo + e_hi) / 2.;


    // double tip_hi = (h_d->GetBinContent(p + 1) + hi - h_b->GetBinContent(p + 1));
    // double tip_lo = (h_d->GetBinContent(p + 1) - lo - h_b->GetBinContent(p + 1));
    // double tipe_hi = (h_b->GetBinContent(p + 1) + e_hi) - h_b->GetBinContent(p+1);
    // double tipe_lo = (h_b->GetBinContent(p + 1) - e_lo) - h_b->GetBinContent(p+1);

    double tip_hi = (h_d->GetBinContent(p + 1) + hi) / h_b->GetBinContent(p+1);
    double tip_lo = (h_d->GetBinContent(p + 1) - lo) / h_b->GetBinContent(p+1);
    double tipe_hi = (h_b->GetBinContent(p + 1) + e_hi) / h_b->GetBinContent(p+1);
    double tipe_lo = (h_b->GetBinContent(p + 1) - e_lo) / h_b->GetBinContent(p+1);

    // double tip_hi = (h_d->GetBinContent(p + 1) + hi - h_b->GetBinContent(p + 1)) / e;
    // double tip_lo = (h_d->GetBinContent(p + 1) - lo - h_b->GetBinContent(p + 1)) / e;
    // double tipe_hi = e / e;
    // double tipe_lo = e / e;

    // double tip_hi = (h_d->GetBinContent(p + 1) + hi - h_b->GetBinContent(p + 1)) / sqrt(h_b->GetBinContent(p+1) + e * e);
    // double tip_lo = (h_d->GetBinContent(p + 1) - lo - h_b->GetBinContent(p + 1)) / sqrt(h_b->GetBinContent(p+1) + e * e);

    // double tip_hi = (h_d->GetBinContent(p + 1) + hi - h_b->GetBinContent(p + 1)) / sqrt(h_b->GetBinContent(p+1));
    // double tip_lo = (h_d->GetBinContent(p + 1) - lo - h_b->GetBinContent(p + 1)) / sqrt(h_b->GetBinContent(p+1));
    // double tipe_hi = (e) / sqrt(h_b->GetBinContent(p+1));
    // double tipe_lo = (e) / sqrt(h_b->GetBinContent(p+1));

    data_graph2->SetPointError(p, 0, 0, fabs(tip_lo - y), fabs(tip_hi - y));
    values_pad2.push_back(tip_lo);
    values_pad2.push_back(tip_hi);

    // b2_graph2->SetPointError(p, 0, 0, fabs(tipe_lo - by), fabs(tipe_hi - by));
    b2_graph2->SetPointError(p, 0, 0, fabs(tipe_lo - by), fabs(tipe_hi - by));
    values_pad2.push_back(-tipe_lo);
    values_pad2.push_back(+tipe_hi);
  }

  h_s12->SetLineColor(myOrange);
  h_s12->SetLineStyle(7);
  h_s2->SetLineColor(myRed);
  h_b02->SetLineColor(kBlack);
  h_b02->SetLineStyle(2);

  for (int i = 0; i < b2_graph2->GetN(); ++i) {
    double x = (bin_boundaries_rebin[i] + bin_boundaries_rebin[i + 1]) / 2;
    double y = h_b->GetBinContent(i + 1);

    double bx, by;
    b2_graph2->GetPoint(i, bx, by);

    b2_graph2->SetPoint(i, x, by);
    b2_graph2->SetPointEXlow(i, fabs(x - bin_boundaries_rebin[i]));
    b2_graph2->SetPointEXhigh(i, fabs(x - bin_boundaries_rebin[i + 1]));

    // b2_graph2->SetPointEYhigh(i, h_b_hi->GetBinError(i + 1));
    // b2_graph2->SetPointEYlow(i, h_b_lo->GetBinError(i + 1));
  }

  b2_graph2->SetLineColor(kBlack);
  b2_graph2->SetFillColor(kBlack);
  b2_graph2->SetFillStyle(3245);

  // plotting
  TCanvas *c;
  TPad *pad1;
  TPad *pad2;
  TPad *pad3;

  if (!showOrigin) {
    c = new TCanvas("c", "canvas", 800, 800);

    pad1 = new TPad("pad1", "pad1", 0, 0.325, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetLogy(1);
    pad1->Draw();

    c->cd();
    pad2 = new TPad("pad2", "pad2", 0, 0.025, 1, 0.325);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.325);
    pad2->Draw();
  } else {
    c = new TCanvas("c", "canvas", 800, 540.+150.+250.);

    pad1 = new TPad("pad1", "pad1", 0, 1 - 540./(540.+150.+250.), 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetLogy(1);
    pad1->Draw();

    c->cd();
    pad3 = new TPad("pad3", "pad3", 0, 1 - 540./(540.+150.+250.) - 150./(540.+150.+250.), 1, 1 - 540./(540.+150.+250.));
    pad3->SetTopMargin(0);
    pad3->SetBottomMargin(0);
    pad3->Draw();

    c->cd();
    pad2 = new TPad("pad2", "pad2", 0, 0.025, 1, 1 - 540./(540.+150.+250.) - 150./(540.+150.+250.));
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(320./(540.+150.+250.));
    pad2->Draw();
  }

  // boundaries
  pad1->cd();
  vector<double> edges1 = {0.5, 1e8};
  TH2D *boundary1 = new TH2D("boundary1", "", nbins, &bin_boundaries[0], 1, &edges1[0]);
  double sizeTick1 = 20.0;
  double padW1 = pad1->GetWw() * pad1->GetAbsWNDC();
  double padH1 = pad1->GetWh() * pad1->GetAbsHNDC();
  double tickScaleX1 = (pad1->GetUxmax() - pad1->GetUxmin()) / (pad1->GetX2() - pad1->GetX1()) * padH1;
  double tickScaleY1 = (pad1->GetUymax() - pad1->GetUymin()) / (pad1->GetY2() - pad1->GetY1()) * padW1;
  boundary1->GetXaxis()->SetTickLength(sizeTick1 / tickScaleX1);
  boundary1->GetYaxis()->SetTickLength(sizeTick1 / tickScaleY1);
  boundary1->GetXaxis()->SetLabelSize(28);
  boundary1->GetYaxis()->SetLabelSize(28);
  boundary1->GetXaxis()->SetTitleSize(28);
  boundary1->GetYaxis()->SetTitleSize(28);
  boundary1->SetTitle("");
  boundary1->GetYaxis()->SetTitle("Events/bin");

  pad2->cd();
  double min2 = *std::min_element(values_pad2.begin(), values_pad2.end());
  double max2 = *std::max_element(values_pad2.begin(), values_pad2.end());
  double delta2 = max2 - min2;

  // vector<double> edges2 = {min2 - delta2 * 0.175, max2 + delta2 * 0.175};
  // vector<double> edges2 = {min2 - delta2 * 0.175, max2 + delta2 * 0.205};
  // vector<double> edges2 = {-1.5, 9};
  vector<double> edges2 = {0.5, 4.75};
  TH2D *boundary2 = new TH2D("boundary2", "", nbins, &bin_boundaries[0], 1, &edges2[0]);
  double sizeTick2 = 20.0;
  double padW2 = pad2->GetWw() * pad2->GetAbsWNDC();
  double padH2 = pad2->GetWh() * pad2->GetAbsHNDC();
  double tickScaleX2 = (pad2->GetUxmax() - pad2->GetUxmin()) / (pad2->GetX2() - pad2->GetX1()) * padH2;
  double tickScaleY2 = (pad2->GetUymax() - pad2->GetUymin()) / (pad2->GetY2() - pad2->GetY1()) * padW2;
  boundary2->GetXaxis()->SetTickLength(sizeTick2 / tickScaleX2);
  boundary2->GetYaxis()->SetTickLength(sizeTick2 / tickScaleY2);
  boundary2->GetXaxis()->SetLabelSize(28);
  boundary2->GetYaxis()->SetLabelSize(28);
  boundary2->GetXaxis()->SetTitleSize(28);
  boundary2->GetYaxis()->SetTitleSize(28);
  boundary2->SetTitle("");
  boundary2->GetXaxis()->SetTitle("log_{10}(S/B)");
  // boundary2->GetYaxis()->SetTitle("Data - Bkgd.");
  boundary2->GetYaxis()->SetTitle("Data/Bkgd.");
  // boundary2->GetYaxis()->SetTitle("S/#sigma_{B}");
  // boundary2->GetYaxis()->SetTitle("S/#sqrt{B}");
  // boundary2->GetYaxis()->SetTitle("Pull");
  // boundary2->GetYaxis()->SetTitle("S/#sqrt{B + #sigma^{2}_{B}}");
  boundary2->GetXaxis()->SetTitleOffset(4);
  boundary2->GetYaxis()->SetNdivisions(010);

  TH2D *boundary3;
  
  if (showOrigin) {
    vector<double> edges3 = {0.0, 1.0};
    boundary3 = new TH2D("boundary3", "", nbins, &bin_boundaries[0], 1, &edges3[0]);
    double sizeTick3 = 20.0;
    double padW3 = pad3->GetWw() * pad3->GetAbsWNDC();
    double padH3 = pad3->GetWh() * pad3->GetAbsHNDC();
    double tickScaleX3 = (pad3->GetUxmax() - pad3->GetUxmin()) / (pad3->GetX2() - pad3->GetX1()) * padH3;
    double tickScaleY3 = (pad3->GetUymax() - pad3->GetUymin()) / (pad3->GetY2() - pad3->GetY1()) * padW3;
    boundary3->GetXaxis()->SetTickLength(sizeTick3 / tickScaleX3);
    boundary3->GetYaxis()->SetTickLength(sizeTick3 / tickScaleY3);
    boundary3->GetXaxis()->SetLabelSize(28);
    boundary3->GetYaxis()->SetLabelSize(28);
    boundary3->GetXaxis()->SetTitleSize(28);
    boundary3->GetYaxis()->SetTitleSize(28);
    boundary3->SetTitle("");
    boundary3->GetYaxis()->SetTitle("S+B frac.");
    boundary3->GetXaxis()->SetTitleOffset(4);
    boundary3->GetYaxis()->SetNdivisions(010);
    boundary1->GetYaxis()->SetTitleOffset(1.75);
    boundary2->GetYaxis()->SetTitleOffset(1.75);
    boundary2->GetXaxis()->SetTitleOffset(5);
    boundary3->GetYaxis()->SetTitleOffset(1.75);
  }

  // TLine *line = new TLine(bin_boundaries.front(), 0.0, bin_boundaries.back(), 0.0);
  TLine *line = new TLine(bin_boundaries.front(), 1.0, bin_boundaries.back(), 1.0);
  // line->SetLineColor(kGray + 1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(1);
  // line->SetLineStyle(kDashed);

  // draw histograms
  pad1->cd();
  boundary1->Draw("hist");
  h_stack->Draw("same");
  h_b->Draw("same");
  h_b0->Draw("same");
  data_graph->DrawClone("EP same");
  // backround_band->Draw("2 SAME");
  pad1->RedrawAxis();

  if (showOrigin) {
    pad3->cd();
    boundary3->Draw("hist");

    h_stack_components->Draw("same");
    // h_stack_components_sig->Draw("same");

    pad3->RedrawAxis();
  }

  // draw ratio
  pad2->cd();
  boundary2->Draw("hist");
  line->Draw("same");
  // b2_graph2->Draw("2 SAME");
  // h_b02->Draw("same");
  h_s12->Draw("same");
  h_s2->Draw("same");
  data_graph2->DrawClone("EP same");
  pad2->RedrawAxis();

  // legend
  pad1->cd();
  TLegend *legend = new TLegend(0.7, 0.55, pad1->GetX2() - 1.5 * pad1->GetRightMargin(), pad1->GetY2() - 2 * pad1->GetTopMargin());

  legend->AddEntry(data_graph, "Data", "EP");
  legend->AddEntry(h_s, "t#bar{t}H (#mu=1.32)", "F");
  legend->AddEntry(h_s1, "t#bar{t}H (#mu=1)", "F");
  legend->AddEntry(h_b, "Background", "F");
  // legend->AddEntry(backround_band, "Bkgd. Unc.", "F");
  legend->AddEntry(h_b0, "Bkgd. (#mu=0)", "L");

  TLegend *legend2 = new TLegend(pad1->GetX1() + 1.25 * pad1->GetLeftMargin(), pad1->GetY1() + 2 * pad1->GetBottomMargin()+0.025, 0.45, 0.3);

  legend2->AddEntry(component_histograms["H4l"], "t#bar{t}H (ZZ)", "F");
  legend2->AddEntry(component_histograms["Hyy"], "t#bar{t}H (#gamma#gamma)", "F");
  legend2->AddEntry(component_histograms["HML"], "t#bar{t}H (multilepton)", "F");
  legend2->AddEntry(component_histograms["Hbb"], "t#bar{t}H (b#bar{b})", "F");

  legend->SetFillStyle(0);
  legend->SetBorderSize(1);
  legend->SetLineWidth(0);
  legend->SetTextFont(43);
  legend->SetTextSize(24);
  legend->SetNColumns(1);
  legend->Draw();

  if (showOrigin) {
  legend2->SetFillStyle(0);
  legend2->SetBorderSize(1);
  legend2->SetLineWidth(0);
  legend2->SetTextFont(43);
  legend2->SetTextSize(24);
  legend2->SetNColumns(1);
  legend2->Draw();
}

  // decorations
  pad1->cd();

  TLatex tex;
  tex.SetTextFont(72);
  tex.SetTextSize(1.4 * tex.GetTextSize());
  tex.SetTextColor(kBlack);
  tex.SetNDC();
  tex.DrawLatex(pad1->GetX1() + 1.25 * pad1->GetLeftMargin(), pad1->GetY2() - 3.0 * pad1->GetTopMargin(), "ATLAS");

  TLatex tex_label;
  tex_label.SetTextFont(42);
  tex_label.SetTextSize(1.4 * tex_label.GetTextSize());
  tex_label.SetTextColor(1);
  tex_label.SetNDC();
  // tex_label.DrawLatex(pad1->GetX1() + 1.25 * pad1->GetLeftMargin() + 0.175, pad1->GetY2() - 3.0 * pad1->GetTopMargin(), "Internal");

  TLatex l;
  l.SetTextFont(42);
  l.SetTextSize(1.15 * l.GetTextSize());
  l.SetNDC();

  string energy = "13";
  string lumi = "36.1 - 79.8";

  // string lumi_label = "#sqrt{s}=" + energy + " TeV, " +  boost::str(boost::format("%.1f") % lumi) + " fb^{-1}";
  string lumi_label = "#sqrt{s}=" + energy + " TeV, " + lumi + " fb^{-1}";
  l.DrawLatex(pad1->GetX1() + 1.25 * pad1->GetLeftMargin(), pad1->GetY2() - 4.5 * pad1->GetTopMargin(), lumi_label.c_str());

  // cover up

  if (showOrigin) {
    c->cd();
    TPad *b1 = new TPad("b1", "b1", 0.12, 0.405, 0.155, 0.43);
    b1->SetFillColor(kWhite);
    b1->SetBorderMode(0);
    b1->Draw();

    c->cd();
    TPad *b2 = new TPad("b2", "b2", 0.12, 0.405 - 0.14, 0.155, 0.43 - 0.15);
    b2->SetFillColor(kWhite);
    b2->SetBorderMode(0);
    b2->Draw();
  }

  // save plot
  c->SaveAs("log_sb.pdf");
  c->SaveAs("log_sb.C");
  c->SaveAs("log_sb.root");
  c->SaveAs("log_sb.eps");
  c->SaveAs("log_sb.png");

  PrintResourcesUsed(thistime);
}

// _____________________________________________________________________________
vector<bin> getBinnedValues(RooAbsPdf *pdfi, RooDataSet *datai, RooRealVar *x, RooRealVar *poi, bool bandFromToys, RooFitResult *fr, RooFitResult *fr_bkg, double muhat, RooArgSet *nuis, double edgeLo, double edgeHi) {
  vector<bin> bins;

  RooArgSet *vars = pdfi->getVariables();
  RooStats::RemoveConstantParameters(vars);

  RooArgSet ucmles = fr->floatParsFinal();
  RooArgSet bkg = fr_bkg->floatParsFinal();

  bool isBinned = true;
  bool mergeBins = false;

  double oldEdgeLo = x->getMin();
  double oldEdgeHi = x->getMax();
  double oldRange = oldEdgeHi - oldEdgeLo;

  if (edgeLo != -999 && edgeHi != -999) {
    // x->setMin(edgeLo);
    // x->setMax(edgeHi);
    // x->setRange(edgeLo, edgeHi);
    // x->setRange( "QTRange", edgeLo, edgeHi );
    // x->setBins(550000);
    isBinned = true;
    mergeBins = true;
  }

  // data
  TH1 *h_data = datai->createHistogram(x->GetName());

  if (edgeLo != -999 && edgeHi != -999) {
    h_data = new TH1D("h_data", "h_data", 1, 0, 1);
    string expression = string(x->GetName()) + ">" + to_string(edgeLo) + "&&" + string(x->GetName()) + "<" + to_string(edgeHi);
    double d = datai->sumEntries(expression.c_str());
    h_data->Fill(0.5, d);
  }

  cout << "---- data" << endl;
  h_data->Print();
  h_data->Print("v");

  // background only, incl. uncertainty
  for (RooLinkedListIter it = ucmles.iterator(); TObject *o = it.Next();) {
    RooRealVar *v = dynamic_cast<RooRealVar *>(o);
    RooRealVar *pdfvar = (RooRealVar *)vars->find(v->GetName());
    if (pdfvar) {
      pdfvar->setVal(v->getVal());
    }
  }

  poi->setVal(0.0);
  poi->setConstant(1);

  unique_ptr<RooDataSet> asimovData_background((RooDataSet *)AsymptoticCalculator::GenerateAsimovData(*pdfi, RooArgSet(*x)));

  if (mergeBins) {
    string expression = string(x->GetName()) + ">" + to_string(edgeLo) + "&&" + string(x->GetName()) + "<" + to_string(edgeHi);
    asimovData_background.reset((RooDataSet *)asimovData_background->reduce(expression.c_str()));
  }

  RooCurve *cenCurve_background = convertAsimovToCurve(asimovData_background.get(), x->GetName(), isBinned, mergeBins);

  // cenCurve_background->Print();
  // cenCurve_background->Print("v");

  vector<vector<double>> h_background_toys;
  vector<double> h_background;

  if (bandFromToys) {
    RooCurve *bandFull_background = createBand(cenCurve_background, pdfi, x, fr, nuis, 1, isBinned, mergeBins, edgeLo, edgeHi);
    h_background_toys = makeHistBand("hist", cenCurve_background, bandFull_background);
  } else {
    h_background = makeHist("hist", cenCurve_background);
  }

  // background only fit

  for (RooLinkedListIter it = bkg.iterator(); TObject *o = it.Next();) {
    RooRealVar *v = dynamic_cast<RooRealVar *>(o);
    RooRealVar *pdfvar = (RooRealVar *)vars->find(v->GetName());
    if (pdfvar) {
      pdfvar->setVal(v->getVal());
    }
  }

  poi->setVal(0.0);
  poi->setConstant(1);

  unique_ptr<RooDataSet> asimovData_background_only((RooDataSet *)AsymptoticCalculator::GenerateAsimovData(*pdfi, RooArgSet(*x)));

  if (mergeBins) {
    string expression = string(x->GetName()) + ">" + to_string(edgeLo) + "&&" + string(x->GetName()) + "<" + to_string(edgeHi);
    asimovData_background_only.reset((RooDataSet *)asimovData_background_only->reduce(expression.c_str()));
  }

  RooCurve *cenCurve_background_only = convertAsimovToCurve(asimovData_background_only.get(), x->GetName(), isBinned, mergeBins);
  vector<double> h_background_only = makeHist("hist", cenCurve_background_only);

  // signal (muhat) + background

  for (RooLinkedListIter it = ucmles.iterator(); TObject *o = it.Next();) {
    RooRealVar *v = dynamic_cast<RooRealVar *>(o);
    RooRealVar *pdfvar = (RooRealVar *)vars->find(v->GetName());
    if (pdfvar) {
      pdfvar->setVal(v->getVal());
    }
  }

  poi->setVal(muhat);

  unique_ptr<RooDataSet> asimovData_splusb((RooDataSet *)AsymptoticCalculator::GenerateAsimovData(*pdfi, RooArgSet(*x)));

  if (mergeBins) {
    string expression = string(x->GetName()) + ">" + to_string(edgeLo) + "&&" + string(x->GetName()) + "<" + to_string(edgeHi);
    asimovData_splusb.reset((RooDataSet *)asimovData_splusb->reduce(expression.c_str()));
  }

  RooCurve *cenCurve_splusb = convertAsimovToCurve(asimovData_splusb.get(), x->GetName(), isBinned, mergeBins);
  vector<double> h_splusb = makeHist("hist", cenCurve_splusb);

  // signal (nominal) + background

  for (RooLinkedListIter it = ucmles.iterator(); TObject *o = it.Next();) {
    RooRealVar *v = dynamic_cast<RooRealVar *>(o);
    RooRealVar *pdfvar = (RooRealVar *)vars->find(v->GetName());
    if (pdfvar) {
      pdfvar->setVal(v->getVal());
    }
  }

  poi->setVal(1.0);

  unique_ptr<RooDataSet> asimovData_s1plusb((RooDataSet *)AsymptoticCalculator::GenerateAsimovData(*pdfi, RooArgSet(*x)));

  if (mergeBins) {
    string expression = string(x->GetName()) + ">" + to_string(edgeLo) + "&&" + string(x->GetName()) + "<" + to_string(edgeHi);
    asimovData_s1plusb.reset((RooDataSet *)asimovData_s1plusb->reduce(expression.c_str()));
  }

  RooCurve *cenCurve_s1plusb = convertAsimovToCurve(asimovData_s1plusb.get(), x->GetName(), isBinned, mergeBins);
  vector<double> h_s1plusb = makeHist("hist", cenCurve_s1plusb);

  // bring us back

  for (RooLinkedListIter it = ucmles.iterator(); TObject *o = it.Next();) {
    RooRealVar *v = dynamic_cast<RooRealVar *>(o);
    RooRealVar *pdfvar = (RooRealVar *)vars->find(v->GetName());
    if (pdfvar) {
      pdfvar->setVal(v->getVal());
    }
  }

  poi->setVal(muhat);
  poi->setConstant(0);

  // collect values
  for (int i = 1; i < h_data->GetNbinsX() + 1; ++i) {
    double y_d = h_data->GetBinContent(i);
    double y_sb = h_splusb[i - 1];
    double y_s1b = h_s1plusb[i - 1];

    double y_b = 0.;
    double y_b_hi = 0.;
    double y_b_lo = 0.;

    if (bandFromToys) {
      y_b = h_background_toys[i - 1][0];
      y_b_hi = h_background_toys[i - 1][1];
      y_b_lo = h_background_toys[i - 1][2];
    } else {
      y_b = h_background[i - 1];
      y_b_hi = 1.2 * y_b;
      y_b_lo = 0.8 * y_b;
    }

    double y_b0 = h_background_only[i - 1];
    double y_s1 = y_s1b - y_b;
    double y_s = y_sb - y_b - y_s1;

    cout << y_d << " " << y_b << " " << y_s + y_s1 << endl;

    bin this_bin;
    this_bin.d = y_d;
    this_bin.b = y_b;
    this_bin.b0 = y_b0;
    this_bin.s = y_s;
    this_bin.s1 = y_s1;
    this_bin.sb = (y_sb - y_b) / y_b;
    this_bin.log_sb = log10(this_bin.sb);
    this_bin.b_hi = fabs(y_b - y_b_hi);
    this_bin.b_lo = fabs(y_b - y_b_lo);

    bins.push_back(this_bin);
  }

  if (edgeLo != -999 && edgeHi != -999) {
    x->setRange(oldEdgeLo, oldEdgeHi);
    x->setRange("QTRange", oldEdgeLo, oldEdgeHi);
  }

  return bins;
}

// _____________________________________________________________________________
void copyDataSet(RooDataSet *datadebit, RooRealVar *obsdebit, RooDataSet *datacredit, RooRealVar *obscredit, RooRealVar *wt, double scale, bool verbose) {
  RooRealVar *obs = (RooRealVar *)datadebit->get()->find(obsdebit->GetName());
  if (verbose)
    cout << datadebit->numEntries() << endl;
  double xmin = obscredit->getMin();
  double xmax = obscredit->getMax();
  // cout<<obs->GetName()<<" "<<obs->GetTitle()<<endl;getchar();
  for (int i = 0; i < datadebit->numEntries(); i++) {
    datadebit->get(i);
    if (verbose)
      cout << obs->getVal() << " " << datadebit->weight() << endl;
    double weight = datadebit->weight();
    weight *= scale;
    double value = obs->getVal();
    if (value < xmin || value > xmax) {
      cerr << "Warning: entry " << i << " with value " << value << " is out of the observable range [" << xmin << "," << xmax << "]!" << endl;
      continue;
    }
    obscredit->setVal(obs->getVal());
    wt->setVal(datadebit->weight());
    datacredit->add(RooArgSet(*obscredit, *wt), weight);
  }
  if (verbose)
    datacredit->Print("v");
}

// _____________________________________________________________________________
RooCurve *convertAsimovToCurve(RooDataSet *asimov, TString observableName, bool isHist, bool mergeBins, double epsilon) {
  RooCurve *curve = new RooCurve();
  RooArgSet *obs = const_cast<RooArgSet *>(asimov->get());
  RooRealVar *xdata = dynamic_cast<RooRealVar *>(obs->find(observableName));
  double binw = 0;
  if (isHist)
    binw = (xdata->getMax() - xdata->getMin()) / double(xdata->numBins());

  if (mergeBins) {
    curve->addPoint(0.5, asimov->sumEntries());
    return curve;
  }

  for (int i = 0; i < asimov->numEntries(); i++) {
    asimov->get(i);
    double mass = xdata->getVal();
    double weight = asimov->weight();

    if (isHist) {
      curve->addPoint(mass - 0.5 * binw + epsilon, weight);
      curve->addPoint(mass + 0.5 * binw - epsilon, weight);
      // cout<<mass<<" "
      //    <<mass-0.5*binw<<" "
      //    <<mass+0.5*binw<<" "
      //    <<weight<<endl; getchar();
    } else {
      curve->addPoint(mass, weight);
    }
  }
  return curve;
}

// _____________________________________________________________________________
void releaseTheGhost(RooDataSet *obsdata, RooRealVar *x, RooRealVar *w, double ghostwt) {
  TH1D h_data("h_data", "", x->numBins(), x->getMin(), x->getMax());
  RooArgSet *obs = const_cast<RooArgSet *>(obsdata->get());
  RooRealVar *xdata = dynamic_cast<RooRealVar *>(obs->find(x->GetName()));
  int nevt1 = obsdata->numEntries();
  for (int i = 0; i < nevt1; i++) {
    obsdata->get(i);
    h_data.Fill(xdata->getVal(), obsdata->weight());
  }
  // Now, reset!
  obsdata->reset();
  for (int ibin = 1; ibin <= x->numBins(); ibin++) {
    if (h_data.GetBinContent(ibin) == 0) {
      cout << "Bin " << ibin << " has 0 content. Inject ghost event..." << endl;
      x->setVal(h_data.GetBinCenter(ibin));
      w->setVal(ghostwt);
      obsdata->add(RooArgSet(*x, *w), ghostwt);
    } else {
      x->setVal(h_data.GetBinCenter(ibin));
      w->setVal(h_data.GetBinContent(ibin));
      obsdata->add(RooArgSet(*x, *w), h_data.GetBinContent(ibin));
    }
  }
}

// _____________________________________________________________________________
RooCurve *createBand(RooCurve *cenCurve, RooAbsPdf *pdf, RooRealVar *x, RooFitResult *res, RooArgSet *params, double Z, bool isHist, bool mergeBins, double edgeLo, double edgeHi) {
  // Generate 100 random parameter points distributed according to fit result covariance matrix

  unique_ptr<RooAbsPdf> cloneFunc((RooAbsPdf *)pdf->cloneTree());

  // cloneFunc->Print();

  RooArgSet *cloneParams = cloneFunc->getObservables(res->floatParsFinal());

  // RooArgSet *errorParams = (RooArgSet *)cloneParams->selectCommon(*params);



  RooArgSet *errorParams = new RooArgSet();
  TIterator *iter = cloneParams->createIterator();
  RooRealVar *arg;
  double val = cloneFunc->expectedEvents(RooArgSet(*x));
  cout << val << endl;
  while ((arg = (RooRealVar *)iter->Next())) {
    double aval = arg->getVal();
    double aerr = ((RooRealVar*)res->floatParsFinal().find(arg->GetName()))->getError();
    arg->setVal(aval + aerr);
    cout << arg->GetName() << " " << aval << endl;
    double newval = cloneFunc->expectedEvents(RooArgSet(*x));
    cout << newval << endl;
    arg->setVal(aval);
    TString name = arg->GetName();
    if (fabs(newval - val) > 1e-3 && !arg->isConstant()) {
      errorParams->add(*arg);
    }
  }
  delete iter;

  errorParams->Print("v");


  unique_ptr<RooAbsPdf> paramPdf(createHessePdf(res, *errorParams));

  unique_ptr<RooDataSet> asimovData;

  Int_t n = Int_t(100. / TMath::Erfc(Z / sqrt(2.))) * 100;
  // Int_t n = Int_t(100./TMath::Erfc(Z/sqrt(2.)));
  if (n < 100)
    n = 100;

  // Generate variation curves with above set of parameter values
  RooRandom::randomGenerator()->SetSeed(0);
  RooDataSet *d = paramPdf->generate(*errorParams, n);
  vector<RooCurve *> cvec;
  for (int i = 0; i < d->numEntries(); i++) {

    *cloneParams = (*d->get(i));
    unique_ptr<RooDataSet> asimovData_tmp((RooDataSet *)AsymptoticCalculator::GenerateAsimovData(*cloneFunc, RooArgSet(*x)));
    RooRealVar w("w", "w", 1);
    asimovData.reset(new RooDataSet("asimovData_tmp", "asimovData_tmp", RooArgSet(*x, w), WeightVar(w)));
    copyDataSet(asimovData_tmp.get(), x, asimovData.get(), x, &w);
    releaseTheGhost(asimovData.get(), x, &w, 1e-9);

    if (mergeBins) {
      string expression = string(x->GetName()) + ">" + to_string(edgeLo) + "&&" + string(x->GetName()) + "<" + to_string(edgeHi);
      asimovData.reset((RooDataSet *)asimovData->reduce(expression.c_str()));
    }

    // if(asimovData_tmp->numEntries()!=x->numBins()){
    //   asimovData->Print("v");
    //   asimovData->get()->Print("v");
    //   convertAsimovToCurve(asimovData.get(), x->GetName(), isHist)->Print("v");
    //   getchar();
    // }
    // convertAsimovToCurve(asimovData.get(), x->GetName(), isHist)->Print("v"); getchar();
    cvec.push_back(convertAsimovToCurve(asimovData.get(), x->GetName(), isHist, mergeBins));
  }
  // RooCurve *band = makeErrorBand(cenCurve,cvec,Z);
  RooCurve *band = cenCurve->makeErrorBand(cvec, Z);
  return band;
}

// _____________________________________________________________________________
RooAbsPdf *createHessePdf(RooFitResult *res, const RooArgSet &params) {
  const TMatrixDSym &V = res->covarianceMatrix();

  EigenMatrixR V_eigen(V.GetNrows(), V.GetNrows());

  for (int i = 0; i < V.GetNrows(); ++i) {
    for (int j = 0; j < V.GetNrows(); ++j) {
      V_eigen(i, j) = V(i, j);
    }
  }

  cpp_dec_float_50 det = V_eigen.determinant();

  if (det <= 0) {
    LOG(logERROR) << "Covariance matrix is not positive definite (|V|=" << det << ") cannot construct p.d.f.";
    exit(-1);
    return 0;
  }

  // Make sure that all given params were floating parameters in the represented
  // fit
  RooArgList params2;
  TIterator *iter = params.createIterator();
  RooAbsArg *arg;
  while ((arg = (RooAbsArg *)iter->Next())) {
    if (res->floatParsFinal().find(arg->GetName())) {
      params2.add(*arg);
    } else {
      LOG(logWARNING) << "Input variable " << arg->GetName() << " was not a floating parameters in fit result and is ignored";
    }
  }
  delete iter;

  // Need to order params in vector in same order as in covariance matrix
  RooArgList params3;
  iter = res->floatParsFinal().createIterator();
  while ((arg = (RooAbsArg *)iter->Next())) {
    if (params2.find(arg->GetName())) {
      params3.add(*arg);
    }
  }
  delete iter;

  // Handle special case of representing full covariance matrix here
  if (params3.getSize() == res->floatParsFinal().getSize()) {

    RooArgList mu;
    for (Int_t i = 0; i < res->floatParsFinal().getSize(); i++) {
      RooRealVar *parclone = (RooRealVar *)res->floatParsFinal().at(i)->Clone(Form("%s_centralvalue", res->floatParsFinal().at(i)->GetName()));
      parclone->setConstant(kTRUE);
      mu.add(*parclone);
    }

    string name = Form("pdf_%s", res->GetName());
    string title = Form("P.d.f of %s", res->GetTitle());

    // Create p.d.f.
    RooAbsPdf *mvg = new RooMultiVarGaussianHighPrecision(name.c_str(), title.c_str(), params3, mu, V);
    mvg->addOwnedComponents(mu);
    return mvg;
  }






  cout << "================================" << endl;
  cout << "getting conditiona pdf????" << endl;
  // exit(-1);



    RooArgList mymu;
    RooArgList params4;
    for (Int_t i = 0; i < res->floatParsFinal().getSize(); i++) {
      RooRealVar *parclone = (RooRealVar *)res->floatParsFinal().at(i)->Clone(Form("%s_centralvalue", res->floatParsFinal().at(i)->GetName()));
      parclone->setConstant(kTRUE);

      if (params3.find(res->floatParsFinal().at(i)->GetName())) {
        mymu.add(*parclone);
        params4.add(*params3.find(res->floatParsFinal().at(i)->GetName()));
      }
    }

    string myname = Form("pdf_%s", res->GetName());
    string mytitle = Form("P.d.f of %s", res->GetTitle());


    const TMatrixDSym &myVred = res->reducedCovarianceMatrix(params4);



    // Create p.d.f.
    RooAbsPdf *mymvg = new RooMultiVarGaussianHighPrecision(myname.c_str(), mytitle.c_str(), params4, mymu, myVred);
    mymvg->addOwnedComponents(mymu);
    return mymvg;











  //                                       -> ->
  // Handle case of conditional p.d.f. MVG(p1|p2) here

  // Find (subset) of parameters that are stored in the covariance matrix
  vector<int> map1, map2;
  for (int i = 0; i < res->floatParsFinal().getSize(); i++) {
    if (params3.find(res->floatParsFinal().at(i)->GetName())) {
      map1.push_back(i);
    } else {
      map2.push_back(i);
    }
  }

  // Rearrange matrix in block form with 'params' first and 'others' last
  // (preserving relative order)
  TMatrixDSym S11, S22;
  TMatrixD S12, S21;
  RooMultiVarGaussianHighPrecision::blockDecompose(V, map1, map2, S11, S12, S21, S22);

  // Calculate offset vectors mu1 and mu2
  RooArgList mu1;
  for (UInt_t i = 0; i < map1.size(); i++) {
    RooRealVar *parclone = (RooRealVar *)res->floatParsFinal().at(map1[i])->Clone(Form("%s_centralvalue", res->floatParsFinal().at(map1[i])->GetName()));
    parclone->setConstant(kTRUE);
    mu1.add(*parclone);
  }

  // Constructed conditional matrix form         -1
  // F(X1|X2) --> CovI --> S22bar = S11 - S12 S22  S21

  // Do eigenvalue decomposition
  TMatrixD S22Inv(TMatrixD::kInverted, S22);

  EigenMatrixR S22_eigen(S22.GetNrows(), S22.GetNrows());

  for (int i = 0; i < S22.GetNrows(); ++i) {
    for (int j = 0; j < S22.GetNrows(); ++j) {
      S22_eigen(i, j) = S22(i, j);
    }
  }

  EigenMatrixR S22Inv_eigen = S22_eigen.inverse();

  for (int i = 0; i < S22.GetNrows(); ++i) {
    for (int j = 0; j < S22.GetNrows(); ++j) {
      S22Inv(i, j) = S22Inv_eigen(i, j).convert_to<double>();
    }
  }

  TMatrixD S22bar = S11 - S12 * (S22Inv * S21);

  // Convert explicitly to symmetric form
  TMatrixDSym Vred(S22bar.GetNcols());
  for (int i = 0; i < Vred.GetNcols(); i++) {
    for (int j = i; j < Vred.GetNcols(); j++) {
      Vred(i, j) = (S22bar(i, j) + S22bar(j, i)) / 2;
      Vred(j, i) = Vred(i, j);
    }
  }
  string name = Form("pdf_%s", res->GetName());
  string title = Form("P.d.f of %s", res->GetTitle());

  // Create p.d.f.
  RooAbsPdf *ret = new RooMultiVarGaussianHighPrecision(name.c_str(), title.c_str(), params3, mu1, Vred);
  ret->addOwnedComponents(mu1);
  return ret;
}

// _____________________________________________________________________________
vector<double> makeHist(TString name, RooCurve *curve) {
  vector<double> values;

  int nbin = curve->GetN() / 2;

  cout << "Create histogram of " << nbin << " bins from " << curve->GetX()[0] << " to " << curve->GetX()[curve->GetN() - 1] << endl;

  for (int i = 0; i < curve->GetN(); i += 2) {
    double central = curve->GetY()[i];
    values.push_back(central);
  }

  return values;
}

// _____________________________________________________________________________
vector<vector<double>> makeHistBand(TString name, RooCurve *curve, RooCurve *band) {
  vector<vector<double>> values;

  int nbin = curve->GetN() / 2;

  cout << "Create histogram of " << nbin << " bins from " << curve->GetX()[0] << " to " << curve->GetX()[curve->GetN() - 1] << endl;

  for (int i = 0; i < curve->GetN(); i += 2) {
    double central = curve->GetY()[i];
    double down = band->GetY()[i];
    double up = band->GetY()[band->GetN() - i - 1];

    vector<double> this_bin = {central, up, down};
    values.push_back(this_bin);
  }

  return values;
}

// _____________________________________________________________________________
// Calculate lower confidence limit
double calcPoissonCLLower(double q, double obs) {
  double LL = 0.;
  if (obs >= 0.) {
    double a = (1. - q) / 2.; // = 0.025 for 95% confidence interval
    LL = TMath::ChisquareQuantile(a, 2. * obs) / 2.;
  }
  return LL;
}

// _____________________________________________________________________________
// Calculate upper confidence limit
double calcPoissonCLUpper(double q, double obs) {
  double UL = 0.;
  if (obs >= 0.) {
    double a = 1. - (1. - q) / 2.; // = 0.025 for 95% confidence interval
    UL = TMath::ChisquareQuantile(a, 2. * (obs + 1.)) / 2.;
  }
  return UL;
}
