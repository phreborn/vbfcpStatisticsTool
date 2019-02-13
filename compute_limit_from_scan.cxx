// Author      : Stefan Gadatsch, Wouter Verkerke
// Email       : stefan.gadatsch@cern.ch, verkerke@nikhef.nl
// Date        : 2019-02-13
// Description : Compute upper limit from profile likelihood scan

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Math/DistFuncMathCore.h"

#include "log.hxx"
#include "utils.hxx"

#include "boost/program_options.hpp"
#include "boost/program_options/cmdline.hpp"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/variables_map.hpp"

using namespace std;

// _____________________________________________________________________________
// Declarations of functions used in this file

void driver(TH1D* h_ll_obs, TH1D* h_ll_asi, char type);

// _____________________________________________________________________________
// Main routine
int main(int argc, char** argv)
{
  TTime thistime = gSystem->Now();

  // Input information
  string inFileName = "path/to/file/with/graphs.root";
  string observedScanName = "name_of_observed_scan";
  string asimovScanName = "name_of_asimov_scan";
  int nbins = 1000;
  double lowEdge = -0.5;
  double highEdge = 0.5;
  double boundary = 0.0;
  string limitType = "upper";
  char style = 't';

  // Misc settings
  string loglevel = "INFO";

  // Bookkeeping

  using namespace boost;
  namespace po = boost::program_options;
  po::options_description desc( "Program options" );
  desc.add_options()
    ( "help      , h"                                                                        , "Print this help message")
    ( "input"    , po::value<string>( &inFileName )->default_value( inFileName )             , "File holding TGraphs with likelihood scans." )
    ( "observed" , po::value<string>( &observedScanName )->default_value( observedScanName ) , "Name of observed likelihood scan." )
    ( "asimov"   , po::value<string>( &asimovScanName )->default_value( asimovScanName )     , "Name of expected likelihood scan." )
    ( "nbins"    , po::value<int>( &nbins )->default_value( nbins )                          , "Number of bins used for sampling histogram from graph." )
    ( "low"      , po::value<double>( &lowEdge )->default_value( lowEdge )                   , "Low edge for booking histogram." )
    ( "high"     , po::value<double>( &highEdge )->default_value( highEdge )                 , "High edge for booking histogram." )
    ( "boundary" , po::value<double>( &boundary )->default_value( boundary )                 , "Value of the physical boundary, if present." )
    ( "type"     , po::value<string>( &limitType )->default_value( limitType )               , "Upper or lower limit." )
    ( "style"    , po::value<char>( &style )->default_value( style )                         , "Test statistics and aymptotics." )
    ( "loglevel" , po::value<string>( &loglevel )->default_value( loglevel )                 , "POIs to use." )
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

  // Create histograms from likelihood scans
  int sign = 1;
  if (limitType == "upper") {
    sign = 1;
  } else if (limitType == "lower") {
    sign = -1;
  } else {
    LOG(logERROR) << "Unknown limit type " << limitType << ". Choose upper or lower.";
  }

  LOG(logINFO) << "Sampling histograms from likelihood scans in the range [" << lowEdge << "," << highEdge << "] with " << nbins << " bins.";
  TH1D* h_ll_obs = new TH1D("h_ll_obs", "h_ll_obs", nbins, lowEdge, highEdge);
  TH1D* h_ll_asi = new TH1D("h_ll_obs", "h_ll_obs", nbins, lowEdge, highEdge);

  // Load graphs from scan and fill histogram
  LOG(logDEBUG) << "Reading observed graph " << observedScanName << " from file " << inFileName << ".";

  TFile* f = TFile::Open(inFileName.c_str());
  TGraph* g = (TGraph*) f->Get(observedScanName.c_str());

  for (int i=1; i<=nbins; i++) {
    Double_t mu = h_ll_obs->GetXaxis()->GetBinCenter(i);
    h_ll_obs->SetBinContent(i, g->Eval(sign * (mu - boundary), 0, "S"));
  }

  delete f;
  LOG(logDEBUG) << "Oberved histograms filled.";

  LOG(logDEBUG) << "Reading expected graph " << asimovScanName << " from file " << inFileName << ".";
  f = TFile::Open(inFileName.c_str());
  g = (TGraph*) f->Get(asimovScanName.c_str());

  for (int i=1; i<=nbins; i++) {
    Double_t mu = h_ll_asi->GetXaxis()->GetBinCenter(i);
    h_ll_asi->SetBinContent(i, g->Eval(sign * (mu - boundary), 0, "S"));
  }

  delete f;
  LOG(logDEBUG) << "Expected histogram filled.";

  // Debug plot for validation
  h_ll_obs->SetLineColor(kRed);
  h_ll_obs->SetLineColor(kAzure);
  h_ll_obs->Draw();
  h_ll_asi->Draw("same");
  h_ll_obs->SetTitle(";POI;-2 ln #Lambda;");

  TLegend* leg = new TLegend(0.1, 0.7, 0.48, 0.9);
  leg->AddEntry(h_ll_obs, "observed", "l");
  leg->AddEntry(h_ll_asi, "expected", "l");
  leg->Draw("same");

  // Run limit computation
  driver(h_ll_obs, h_ll_asi, style);
}


