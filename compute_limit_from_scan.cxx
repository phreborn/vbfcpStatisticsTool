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
    ( "style"    , po::value<char>( &style )->default_value( style )                         , "Test statistics and aymptotics. Choices: t: t_mu; T: t~_mu; y: t_mu(t~_mu asymptotics); Y: t~_mu(t_mu asymptotics); u: q_mu; U: q~_mu." )
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

// _____________________________________________________________________________
void driver(TH1D* h_ll_obs, TH1D* h_ll_asi, char type) {
  Int_t nbin = h_ll_obs->GetNbinsX();

  Bool_t haveTilde(kFALSE);
  if (type=='U' || type=='T' || type=='Y') haveTilde = kTRUE;

  // Harvest for OBS and ASI: mu_hat, L(mu_hat) and L(mu=0)
  Double_t ll_obs_unc(1e30), ll_asi_unc(1e30);
  Double_t ll_obs_mu0, ll_asi_mu0;
  Double_t muhat_obs, muhat_asi;

  for (int i=1; i<=nbin; i++) {
    Double_t mu = h_ll_obs->GetXaxis()->GetBinCenter(i);
    Double_t ll_obs = h_ll_obs->GetBinContent(i);
    Double_t ll_asi = h_ll_asi->GetBinContent(i);

    if (ll_obs<ll_obs_unc) {
      ll_obs_unc = ll_obs;
      muhat_obs = mu;
    }

    if (ll_asi<ll_asi_unc) {
      ll_asi_unc = ll_asi;
      muhat_asi = mu;
    }

    if (mu>0 && h_ll_obs->GetXaxis()->GetBinCenter(i-1)<0) {
      ll_obs_mu0 = ll_obs;
      ll_asi_mu0 = ll_asi;
    }
  }

  LOG(logDEBUG) << "poi_hat observed: " << muhat_obs;
  LOG(logDEBUG) << "L(poi=poi_hat) observed: " << ll_obs_unc;
  LOG(logDEBUG) << "L(poi=0) observed: " << ll_obs_mu0;

  LOG(logDEBUG) << "poi_hat expected: " << muhat_asi;
  LOG(logDEBUG) << "L(poi=poi_hat) expected: " << ll_asi_unc;
  LOG(logDEBUG) << "L(poi=0) expected: " << ll_asi_mu0;

  Double_t z_obs = sqrt((ll_obs_mu0 - ll_obs_unc));
  Double_t p_obs = 1 - ROOT::Math::gaussian_cdf(z_obs);

  Double_t z_asi = sqrt((ll_asi_mu0 - ll_asi_unc));
  Double_t p_asi = 1 - ROOT::Math::gaussian_cdf(z_asi);

  LOG(logDEBUG) << "z observed: " << z_obs;
  LOG(logDEBUG) << "p observed: " << p_obs;

  LOG(logDEBUG) << "z expected " << z_asi;
  LOG(logDEBUG) << "p expected " << p_asi;

  // Axis edges
  Double_t axlo = h_ll_obs->GetXaxis()->GetXmin();
  Double_t axhi = h_ll_obs->GetXaxis()->GetXmax();

  // Setup histograms for plotting
  TH1D* h_qmu   = new TH1D("h_qmu","h_qmu",nbin,axlo,axhi);
  TH1D* h_qmu_A = new TH1D("h_qmu_A","h_qmu_A",nbin,axlo,axhi);
  h_qmu_A->SetLineColor(kRed);

  TH1D* h_sqrtqmu   = new TH1D("h_sqrtqmu","h_sqrtqmu",nbin,axlo,axhi);
  TH1D* h_sqrtqmu_A = new TH1D("h_sqrtqmu_A","h_sqrtqmu_A",nbin,axlo,axhi);
  h_sqrtqmu_A->SetLineColor(kRed);

  // Test statistic
  const char* tsname="unknown";

  switch(type) {
    case 't': tsname = "t_mu"; break;
    case 'T': tsname = "t~_mu"; break;
    case 'y': tsname = "t_mu(t~_mu asymptotics)"; break;
    case 'Y': tsname = "t~_mu(t_mu asymptotics)"; break;
    case 'u': tsname = "q_mu"; break;
    case 'U': tsname = "q~_mu"; break;
  }

  // Book histograms for limits
  LOG(logINFO) << "Booking histograms";

  TH1D* h_pnull = new TH1D("h_pnull", Form("h_pmu (CLs+b) TS=%s", tsname), nbin, axlo, axhi);
  TH1D* h_palt = new TH1D("h_palt", "h_pdom (1-CLb)", nbin, axlo, axhi);
  TH1D* h_pcls = new TH1D("h_pcls", "h_pcls (CLs)", nbin, axlo, axhi);

  h_palt->SetLineColor(kRed);
  h_pcls->SetLineColor(kCyan);
  h_pcls->SetLineWidth(3);

  // Calculation loop
  LOG(logINFO) << "Calculate p-values";

  for (int i=1; i<=nbin; i++) {
    Double_t mu = h_ll_obs->GetXaxis()->GetBinCenter(i);
    Double_t ll_obs = h_ll_obs->GetBinContent(i);
    Double_t ll_asi = h_ll_asi->GetBinContent(i);

    // Starting point q(~)mu = t(~)_mu (eq 8)

    // Denominator is L(mu_hat) unless have tilde and mu_hat<0, in which case it is L(mu=0)
    Double_t ll_denom = ll_obs_unc;

    if ((type=='T' || type=='Y') && muhat_obs<0) ll_denom = ll_obs_mu0;
    if (type=='U' && muhat_obs<0) ll_denom = ll_obs_mu0;

    Double_t qmu = (ll_obs - ll_denom); // SG removing factor 2

    Double_t ll_denom_A = ll_asi_unc;

    if ((type=='T' || type=='Y') && muhat_asi<0) ll_denom_A = ll_asi_mu0;
    if (type=='U' && muhat_asi<0) ll_denom_A = ll_asi_mu0;

    Double_t qmu_A = (ll_asi - ll_denom_A); // SG removing factor 2

    // Now modify qmu for 2-sided (i.e. t(~)), or 1-sided upper/lower type (i.e. q(~))
    switch (type) {
      case 'T':
      case 't':
      case 'Y':
      case 'y':
        // Starting point qmu = t_mu (eq 8)
      break;
      case 'u':
      case 'U':
        // Upper limit qmu = q_mu (eq 14)
        if (muhat_obs > mu) qmu = 0;
        if (muhat_asi > mu) qmu_A = 0;
      break;
    }

    LOG(logDEBUG) << "Bin " << i << ": qmu (observed) = " << qmu;
    LOG(logDEBUG) << "Bin " << i << ": qmu (expected) = " << qmu_A;

    // Histogramming
    h_qmu->SetBinContent(i, qmu);
    h_qmu_A->SetBinContent(i, qmu_A);

    // Calculate sqrt(qmu)
    if (qmu < 0) qmu = 0;
    if (qmu_A < 0) qmu_A = 0;

    double sqrtqmu = sqrt(qmu);
    double sqrtqmu_A = sqrt(qmu_A); // is 'mu/sigma'

    if (!haveTilde || mu>0) {
      h_sqrtqmu->SetBinContent(i, sqrtqmu);
      h_sqrtqmu_A->SetBinContent(i, sqrtqmu_A);
    }

    LOG(logDEBUG) << "Bin " << i << ": sqrt(qmu) (observed) = " << sqrtqmu;
    LOG(logDEBUG) << "Bin " << i << ": sqrt(qmu) (expected) = " << sqrtqmu_A;

    // Asymptotic calculations of p-values
    Double_t pnull, palt;

    switch(type) {
      case 't':
        // tmu: 2-sided limits, eq 35,36 in asymptotics paper ---
        pnull = 2.*ROOT::Math::normal_cdf_c( sqrtqmu, 1.); // WVE CHECKED
        palt = ROOT::Math::normal_cdf_c( sqrtqmu + sqrtqmu_A, 1.) + ROOT::Math::normal_cdf_c( sqrtqmu - sqrtqmu_A, 1.); // WVE CHECKED
      break;

      case 'T':
        // t~mu: F-C limits, eq 43,44 in asymptotics paper
        if (qmu<qmu_A) {
          pnull = 2 * ROOT::Math::normal_cdf_c(sqrtqmu,1); // WVE CHECKED
          palt = ROOT::Math::normal_cdf_c( sqrtqmu + sqrtqmu_A, 1.) + ROOT::Math::normal_cdf_c( sqrtqmu - sqrtqmu_A, 1.); // WVE CHECKED
        } else {
          if (sqrtqmu_A > 1e-30) {
            pnull = ROOT::Math::normal_cdf_c(sqrtqmu,1.) + ROOT::Math::normal_cdf_c( (qmu + qmu_A) / (2 * sqrtqmu_A), 1.); // WVE CHECKED
            palt = ROOT::Math::normal_cdf_c( sqrtqmu_A + sqrtqmu, 1.) + ROOT::Math::normal_cdf_c((qmu - qmu_A) / (2 * sqrtqmu_A), 1.); // WVE CHECKED
          } else {
            pnull = 0;
            palt = 1e-30;
          }
        }
      break;

      case 'Y':
        // tmu: 2-sided limits, eq 35,36 in asymptotics paper ---
        pnull = 2.*ROOT::Math::normal_cdf_c( sqrtqmu, 1.); // WVE CHECKED
        palt = ROOT::Math::normal_cdf_c( sqrtqmu + sqrtqmu_A, 1.) + ROOT::Math::normal_cdf_c( sqrtqmu - sqrtqmu_A, 1.); // WVE CHECKED
      break;

      case 'y':
        // t~mu: F-C limits, eq 43,44 in asymptotics paper
        if (qmu<qmu_A) {
          pnull = 2*ROOT::Math::normal_cdf_c(sqrtqmu,1); // WVE CHECKED
          palt = ROOT::Math::normal_cdf_c( sqrtqmu + sqrtqmu_A, 1.) + ROOT::Math::normal_cdf_c( sqrtqmu - sqrtqmu_A, 1.);  // WVE CHECKED
        } else {
          if (sqrtqmu_A>1e-30) {
            pnull = ROOT::Math::normal_cdf_c(sqrtqmu,1.) + ROOT::Math::normal_cdf_c( (qmu + qmu_A)/(2 * sqrtqmu_A), 1.); // WVE CHECKED
            palt = ROOT::Math::normal_cdf_c( sqrtqmu_A + sqrtqmu, 1.) + ROOT::Math::normal_cdf_c( (qmu - qmu_A)/(2 * sqrtqmu_A), 1.); // WVE CHECKED
          } else {
            pnull = 0;
            palt = 1e-30;
          }
        }
      break;

      case 'u':
        // qmu: 1-sided upper limit, eq 56,57 in asymptotics paper
        pnull = ROOT::Math::normal_cdf_c( sqrtqmu, 1.); // WVE CHECKED
        palt = ROOT::Math::normal_cdf( sqrtqmu_A - sqrtqmu, 1.); // WVE CHECKED
      break;

      case 'U':
        // qmuq~: 1-sided upper limit with bound, eq
        if (qmu<qmu_A) {
          pnull = ROOT::Math::normal_cdf_c( sqrtqmu, 1.); // WVE CHECKED
          palt = ROOT::Math::normal_cdf(sqrtqmu_A - sqrtqmu, 1.); // WVE CHECKED
        } else {
          if (sqrtqmu_A>1e-30) {
            pnull = ROOT::Math::normal_cdf_c( (qmu+qmu_A)/(2*sqrtqmu_A), 1.); // WVE CHECKED
            palt = ROOT::Math::normal_cdf_c( (qmu - qmu_A)/(2 * sqrtqmu_A), 1.); // WVE CHECKED
          } else {
            pnull = 0;
            palt = 1e-30;
          }
        }
      break;
    }

    if (palt < 1e-30) palt = 1e-30;

    LOG(logDEBUG) << "Bin " << i << ": p_null = " << pnull;
    LOG(logDEBUG) << "Bin " << i << ": p_alt = " << palt;
    LOG(logDEBUG) << "Bin " << i << ": p_null / p_alt = " << pnull / palt;

    if (!haveTilde || mu>0) {
      h_pnull->SetBinContent(i, pnull);
      h_palt->SetBinContent(i, palt);
      h_pcls->SetBinContent(i, pnull / palt);
    }
  }

  // Scan for CLs+b transitions over alpha=0.05
  LOG(logINFO) << "Scan for CLs+b transitions over alpha=0.05";

  bool have_up(false);
  bool have_dn(false);
  double mu_05 = 0.0;

  for (int i = 2; i <= nbin; i++) {
    Double_t mu = h_ll_obs->GetXaxis()->GetBinCenter(i);
    Double_t mu_prev = h_ll_obs->GetXaxis()->GetBinCenter(i-1);

    Double_t clsb = h_pnull->GetBinContent(i);
    Double_t clsb_prev = h_pnull->GetBinContent(i-1);

    // Upward transition
    if (clsb>0.05 && clsb_prev<0.05) {
      double grad = (clsb-clsb_prev)/(mu-mu_prev);
      mu_05 = mu_prev + (0.05-clsb_prev)/grad;

      if (!have_up && type != 'u' && type != 'U') {
        // Print first upward transition
        LOG(logWARNING) << "Lower limit on CL(S+B) (95%CL) at mu=" << mu_05;
        have_up = true;
      }
    }

    // Downward transition
    if (clsb_prev>0.05 && clsb<0.05) {
      double grad = (clsb-clsb_prev)/(mu-mu_prev);
      mu_05 = mu_prev + (0.05-clsb_prev)/grad;
      have_dn = true;
    }
  }

  if (have_dn) {
  // Print last downared transition
    LOG(logWARNING) << "Upper limit on CL(S+B) (95%CL) at mu=" << mu_05;
  }

  // Scan for CLs transitions over alpha=0.05
  LOG(logINFO) << "Scan for CLs transitions over alpha=0.05";

  have_up = false;
  have_dn = false;
  mu_05 = 0.0;

  for (int i=2; i<=nbin; i++) {
    Double_t mu = h_ll_obs->GetXaxis()->GetBinCenter(i);
    Double_t mu_prev = h_ll_obs->GetXaxis()->GetBinCenter(i-1);

    Double_t cls = h_pcls->GetBinContent(i);
    Double_t cls_prev = h_pcls->GetBinContent(i-1);

    // Upward transition
    double mu_05;
    if (cls>0.05 && cls_prev<0.05) {
      double grad = (cls-cls_prev)/(mu-mu_prev);
      mu_05 = mu_prev + (0.05-cls_prev)/grad;

      if (!have_up && type != 'u' && type != 'U') {
        // Print first upward transition
        LOG(logWARNING) << "Lower limit on CLs (95%CL) at mu=" << mu_05;
        have_up = true;
      }
    }

    // Downward transition
    if (cls_prev>0.05 && cls<0.05) {
      double grad = (cls-cls_prev)/(mu-mu_prev);
      mu_05 = mu_prev + (0.05-cls_prev)/grad;
      have_dn=true;
    }
  }

  if (have_dn) {
  // Print last downared transition
    LOG(logWARNING) << "Upper limit on CLs (95%CL) at mu=" << mu_05;
  }

  // Debug plot
  TCanvas* c = new TCanvas("c","c",800,800);
  c->Divide(1,2);
  c->cd(1);
  gPad->SetLogy();
  h_pnull->SetMinimum(1e-3);
  h_pnull->SetMaximum(1.1);
  h_pnull->Draw();
  h_pcls->Draw("same");
  h_palt->Draw("same");
  c->cd(2);
  h_sqrtqmu->Draw();
  h_sqrtqmu_A->Draw("same");
}
