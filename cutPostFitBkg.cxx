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
#include "RooPlot.h"
#include "RooConstVar.h"
#include "RooSimultaneous.h"
#include "RooAbsCategoryLValue.h"

#include "RooMinimizer.h"
#include "Math/MinimizerOptions.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooLognormal.h"
#include "RooGamma.h"
#include "RooPoisson.h"
#include "RooBifurGauss.h"

#include "RooChebychev.h"
#include "RooCategory.h"
#include "TAxis.h"

#include "TH1.h"
#include "TList.h"

#include "ExtendedModel.hxx"
#include "ExtendedMinimizer.hxx"

#include "log.hxx"
#include "utils.hxx"

#include "atlasrootstyle/AtlasUtils.h"
#include "atlasrootstyle/AtlasLabels.h"
#include "atlasrootstyle/AtlasStyle.h"

using namespace RooFit;

int main(int argc, char** argv)
{
  char *fpath = argv[1];
  const char *cat = argv[2];
  char *dScan = argv[3];
  char *dataset = argv[4];

  char *outName = Form("bkgSRFrac/bkgSRFrac_%s.txt", dScan);
  TString tspngName = TString(outName);

  ofstream ofFrac(outName, ios::out);
  if(!ofFrac){
    ofFrac.close();
    std::cout<<"can't open "<<outName<<std::endl;
  }

  map<TString, pair<float, float>> bins;
  bins["b1"] = make_pair(-999999999, -2);
  bins["b2"] = make_pair(-2, -1);
  bins["b3"] = make_pair(-1, 0);
  bins["b4"] = make_pair(0, 1);
  bins["b5"] = make_pair(1, 2);
  bins["b6"] = make_pair(2, 99999999);

  std::map<TString, std::vector<float>> cats;
  cats["TT"] = {0.14, 1., 0.23, 1.};
  cats["TL"] = {0.14, 1., -1., 0.23};
  cats["LT"] = {-1, 0.14, 0.05, 1.};

  for(auto bin : bins){
    TString oobin = bin.first;
    for(auto c : cats){
      TString bdtcat = c.first;
      cat = (bdtcat+"_"+oobin).Data();
      std::cout<<cat<<std::endl;

      //cat = "OO_TT_b1";
      const char* oricat = cat;
      cat = Form("OO_%s", cat);
      TString tscat = TString(cat);

      TFile *f_ws = new TFile(fpath, "read");
      //TFile *f_ws = new TFile("../xmlAnaWSBuilder/run/workspace/vbf_cp_m08/vbf_cp_m08.root", "read");
      RooWorkspace *w = (RooWorkspace*) f_ws->Get("combWS");
      //w->Print();
      RooRealVar *myy = w->var(Form("atlas_invMass_%s", tscat.Data()));
      //RooRealVar *myy = w->var("atlas_invMass_OO_TT_b1");

      //***** refer to xmlAnaWSBuilder/src/xmlAnaWSBuilder.cc: xmlAnaWSBuilder::Summary()
      RooSimultaneous *m_pdf = dynamic_cast<RooSimultaneous*>(w->pdf("CombinedPdf")); assert (m_pdf);
      RooCategory *channelCat = (RooCategory*) (&m_pdf->indexCat());
      RooAbsCategoryLValue* m_cat = const_cast<RooAbsCategoryLValue*>(&m_pdf->indexCat());

      //RooAbsPdf *pdf = m_pdf->getPdf("OO_TT_b1");
      RooAbsPdf *pdf = m_pdf->getPdf(tscat.Data());
      RooRealVar *mu = w->var("mu");
      RooRealVar *mu_vbf_sm = w->var("mu_VBF_SM");
      RooRealVar *mu_vbf_rw = w->var("mu_VBF_RW");
      RooRealVar *mu_ggF_sm = w->var("mu_ggH_SM");
      RooRealVar *mu_ggF = w->var("mu_ggH");

      mu->setVal(1.); mu->setConstant(1);
      mu_ggF->setVal(1.); mu_ggF->setConstant(1);
      mu_vbf_sm->setVal(0.); mu->setConstant(1);
      mu_ggF_sm->setVal(0.); mu->setConstant(1);

      mu_vbf_rw->setConstant(0); //mu_vbf_rw->setVal(1.);

      unique_ptr<TIterator> iter(w->set("nuisanceParameters")->createIterator());
      RooRealVar *parg=NULL;
      while((parg=dynamic_cast<RooRealVar*>(iter->Next()))){
        if(TString(parg->GetName()).Contains("ATLAS_")) parg->setConstant(true);
      }

      //checking fitResult and post-fit WS
      RooRealVar *bpar0 = w->var(Form("p1_%s", oricat));
      std::cout<<"post-fit WS: "<<bpar0->getVal()<<std::endl;

      RooFitResult *frlt = (RooFitResult*)f_ws->Get("fitResult");
      RooArgList parlist = frlt->floatParsFinal();
      RooAbsReal *par0 = (RooAbsReal*)parlist.find(Form("p1_%s", oricat));
      std::cout<<"fitResult: "<<par0->getVal()<<std::endl;

      mu_vbf_rw->setVal(0);

      myy->setRange("full", 105000, 160000);
      myy->setRange("SBup", 105000, 118000);
      myy->setRange("SR", 118000, 132000);
      myy->setRange("SBdn", 132000, 160000);

      RooAbsReal *int_full = pdf->createIntegral(*myy, NormSet(*myy), Range("full"));
      RooAbsReal *int_SR = pdf->createIntegral(*myy, NormSet(*myy), Range("SR"));
      RooAbsReal *int_SBup = pdf->createIntegral(*myy, NormSet(*myy), Range("SBup"));
      RooAbsReal *int_SBdn = pdf->createIntegral(*myy, NormSet(*myy), Range("SBdn"));
      double fracSR = int_SR->getVal()/int_full->getVal();
      double fracSBup = int_SBup->getVal()/int_full->getVal();
      double fracSBdn = int_SBdn->getVal()/int_full->getVal();

      std::cout<<fracSBup<<"+"<<fracSR<<"+"<<fracSBdn<<"="<<fracSBup+fracSR+fracSBdn<<std::endl;

      ofFrac<<oricat<<","<<fracSR<<std::endl;
    }
  }
}
