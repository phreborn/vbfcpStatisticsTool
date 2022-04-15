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
  char *dataset = argv[2];
  char *outsuffix = argv[3];

  TString outDir = "plotsCombMyy";

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

  TFile *f_ws = new TFile(fpath, "read");
  //TFile *f_ws = new TFile("../xmlAnaWSBuilder/run/workspace/vbf_cp_m08/vbf_cp_m08.root", "read");
  RooWorkspace *w = (RooWorkspace*) f_ws->Get("combWS");
  //w->Print();

  //***** refer to xmlAnaWSBuilder/src/xmlAnaWSBuilder.cc: xmlAnaWSBuilder::Summary()
  RooSimultaneous *m_pdf = dynamic_cast<RooSimultaneous*>(w->pdf("CombinedPdf")); assert (m_pdf);
  RooCategory *channelCat = (RooCategory*) (&m_pdf->indexCat());
  RooAbsCategoryLValue* m_cat = const_cast<RooAbsCategoryLValue*>(&m_pdf->indexCat());
  RooDataSet *m_data=dynamic_cast<RooDataSet*>(w->data(dataset));
  TList *m_dataList = m_data->split( *m_cat, true );

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
    //if(TString(parg->GetName()).Contains("ATLAS_")) parg->setConstant(false);
    parg->setConstant(true);
  }

//  m_pdf->fitTo(*m_data);

  mu_vbf_rw->setConstant(1);

//  while((parg=dynamic_cast<RooRealVar*>(iter->Next()))){
//    parg->setConstant(true);
//  }

  for(auto bin : bins){
    TString oobin = bin.first;
    for(auto c : cats){
      TString bdtcat = c.first;
      const char *cat = (bdtcat+"_"+oobin).Data();
      std::cout<<cat<<std::endl;

      TString tspngName = outDir+"/combfit_myy_"+TString(cat)+TString(outsuffix);

      //cat = "OO_TT_b1";
      const char* oricat = cat;
      cat = Form("OO_%s", cat);
      TString tscat = TString(cat);

      RooRealVar *myy = w->var(Form("atlas_invMass_%s", tscat.Data()));
      //RooRealVar *myy = w->var("atlas_invMass_OO_TT_b1");

      RooDataSet *dataAsi = (RooDataSet *)(m_dataList->FindObject(tscat));

      //RooAbsPdf *pdf = m_pdf->getPdf("OO_TT_b1");
      RooAbsPdf *pdf = m_pdf->getPdf(tscat.Data());

      TCanvas *canv = new TCanvas("canv", "canvas", 800, 600);

      // refer to tutorials/roofit/rf501_simultaneouspdf.C
//      RooPlot *myyfr = myy->frame(Bins(55));
//      myyfr->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
//      myyfr->SetTitle("Asimov fit in "+tscat);
//      m_data->plotOn(myyfr, Cut("channellist==channellist::"+tscat), DataError(RooAbsData::Poisson));
//      m_pdf->plotOn(myyfr, Slice(*channelCat, tscat), ProjWData(*channelCat, *m_data), LineColor(kRed), Invisible());
//      m_pdf->plotOn(myyfr, Slice(*channelCat, tscat), Components("pdf__background_"+tscat), ProjWData(*channelCat, *m_data), LineColor(kBlue));
//      m_pdf->plotOn(myyfr, Slice(*channelCat, tscat), ProjWData(*channelCat, *m_data), LineColor(kRed));

      RooPlot *myyfr = myy->frame(Bins(55));
      myyfr->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
      myyfr->SetTitle("Asimov fit in "+tscat);
      //myyfr->SetTitle(0);
      dataAsi->plotOn(myyfr, DataError(RooAbsData::Poisson));
      pdf->plotOn(myyfr, Components("pdf__background_"+tscat), LineColor(kBlue));
      pdf->plotOn(myyfr, LineColor(kRed));

//      mu->setVal(0.);
//      double nbkg_hat = w->function("yield__background_"+tscat)->getVal();
//      pdf->plotOn(myyfr, Normalization(nbkg_hat, RooAbsReal::NumEvent), LineColor(kBlue));
//      double nVBF_hat = w->function("yield__VBF_RW_"+tscat)->getVal();
//      double nggH_hat = w->function("yield__ggH_"+tscat)->getVal();
//      double nSpur_hat = w->function("yield__spurious_"+tscat)->getVal();
//      double ntot = nVBF_hat+nggH_hat+nSpur_hat+nbkg_hat;
//      mu->setVal(1.);
//      pdf->plotOn(myyfr, Normalization(ntot, RooAbsReal::NumEvent), LineColor(kRed));

//    int obsNBins = 220;
//    unique_ptr<TH1D> hpdf;
//    hpdf.reset(new TH1D("hpdf", "hpdf", obsNBins, myy->getMin(), myy->getMax()));
//    hpdf.reset((TH1D*)pdf->createHistogram("hpdf", *myy));
//    hpdf->Rebin(1);
//    hpdf->Scale(pdf->expectedEvents(RooArgSet(*myy))/hpdf->Integral());
//    for( int ibin = 1 ; ibin <= obsNBins; ibin ++ ) hpdf->SetBinError(ibin, 0);
//    hpdf->Draw();

      myyfr->Draw();
      myyfr->SetMinimum(1e-3); // must be after calling Draw()?

      canv->SaveAs(tspngName+".png");
      canv->SaveAs(tspngName+".pdf");

// ---------------------------------

//      //checking fitResult and post-fit WS
//      RooRealVar *bpar0 = w->var(Form("p1_%s", oricat));
//      std::cout<<"post-fit WS: "<<bpar0->getVal()<<std::endl;
//
//      RooFitResult *frlt = (RooFitResult*)f_ws->Get("fitResult");
//      RooArgList parlist = frlt->floatParsFinal();
//      RooAbsReal *par0 = (RooAbsReal*)parlist.find(Form("p1_%s", oricat));
//      std::cout<<"fitResult: "<<par0->getVal()<<std::endl;
//
//      mu_vbf_rw->setVal(0);
//
//      myy->setRange("full", 105000, 160000);
//      myy->setRange("SBup", 105000, 118000);
//      myy->setRange("SR", 118000, 132000);
//      myy->setRange("SBdn", 132000, 160000);
//
//      RooAbsReal *int_full = pdf->createIntegral(*myy, NormSet(*myy), Range("full"));
//      RooAbsReal *int_SR = pdf->createIntegral(*myy, NormSet(*myy), Range("SR"));
//      RooAbsReal *int_SBup = pdf->createIntegral(*myy, NormSet(*myy), Range("SBup"));
//      RooAbsReal *int_SBdn = pdf->createIntegral(*myy, NormSet(*myy), Range("SBdn"));
//      double fracSR = int_SR->getVal()/int_full->getVal();
//      double fracSBup = int_SBup->getVal()/int_full->getVal();
//      double fracSBdn = int_SBdn->getVal()/int_full->getVal();
//
//      std::cout<<fracSBup<<"+"<<fracSR<<"+"<<fracSBdn<<"="<<fracSBup+fracSR+fracSBdn<<std::endl;
//
//      ofFrac<<oricat<<","<<fracSR<<std::endl;
    }
  }
}
