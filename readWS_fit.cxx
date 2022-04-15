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
  char *cat = argv[2];
  char *dScan = argv[3];
  char *dataset = argv[4];

  char *pngName = Form("asimovFit/asimovFit_%s_%s.png", dScan, cat);
  TString tspngName = TString(pngName);

  //cat = "OO_TT_b1";
  cat = Form("OO_%s", cat);
  TString tscat = TString(cat);

  TCanvas *c = new TCanvas("c", "canvas", 800, 600);

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
  RooDataSet *m_data=dynamic_cast<RooDataSet*>(w->data(dataset));
  TList *m_dataList = m_data->split( *m_cat, true );

  RooDataSet *dataAsi = (RooDataSet *)(m_dataList->FindObject(tscat));

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

  pdf->fitTo(*dataAsi);

  RooPlot *myyfr = myy->frame(Bins(55));
  myyfr->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
  myyfr->SetTitle("Asimov fit in "+tscat);
  //myyfr->SetTitle(0);
  dataAsi->plotOn(myyfr, DataError(RooAbsData::Poisson));
  mu->setVal(0.);
  double nbkg_hat = w->var("nbkg_"+tscat)->getVal();
  pdf->plotOn(myyfr, Normalization(nbkg_hat, RooAbsReal::NumEvent), LineColor(kBlue));
  mu->setVal(1.);
  pdf->plotOn(myyfr, LineColor(kRed));

  std::cout<<"ggH yield: "<<w->var("yield_ggH_"+tscat)->getVal()<<std::endl;
  std::cout<<"Reweighted VBF yield: "<<w->var("yield_VBF_RW_"+tscat)->getVal()<<std::endl;

  //m_data->plotOn(myyfr, Cut("channellist==channellist::OO_TT_b1"), DataError(RooAbsData::Poisson));
  //mu->setVal(0.);
  //myy->setRange("SB1", 105000, 120000);
  //myy->setRange("SB2", 130000, 160000);
  //double count_data = m_data->sumEntries("channellist==channellist::OO_TT_b1"); cout<<"yield of Asimov data:"<<count_data<<endl;
  RooAbsData *ds_asi = w->data(dataset);
  double count_asi = ds_asi->sumEntries(Form("channellist==channellist::%s", tscat.Data())); cout<<"yield of Asimov data:"<<count_asi<<endl;
  RooDataSet *m_sbdata=dynamic_cast<RooDataSet*>(w->data("combDatabinned"));
  double count_data = m_sbdata->sumEntries(Form("channellist==channellist::%s", tscat.Data())); cout<<"yield of sideband data:"<<count_data<<endl;

  myyfr->Draw();

  c->SaveAs(tspngName);
  c->SaveAs(tspngName.ReplaceAll("png", "pdf"));
}
