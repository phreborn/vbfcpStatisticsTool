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

  char *pngName = Form("asimovFit/asimovFit_%s_%s.png", dScan, cat);

  //cat = "OO_TT_b1";
  cat = Form("OO_%s", cat);

  TCanvas *c = new TCanvas("c", "canvas", 800, 800);

  TFile *f_ws = new TFile(fpath, "read");
  //TFile *f_ws = new TFile("../xmlAnaWSBuilder/run/workspace/vbf_cp_m08/vbf_cp_m08.root", "read");
  RooWorkspace *w = (RooWorkspace*) f_ws->Get("combWS");
  //w->Print();
  RooRealVar *myy = w->var(Form("atlas_invMass_%s", cat));
  //RooRealVar *myy = w->var("atlas_invMass_OO_TT_b1");
  RooSimultaneous *m_pdf = dynamic_cast<RooSimultaneous*>(w->pdf("CombinedPdf")); assert (m_pdf);
  RooCategory *channelCat = (RooCategory*) (&m_pdf->indexCat());
  RooAbsCategoryLValue* m_cat = const_cast<RooAbsCategoryLValue*>(&m_pdf->indexCat());
  RooDataSet *m_data=dynamic_cast<RooDataSet*>(w->data("combDatabinned"));
  TList *m_dataList = m_data->split( *m_cat, true );

  //RooAbsPdf *pdf = m_pdf->getPdf("OO_TT_b1");
  RooAbsPdf *pdf = m_pdf->getPdf(cat);
  RooRealVar *mu = w->var("mu");
  RooRealVar *mu_vbf_sm = w->var("mu_VBF_SM");
  RooRealVar *mu_vbf_rw = w->var("mu_VBF_RW");

  mu->setVal(1.);
  mu_vbf_sm->setVal(1.);
  mu_vbf_rw->setVal(0.);


  RooAbsData *ds_asi = w->data("asimovData_SB_SM");

  RooPlot *myyfr = myy->frame();

//  RooDataSet* datai = ( RooDataSet* )( m_dataList->At( 1 ) ); // not work
//  //datai->plotOn(myyfr);

  //const char *channelCut = Form("channellist==channellist::%s", cat);

		 //ds_asi->plotOn(myyfr, Cut("channellist==channellist::OO_TT_b1"));
		 ds_asi->plotOn(myyfr, Cut(Form("channellist==channellist::%s", cat)));
		
		 //RooCategory channellist("channellist", "channellist"); // a)
		 //channellist.defineType("OO_TT_b1");
		 //channellist.defineType("OO_TT_b2");
		 //channellist.defineType("OO_TT_b3");
		 //channellist.defineType("OO_TT_b4");
		 //channellist.defineType("OO_TT_b5");
		 //channellist.defineType("OO_TT_b6");
		 //pdf->plotOn(myyfr, Slice(channellist, "OO_TT_b1"), ProjWData(channellist, *ds_asi));
		
		 //pdf->plotOn(myyfr, Slice(*channelCat, "OO_TT_b1"), ProjWData(*m_cat, *ds_asi)); // b)
		
		 //double count_asi = ds_asi->sumEntries("channellist==channellist::OO_TT_b1"); cout<<"yield of Asimov data:"<<count_asi<<endl; // c)
		 double count_asi = ds_asi->sumEntries(Form("channellist==channellist::%s", cat)); cout<<"yield of Asimov data:"<<count_asi<<endl; // c)
		 pdf->plotOn(myyfr, Normalization(count_asi, RooAbsReal::NumEvent));


  //m_data->plotOn(myyfr, Cut("channellist==channellist::OO_TT_b1"), DataError(RooAbsData::Poisson));
  //mu->setVal(0.);
  //myy->setRange("SB1", 105000, 120000);
  //myy->setRange("SB2", 130000, 160000);
  //double count_data = m_data->sumEntries("channellist==channellist::OO_TT_b1"); cout<<"yield of Asimov data:"<<count_data<<endl;
  double count_data = m_data->sumEntries(Form("channellist==channellist::%s", cat)); cout<<"yield of Asimov data:"<<count_data<<endl;
  //pdf->plotOn(myyfr);
  //pdf->plotOn(myyfr, Normalization(count_data, RooAbsReal::NumEvent), ProjectionRange("SB1,SB2"));
  //pdf->plotOn(myyfr, ProjectionRange("SB1,SB2"));
  //pdf->plotOn(myyfr, NormRange("SB1,SB2"));

  myyfr->Draw();

  c->SaveAs(pngName);
}
