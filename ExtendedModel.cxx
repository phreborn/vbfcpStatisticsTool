// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Load models from ROOT file and prepare them for fits

#include "RooRealSumPdf.h"
#include "RooAddPdf.h"
#include "RooStarMomentMorph.h"
#include "RooMultiPdf.h"

#include "ExtendedModel.hxx"
#include "utils.hxx"

// _____________________________________________________________________________
// Constructor
ExtendedModel::ExtendedModel( string ModelName, string FileName, string WsName, string ModelConfigName, string DataName, string SnapshotName, bool binnedLikelihood, string TagAsMeasurement, bool FixCache, bool FixMulti )
  :
  TNamed( ModelName.c_str(), ModelName.c_str() ),
  fFileName( FileName ),
  fWsName( WsName ),
  fModelConfigName( ModelConfigName ),
  fDataName( DataName ),
  fSnapshotName( SnapshotName ),
  fBinnedLikelihood( binnedLikelihood ),
  fTagAsMeasurement( TagAsMeasurement ),
  fFixCache( FixCache ),
  fFixMulti( FixMulti )
{

  initialise();

  coutP(InputArguments) << "ExtendedModel::ExtendedModel(" << fName <<") created" << endl;
}

// _____________________________________________________________________________
// Destructor
ExtendedModel::~ExtendedModel()
{
  // TODO
}

// _____________________________________________________________________________
// Load all model information from specified file
void ExtendedModel::initialise()
{
  coutP(InputArguments) << "Opening file " << fFileName << endl;
  fFile = new TFile(fFileName.c_str());

  coutP(InputArguments) << "Loading workspace " << fWsName << endl;
  fWorkspace = (RooWorkspace*)(fFile->Get(fWsName.c_str()));
  if (!fWorkspace) {
    coutE(InputArguments) << "Something went wrong when loading the workspace " << fWsName << endl;
    exit(-1);
  }

  // Fixes for known features
  if (fBinnedLikelihood) {
    coutP(InputArguments) << "Activating binned likelihood evaluation" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
        arg->setAttribute("BinnedLikelihood");
        coutI(InputArguments) << "Activating binned likelihood attribute for " << arg->GetName() << endl;
      }
    }
  }

  if (fTagAsMeasurement != "") {
    coutP(InputArguments) << "Tagging CMS main measurements to reduce memory consumption" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator() ;
    RooAbsArg* arg ;
    while ((arg = iter.next())) {
      if (arg->IsA()==RooAddPdf::Class() && TString(arg->GetName()).BeginsWith(fTagAsMeasurement.c_str())) {
      arg->setAttribute("MAIN_MEASUREMENT") ;
      coutI(InputArguments) << "Component " << arg->GetName() << " is a CMS main measurement";
      }
    }
  }

  if (fFixCache) {
    coutP(InputArguments) << "Fixing cache of early RooStarMomentMorph pdf" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooStarMomentMorph::Class()) {
        ((RooStarMomentMorph*)arg)->fixCache();
        coutI(InputArguments) << "Fixing cache of " << arg->GetName() << endl;
      }
    }
  }

  if (fFixMulti) {
    coutP(InputArguments) << "De-activating level 2 constant term optimization for RooMultiPdf" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooMultiPdf::Class()) {
        arg->setAttribute("NOCacheAndTrack");
        coutI(InputArguments) << "De-activation of level 2 constant term optimization for " << arg->GetName() << endl;
      }
    }
  }

  if (kTRUE) {
    coutP(InputArguments) << "De-activating level 2 constant term optimization for specified pdfs" << endl;
    RooFIter iter = fWorkspace->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      TString aname(arg->GetName());
      if (arg->InheritsFrom(RooAbsPdf::Class()) && (aname.EndsWith("_mm") || aname.Contains("mumu_atlas"))) {
        arg->setAttribute("NOCacheAndTrack");
        coutI(InputArguments) << "De-activation of level 2 constant term optimization for " << arg->GetName();
      }
    }
  }

  // Continue loading the model
  coutP(InputArguments) << "Loading ModelConfig " << fModelConfigName << endl;
  fModelConfig = (ModelConfig*)(fWorkspace->obj(fModelConfigName.c_str()));
  if (!fModelConfig) {
    coutE(InputArguments) << "Something went wrong when loading the ModelConfig " << fModelConfigName << endl;
    exit(-1);
  }

  coutP(InputArguments) << "Grabbing the pdf from the ModelConfig" << endl;
  fPdf = (RooAbsPdf*)fModelConfig->GetPdf();
  if (!fPdf) {
    coutE(InputArguments) << "Something went wrong when loading the pdf" << endl;
    exit(-1);
  }

  coutP(InputArguments) << "Loading ModelConfig " << fModelConfigName << endl;
  fData = (RooAbsData*)(fWorkspace->data(fDataName.c_str()));
  if (!fData) {
    coutE(InputArguments) << "Something went wrong when loading the data set " << fDataName << endl;
    exit(-1);
  }

  coutP(InputArguments) << "Loading the nuisance parameters" << endl;
  fNuis = (RooArgSet*)fModelConfig->GetNuisanceParameters();
  if (!fNuis) {
    coutE(InputArguments) << "Something went wrong when loading the nuisance parameters" << endl;
    exit(-1);
  }

  coutP(InputArguments) << "Loading the global observables" << endl;
  fGlobs = (RooArgSet*)fModelConfig->GetGlobalObservables();
  if (!fGlobs) {
    coutE(InputArguments) << "Something went wrong when loading the global observables" << endl;
    exit(-1);
  }

  coutP(InputArguments) << "Loading the parameters of interest" << endl;
  fPOIs = (RooArgSet*)fModelConfig->GetParametersOfInterest();
  if (!fPOIs) {
    coutE(InputArguments) << "Something went wrong when loading the parameters of interest" << endl;
    exit(-1);
  }

  coutP(InputArguments) << "Loading the observables" << endl;
  fObs = (RooArgSet*)fModelConfig->GetObservables();
  if (!fObs) {
    coutE(InputArguments) << "Something went wrong when loading the observables" << endl;
    exit(-1);
  }

  if (fSnapshotName != "") {
    coutP(InputArguments) << "Loading snapshots" << endl;
    vector<string> parsedSnapshots = parseString(fSnapshotName, ",");
    for (size_t i_snapshot = 0; i_snapshot < parsedSnapshots.size(); ++i_snapshot) {
      string thisSnapshot = parsedSnapshots[i_snapshot];
      coutI(InputArguments) << "Loading snapshot " << thisSnapshot << endl;
      fWorkspace->loadSnapshot(thisSnapshot.c_str());
    }
  }
}

// _____________________________________________________________________________
// Fix all nuisance parameters
void ExtendedModel::fixNuisanceParameters()
{
  for (RooLinkedListIter it = fNuis->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    Double_t value = v->getVal();
    string name = v->GetName();
    coutI(ObjectHandling) << "Fixing nuisance parameter " << name << " at value " << value << endl;
    v->setConstant(1);
  }
}
// _____________________________________________________________________________
// Fix a subset of the nuisance parameters at the specified values
void ExtendedModel::fixNuisanceParameters( string fixName )
{
  vector<string> parsed = parseString(fixName, ",");

  for (size_t i = 0; i < parsed.size(); i++) {
     TString thisName = parsed[i].c_str();
     TString thisVal;
     if (thisName.Contains("[")) {
       assert(thisName.Contains("]"));
       TObjArray* thisNameArray = thisName.Tokenize("[");
       thisName = ((TObjString*)thisNameArray->At(0))->GetString();
       thisVal = ((TObjString*)thisNameArray->At(1))->GetString();
       thisVal.ReplaceAll("]","");
     }

     RooRealVar* par = (RooRealVar*)fWorkspace->var(thisName.Data());
     if (!par) {
       coutE(ObjectHandling) << "Nuisance parameter " << thisName.Data() << " does not exist." << endl;
       exit(-1);
     }

     double value = par->getVal();
     if (thisVal.IsFloat()) {
       value = thisVal.Atof();
       par->setVal(value);
     }

     coutI(ObjectHandling) << "Fixing nuisance parameter " << thisName.Data() << " at value " << value << endl;
     par->setConstant(1);
   }

}

