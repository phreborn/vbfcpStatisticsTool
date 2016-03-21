// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Load models from ROOT file and prepare them for fits

#ifndef EXTENDEDMODEL
#define EXTENDEDMODEL

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TNamed.h"
#include "TFile.h"
#include "Math/MinimizerOptions.h"

#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"

#include "RooStats/ModelConfig.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

#ifdef __MAKECINT__
#pragma link C++ class ExtendedModel;
#endif

class ExtendedModel : public TNamed {

// _____________________________________________________________________________
public:

  // Constructor and destructor
  ExtendedModel( string ModelName, string FileName, string WsName,
                 string ModelConfigName, string DataName, string SnapshotName,
                 bool binnedLikelihood = true, string TagAsMeasurement = "pdf_",
                 bool FixCache = true, bool FixMulti = true );
  virtual ~ExtendedModel();
  RooWorkspace* GetWorkspace() { return fWorkspace; }
  ModelConfig* GetModelConfig() { return fModelConfig; }
  RooAbsPdf* GetPdf() { return fPdf; }
  RooAbsData* GetData() { return fData; }
  RooArgSet* GetNuisanceParameters() { return fNuis; }
  RooArgSet* GetGlobalObservables() { return fGlobs; }
  RooArgSet* GetParametersOfInterest() { return fPOIs; }
  RooArgSet* GetObservables() { return fObs; }

  void fixNuisanceParameters();
  void fixNuisanceParameters( string fixName );

// _____________________________________________________________________________
protected:

  void initialise();

// _____________________________________________________________________________
private:

  string fFileName;
  string fWsName;
  string fModelConfigName;
  string fDataName;
  string fSnapshotName;
  bool fBinnedLikelihood;
  string fTagAsMeasurement;
  bool fFixCache;
  bool fFixMulti;

  TFile* fFile;
  RooWorkspace* fWorkspace;
  ModelConfig* fModelConfig;
  RooAbsPdf* fPdf;
  RooAbsData* fData;
  RooArgSet* fNuis;
  RooArgSet* fGlobs;
  RooArgSet* fPOIs;
  RooArgSet* fObs;

// _____________________________________________________________________________
protected:

  // ClassDef(ExtendedModel, 1)

};

#endif
