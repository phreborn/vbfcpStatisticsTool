// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2014-10-23
// Description : Rank parameters iterative based on a reduced Hesse matrix
//               - Find the parameter which changes the Hesse error on specified
//                 POIs the least and removes it from the covariance matrix
//               - Repeat procedure until only one parameter remains or the
//                 change is above a specified threshold
//               - Global ranking uses error propagation in case multiple POIs are
//                 specified to define order of parameters. The final error is
//                 computed for weighted average of POIs with the passed weights
//                 (default is equal weights, other option could be SM Higgs
//                 cross-section)
//               - Filter on the parameter name in order to run the pruning only
//                 on a subset of all parameters
//               - Possibility to specify multiple filters that define sets of
//                 parameters that can be pruned in given order, for example first
//                 MC statistics, then all others
//               - Only parameters that change individual POIs below a given
//                 threshold are marked for pruning. If any POI passes the threshold,
//                 the parameter is kept
//               - 'auto' triggers the automatic determination of thresholds for
//                 pruning based on th PDG rules for rounding.
//               - Option to keep more than just the significant digits (default is
//                 one more)
//               - If only one weight or threshold is specified, it will be used for
//                 all POIs
//               - Example configuration: 'auto:5' and '2' additional digits define the
//                 threshold for pruning. auto tells the algorithm to base the decision on
//                 the actual precision quoted after rounding following the PDG rules. 5 is
//                 the threshold that the selected digit is allowed to vary and 2 defines the
//                 (non-)significant digit to check. So in this example the 2nd non-significant
//                 digit after rounding is allowed to vary by up to 5 units.
//               - Pass a list of parameters that can be pruned to not start from scratch
//                 in case a tighter threshold should be probed

#include <iomanip>
#include <stdlib.h>
#include <list>
#include<iostream>
#include<sstream>

#include "TFile.h"
#include "TMatrixDSym.h"
#include "TRegexp.h"
#include "TArrayI.h"
#include "Math/MinimizerOptions.h"

#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"

#include "RooStats/ModelConfig.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

// ____________________________________________________________________________|__________
// Global options, touch only if you know what you're doing!
string minimizerType   = "Minuit2";
string minimizerAlgo   = "Migrad";
int defaultStrategy    = 0;
bool binnedLikelihood  = 1;
bool offsetting        = 1;
int constOpt           = 2;
double eps             = 1.0;
int printLevel         = 1;

// ____________________________________________________________________________|__________
void runPruning();
void runPruning( const string& inFileName,
                 const string& wsName = "combined",
                 const string& modelConfigName = "ModelConfig",
                 const string& dataName = "obsData",
                 const string& pruningPoi = "mu_ggF,mu_VBF,mu_WH,mu_ZH,mu_ttH",
                 const string& pruningFilter = ".*",
                 const string& pruningWeight = "1.0",
                 const string& pruningThreshold = "auto:5",
                 int pruningAdditionalDigit = 2,
                 list< string > prePrunedParameters = list< string >() );
void runPruningOnResult( const string& inFileName,
                         const string& resultName,
                         const string& pruningPoi = "mu_ggF,mu_VBF,mu_WH,mu_ZH,mu_ttH",
                         const string& pruningFilter = ".*",
                         const string& pruningWeight = "1.0",
                         const string& pruningThreshold = "auto:5",
                         int pruningAdditionalDigit = 2,
                         list< string > prePrunedParameters = list< string >() );
list< string > PruneNuisanceParameters(const TMatrixDSym chesse,
                                       RooFitResult* fitresult,
                                       const string& poi = "mu_ggF,mu_VBF,mu_WH,mu_ZH,mu_ttH",
                                       const string& filter = ".*",
                                       const string& weight = "1.0",
                                       const string& threshold = "auto:5",
                                       int additionalDigit = 2,
                                       list< string > prePrunedParameters = list< string >());
void RemoveParameter( TMatrixDSym& hes,
                      RooArgList& pars,
                      list< string > names );
void PrintRanking( set< pair< double, string > > uncerts,
                   double initTotalError );
pair< double, double > PDGrounding( double value,
                                    double error,
                                    int additionalDigit = 1 );
int GetThreeDigits( double error );
int GetNSigDigits( int threeDigits );
double frexp10( double x,
                int* exp );
double FormatValue( double value,
                    int exponent,
                    int nDigits,
                    int extraRound = 0 );
bool AlmostEqualUlpsAndAbs( float A,
                            float B,
                            float maxDiff,
                            int maxUlpsDiff );

// ____________________________________________________________________________|__________
// Dummy function for compilation
void runPruning()
{

}

// ____________________________________________________________________________|__________
// Wrapper that loads a model from a workspace, performs a fit to find the
// Hesse matrix and initiates the pruning procedure
void runPruning( const string& inFileName,
                 const string& wsName,
                 const string& modelConfigName,
                 const string& dataName,
                 const string& pruningPoi,
                 const string& pruningFilter,
                 const string& pruningWeight,
                 const string& pruningThreshold,
                 int pruningAdditionalDigit,
                 list< string > prePrunedParameters )
{
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerType.c_str(), minimizerAlgo.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

  cout << "Loading model from file." << endl;

  TFile f(inFileName.c_str());
  RooWorkspace* ws = (RooWorkspace*)(f.Get(wsName.c_str()));
  if (!ws) {
    cout << "ERROR::Workspace " << wsName << " doesn't exist!" << endl;
    exit(1);
  }

  if (binnedLikelihood) {
    RooFIter iter = ws->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
        arg->setAttribute("BinnedLikelihood");
        cout << "Activating binned likelihood for " << arg->GetName() << endl;
      }
    }
  }

  ModelConfig* mc = (ModelConfig*)(ws->obj(modelConfigName.c_str()));
  if (!mc) {
    cout << "ERROR::ModelConfig " << modelConfigName << " doesn't exist!" << endl;
    exit(1);
  }

  RooAbsPdf* pdf = (RooAbsPdf*)mc->GetPdf();
  if (!pdf) {
    cout << "ERROR::PDF not found!" << endl;
    exit(1);
  }

  RooAbsData* data = (RooAbsData*)(ws->data(dataName.c_str()));
  if (!data) {
    cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
    exit(1);
  }

  RooArgSet* nuis = (RooArgSet*)mc->GetNuisanceParameters();
  if (!nuis) {
    cout << "ERROR::Nuisance parameter set doesn't exist!" << endl;
    exit(1);
  }

  RooArgSet* globs = (RooArgSet*)mc->GetGlobalObservables();
  if (!globs) {
    cout << "ERROR::Global observables set doesn't exist!" << endl;
    exit(1);
  }

  cout << "Performing fit to compute Hesse matrix for pruning." << endl;

  RooArgSet* allParameters = pdf->getParameters(*data);
  ws->saveSnapshot("tmpBeforePruning", *allParameters);

  for (RooLinkedListIter it = mc->GetParametersOfInterest()->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    v->setConstant(1);
  }

  TString allPoi = pruningPoi;
  allPoi.ReplaceAll(" ", "");
  TObjArray* allPoiArray = allPoi.Tokenize(",");
  unsigned int numPoi = allPoiArray->GetEntries();
  for (unsigned int itrPoi = 0; itrPoi < numPoi; ++itrPoi) {
    TString thisPoi = ((TObjString*)allPoiArray->At(itrPoi))->GetString();
    ws->var(thisPoi.Data())->setConstant(0);
  }

  // Fit of pdf to data to obtain the Hesse matrix
  RooNLLVar* nll = (RooNLLVar*)pdf->createNLL(*data, Constrain(*nuis), GlobalObservables(*globs), Offset(1), NumCPU(1, RooFit::Hybrid));
  RooMinimizer* minimizer = new RooMinimizer(*nll);
  minimizer->setPrintLevel(printLevel);
  minimizer->optimizeConst(constOpt);
  minimizer->setMinimizerType(minimizerType.c_str());
  minimizer->setProfile(1);
  minimizer->setStrategy(defaultStrategy);
  minimizer->setEps(eps);
  minimizer->minimize(minimizerType.c_str(), minimizerAlgo.c_str());
  minimizer->hesse();

  TMatrixDSym hesse = ((TMatrixDSym)minimizer->lastMinuitFit()->covarianceMatrix()).Invert();

  string name = Form("fitresult_%s_%s", pdf->GetName(), data->GetName());
  string title = Form("Result of fit of p.d.f. %s to dataset %s", pdf->GetName(), data->GetName());
  RooFitResult* fitresult = minimizer->save(name.c_str(), title.c_str());

  // TFile of("result.root", "recreate");
  // fitresult->Write("",TObject::kOverwrite);
  // of.Close();

  list< string > prunedNuisanceParameters = PruneNuisanceParameters(hesse, fitresult, pruningPoi, pruningFilter, pruningWeight, pruningThreshold, pruningAdditionalDigit, prePrunedParameters);

  ws->loadSnapshot("tmpBeforePruning");
}

// ____________________________________________________________________________|__________
// Wrapper function to load a RooFitResult from a file and pass it to the
// ranking function
void runPruningOnResult( const string& inFileName,
                         const string& resultName,
                         const string& pruningPoi,
                         const string& pruningFilter,
                         const string& pruningWeight,
                         const string& pruningThreshold,
                         int pruningAdditionalDigit,
                         list< string > prePrunedParameters )
{
  TFile f(inFileName.c_str());

  RooFitResult* fitresult = NULL;
  f.GetObject(resultName.c_str(), fitresult);
  if (!fitresult) {
    cout << "Error::Could not find fitresult " << resultName << endl;
    exit(-1);
  }

  const TMatrixDSym ccov = fitresult->covarianceMatrix();
  TMatrixDSym cov(ccov);
  TMatrixDSym hesse(cov.Invert());

  list< string > prunedNuisanceParameters = PruneNuisanceParameters(hesse, fitresult, pruningPoi, pruningFilter, pruningWeight, pruningThreshold, pruningAdditionalDigit, prePrunedParameters);
}

// ____________________________________________________________________________|__________
// Ranking/pruning function. Compute symmetric uncertainty on POIs using the
// Hessian matrix. Remove iteratively a single NP and compare uncertainty from
// reduced Hessian matrix to initial one. Multiple POIs are combined by computing
// the weighted average and propagating the partial uncertainties and taking the
// OR of the individual rankings.
list< string > PruneNuisanceParameters( const TMatrixDSym chesse,
                                        RooFitResult* fitresult,
                                        const string& poi,
                                        const string& filter,
                                        const string& weight,
                                        const string& threshold,
                                        int additionalDigit,
                                        list< string > prePrunedParameters )
{
  // Decompose POI names and associated weights and thresholds
  TString allPoi = poi;
  TString allWeight = weight;
  TString allThreshold = threshold;
  TString allFilter = filter;
  allPoi.ReplaceAll(" ", "");
  allWeight.ReplaceAll(" ", "");
  allThreshold.ReplaceAll(" ", "");
  allFilter.ReplaceAll(" ", "");
  TObjArray* allPoiArray = allPoi.Tokenize(",");
  TObjArray* allWeightArray = allWeight.Tokenize(",");
  TObjArray* allThresholdArray = allThreshold.Tokenize(",");
  TObjArray* allFilterArray = allFilter.Tokenize(",");
  unsigned int numPoi = allPoiArray->GetEntries();
  unsigned int numWeight = allWeightArray->GetEntries();
  unsigned int numThreshold = allThresholdArray->GetEntries();
  unsigned int numFilter = allFilterArray->GetEntries();

  // Check if number of POIs and weights match
  if (numPoi != numWeight && numWeight != 1) {
    cout << "PruneNuisanceParameters() Number of POIs does not match number of specified weights. No additional pruning." << endl;
    return prePrunedParameters;
  }

  // Check if number of POIs and thresholds match
  if (numPoi != numThreshold && numThreshold !=1) {
    cout << "PruneNuisanceParameters() Number of POIs does not match number of specified thresholds. No additional pruning." << endl;
    return prePrunedParameters;
  }

  // Store information in easier usable objects
  vector<TString> PoiList;
  map<string, double> WeightList;
  map<string, TString> ThresholdList;
  map<string, double> ThresholdCutList;
  vector<TString> FilterList;
  double totalWeights = 0.0;

  for (unsigned int itrPoi = 0; itrPoi < numPoi; ++itrPoi) {
    TString thisPoi = ((TObjString*)allPoiArray->At(itrPoi))->GetString();
    TString thisWeight = (numWeight == 1) ? ((TObjString*)allWeightArray->At(0))->GetString() : ((TObjString*)allWeightArray->At(itrPoi))->GetString();
    TString thisThreshold = (numThreshold == 1) ? ((TObjString*)allThresholdArray->At(0))->GetString() : ((TObjString*)allThresholdArray->At(itrPoi))->GetString();
    double tmpWeight = atof(thisWeight.Data());

    TObjArray* thisThresholdArray = thisThreshold.Tokenize(":");
    TString thisSaveThreshold = ((TObjString*)thisThresholdArray->At(0))->GetString();
    double thisThresholdCut = (thisThresholdArray->GetEntries() == 1) ? 0.0 : atof(((TObjString*)thisThresholdArray->At(1))->GetString().Data());

    PoiList.push_back(thisPoi);
    WeightList[thisPoi.Data()] = tmpWeight;
    ThresholdList[thisPoi.Data()] = thisSaveThreshold;
    ThresholdCutList[thisPoi.Data()] = thisThresholdCut;

    totalWeights += tmpWeight;

    cout << "PruneNuisanceParameters() Parsed POI " << thisPoi.Data()
         << ", assigning weight " << tmpWeight
         << " and threshold " << thisThreshold.Data() << endl;
  }

  if (numFilter > 0) {
    cout << "PruneNuisanceParameters() Nuisance parameters will be pruned in the following order:";
    for (unsigned int itrSet = 0; itrSet < numFilter; ++itrSet) {
      TString thisSet = ((TObjString*)allFilterArray->At(itrSet))->GetString();
      FilterList.push_back(thisSet);
      cout << " " << thisSet.Data() << ";";
    }
    cout << endl;
  } else {
    cout << "PruneNuisanceParameters() Include all parameters in ranking" << endl;
    FilterList.push_back(".*");
  }

  // Clean up
  delete allPoiArray;
  delete allWeightArray;
  delete allThresholdArray;
  delete allFilterArray;

  // Normalize weights
  for (map<string, double>::iterator itrWeight = WeightList.begin(); itrWeight != WeightList.end(); ++itrWeight) {
    WeightList[itrWeight->first] = itrWeight->second / totalWeights;
  }

  // Print the RooFitResult used for pruning
  fitresult->Print();

  // Get initial information, such as
  //   - covariance matrix
  //   - floating parameters
  //   - position of POIs in the list of floating parameters
  //   - best fit and initial Hesse error
  // TMatrixDSym cov = fitresult->covarianceMatrix();
  TMatrixDSym hes(chesse);
  TMatrixDSym cov(hes.Invert());
  RooArgList pars = fitresult->floatParsFinal();

  map<string, int> index;
  map<string, double> initErr;
  map<string, double> bestFit;
  map<string, double> initErrSig;

  // map<string, double> ThresholdList;
  for (vector<TString>::const_iterator itr = PoiList.begin(), end = PoiList.end(); itr != end; ++itr) {
    int thisIndex = pars.index(itr->Data());
    double thisInitErr = sqrt(cov[thisIndex][thisIndex]);
    double thisBestFit = ((RooRealVar*)pars.at(thisIndex))->getVal();
    pair<double, double> tmpVals = PDGrounding(thisBestFit, thisInitErr, additionalDigit);

    index[itr->Data()] = thisIndex;
    initErr[itr->Data()] = thisInitErr;
    initErrSig[itr->Data()] = tmpVals.second;
    bestFit[itr->Data()] = thisBestFit;

    cout << "PruneNuisanceParameters() Harvest parameter (" << thisIndex << "): " << *itr
         << " with Hesse error " << tmpVals.first << " (" << thisBestFit
         << ") +/- " << tmpVals.second << " (" << thisInitErr << ")" << endl;
  }

  // Compute the weighted average of POIs
  double totalBestFit = 0.0;
  for (vector<TString>::const_iterator itr = PoiList.begin(), end = PoiList.end(); itr != end; ++itr) {
    totalBestFit += WeightList[itr->Data()] * bestFit[itr->Data()];
  }

  // Compute uncertainty on the weighted average of the POIs using error propagation.
  // This takes into account the two-point correlations between the POIs.
  double initTotalErrorSquared = 0.0;

  for (vector<TString>::const_iterator itr = PoiList.begin(), end = PoiList.end(); itr != end; ++itr) {
    for (vector<TString>::const_iterator jtr = PoiList.begin(), end = PoiList.end(); jtr != end; ++jtr) {
      initTotalErrorSquared += WeightList[itr->Data()] * WeightList[jtr->Data()] * cov[index[itr->Data()]][index[jtr->Data()]];
    }
  }
  double initTotalError = sqrt(initTotalErrorSquared);

  cout << "PruneNuisanceParameters() Initial averaged POI with propagated uncertainties: "
       <<  totalBestFit <<  " +/- " << initTotalError << endl;;

  // Ranking starts here.
  // Iteratively remove parameters and recompute Hesse uncertainty on the POI and on
  // the averaged parameter used for the final ranking. Maybe filter list of parameters
  // to prune base on regexp.
  list< string > order = prePrunedParameters;
  unsigned int allParams = pars.getSize();
  for (list<string>::iterator itr = order.begin(); itr != order.end(); ++itr) {
    string name(*itr);
    int par2rem_index = pars.index(name.c_str());
    pars.remove(*pars.at(par2rem_index));
  }
  set< pair< double, string > > uncertsPmO;

  for (vector<TString>::const_iterator itrFilter = FilterList.begin(), end = FilterList.end(); itrFilter != end; ++itrFilter) {
    TString thisFilter = *itrFilter;
    string stringFilter = thisFilter.Data();

    TRegexp reg(thisFilter);
    Ssiz_t dummy(0);

    while ((unsigned int)order.size() < allParams - (numPoi + 1)) {
      set< pair< double, string > > uncerts;
      map< string, set< pair< double, string > > > part_uncerts;

      for (RooLinkedListIter it = pars.iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
        string name = v->GetName();

        // Skip removing POIs
        if (allPoi.Contains(name)) {
          continue;
        }

        // Apply filtering on a subset of parameters
        if (stringFilter != "" && reg.Index(TString(v->GetName()), &dummy, 0) == -1) {
          continue;
        }

        // Attach current parameter to a temporary list of parameters that should be removed from
        // Hesse matrix. This includes the previously ranked, and already removed parameters.
        list< string > names = order;
        names.push_back(name);

        // Grab covariance matrix and floating parameters for doing calculations. These objects will
        // be modified for every parameter, and thus need to be reverted to the original ones in
        // in every iteration.
        // TMatrixDSym tmp_cov = fitresult->covarianceMatrix();
        TMatrixDSym tmp_hesse(chesse);
        RooArgList tmp_pars = fitresult->floatParsFinal();

        // Reduce temporary Hesse matrix by a set of parameters and keep track of the parameters left
        RemoveParameter(tmp_hesse, tmp_pars, names);
        TMatrixDSym tmp_cov = tmp_hesse.Invert();

        // Harvest and store reduced components
        map<string, int > index_red;
        // map<string, double > errred;
        for (vector<TString>::const_iterator itr = PoiList.begin(), end = PoiList.end(); itr != end; ++itr) {
          int thisIndex = tmp_pars.index(itr->Data());
          double thisErr = sqrt(tmp_cov[thisIndex][thisIndex]);

          index_red[itr->Data()] = thisIndex;
          // errred[itr->Data()] = thisErr;

          // Store individual uncertainties for reference
          part_uncerts[itr->Data()].insert(make_pair(thisErr, name));
        }

        // Compute uncertainty on the weighted average of the POIs using error propagation.
        // This takes into account the two-point correlations between the POIs.
        double redTotalErrorSquared = 0.0;
        for (vector<TString>::const_iterator i = PoiList.begin(), end = PoiList.end(); i != end; ++i) {
          for (vector<TString>::const_iterator j = PoiList.begin(), end = PoiList.end(); j != end; ++j) {
            redTotalErrorSquared += WeightList[i->Data()] * WeightList[j->Data()] * tmp_cov[index_red[i->Data()]][index_red[j->Data()]];
          }
        }
        double redTotalError = sqrt(redTotalErrorSquared);

        // Store total uncertainty for the ranking for remaining parameters
        uncerts.insert(make_pair(redTotalError, name));
      }
      if (uncerts.empty()) break;

      // Print ranking of remaining parameters based on uncertainty on combined inclusive POI
      cout << endl;
      cout << "PruneNuisanceParameters() Ranking of remaining parameters based on combined inclusive POI" << endl;
      PrintRanking(uncerts, initTotalError);

      // Find the lowest ranked parameter which changes none of the uncertainties of the individual
      // parameters above the given threshold
      bool foundPar2Rem = false;
      string par2rem = uncerts.rbegin()->second;
      double par2rem_err = uncerts.rbegin()->first;
      unsigned int tmpSum = pars.getSize() + 1 - numPoi;

      // Loop over all remaining parameters in the global ranking
      for (set<pair<double, string> >::reverse_iterator rankitr = uncerts.rbegin(); rankitr != uncerts.rend(); ++rankitr) {
        par2rem = rankitr->second;
        par2rem_err = rankitr->first;
        bool passThreshold = true;
        int tmpIndex = 0;

        cout << "PruneNuisanceParameters() --> Testing " << par2rem << endl;;

        // Loop over the individual POIs to test effect of removing proposed parameter on each of them
        for(map<string, set< pair< double, string > > >::iterator iterator = part_uncerts.begin(); iterator != part_uncerts.end(); ++iterator) {
          string thisPoi = iterator->first;

          cout << "PruneNuisanceParameters()  | " << thisPoi;

          set< pair< double, string > > thisRanking = iterator->second;
          unsigned int thisPosition = 1;
          double thisReduction = 1.0;
          double thisInitError = initErr[thisPoi];
          TString thisThreshold = ThresholdList[thisPoi];

          // Find the proposed parameter, the position in the individual ranking and the estimated uncertainty reduction
          for (set<pair<double, string> >::reverse_iterator subitr = thisRanking.rbegin(); subitr != thisRanking.rend(); ++subitr) {
            string varName(subitr->second);
            double thisUncert = subitr->first;

            // Check the individual uncertainty reduction if the parameter is found and exit search
            if (varName == par2rem) {
              thisReduction = thisUncert/thisInitError - 1;

              // Check if the threshold is passed, either the relative reduction or check
              // whether quoted uncertainty would be changed in case of 'auto'
              if (thisThreshold != "auto") {
                double tmpThreshold = atof(thisThreshold.Data());
                if (fabs(thisReduction) > tmpThreshold) {
                  passThreshold = false;
                }
              } else {
                pair<double, double> tmpVals = PDGrounding(bestFit[thisPoi], thisUncert, additionalDigit);

                ostringstream streamInit, streamRed;
                streamInit << initErrSig[thisPoi];
                streamRed << tmpVals.second;
                TString strInit(streamInit.str());
                TString strRed(streamRed.str());

                strInit.ReplaceAll(".", "");
                strRed.ReplaceAll(".", "");

                while (strInit.Length() < strRed.Length()) {
                  strInit.Append("0");
                }

                while (strRed.Length() < strInit.Length()) {
                  strRed.Append("0");
                }

                double lastInit = atof(strInit.Data());
                double lastRed = atof(strRed.Data());

                cout << ": absolute change " << initErrSig[thisPoi] << " --> " << tmpVals.second;

                // Allow the specified non-significant digit to change by a specified value
                if (tmpVals.second < initErrSig[thisPoi]) {
                  if (additionalDigit > 0) {
                    if (!AlmostEqualUlpsAndAbs(lastInit, lastRed, ThresholdCutList[thisPoi], 4)) {
                      cout << " reduced too much: " << lastInit << " - " << lastRed << " > " << ThresholdCutList[thisPoi];
                      passThreshold = false;
                    } else {
                      cout << " tolerated: " << lastInit << " - " << lastRed << " < " << ThresholdCutList[thisPoi];
                    }
                  } else {
                    passThreshold = false;
                  }
                }
              }
              break;
            }
            thisPosition++;
          }

          tmpIndex++;

          cout << ", position " << thisPosition << " of " << tmpSum
               << " (" << Form("%09f%%", thisReduction*100) << ")" << endl;
        }
        // cout << endl;

        // All individual changes of the uncertainties on the POIs to combine are below
        // the specified thresholds. Accept the found parameter and leave loop
        if (!passThreshold) {
          cout << "PruneNuisanceParameters()  └-> Vetoed parameter because uncertainty on one POI reduced too much!" << endl;
        } else {
          cout << "PruneNuisanceParameters()  └-> Removing parameter!" << endl;
          foundPar2Rem = true;
          break;
        }
      }

      // If no parameter is found, break the pruning procedure
      if (!foundPar2Rem) {
        cout << "PruneNuisanceParameters()  └-> No more parameters found to prune!" << endl;
        break;
      }

      // Add the found parameter to the global ranking and remove it from
      // the list of parameters present in the covariance matrix
      order.push_back(par2rem);
      int par2rem_index = pars.index(par2rem.c_str());
      pars.remove(*pars.at(par2rem_index));
      uncertsPmO.insert(make_pair(par2rem_err, par2rem));
    }
  }

  // Print final ranking of non-filtered parameters based on the uncertainty on combined inclusive POI
  cout << "PruneNuisanceParameters() Ranking of all non-filtered parameters based on combined inclusive POI" << endl;
  PrintRanking(uncertsPmO, initTotalError);
  cout << endl;

  double uncert((--uncertsPmO.rend())->first);
  double delta = uncert/initTotalError - 1;
  cout << "PruneNuisanceParameters() Removing " << uncertsPmO.size() << " of " << allParams - numPoi
       << " parameters changes the uncertainty approximately by "
       << Form("%09f%% (%05f -> %05f)", delta*100, initTotalError, uncert) << endl;

  // Return the set of parameters that can be pruned safely
  return order;
}

// ____________________________________________________________________________|__________
// Remove a set of parameters from a matrix and thus reduce it
void RemoveParameter( TMatrixDSym& hes,
                      RooArgList& pars,
                      list<string> names )
{
  // Find rows and columns to keep and remove
  set<int> removeRows;
  vector<int> keepRows;

  for (list<string>::iterator it = names.begin(); it != names.end(); ++it) {
    string name = *it;
    int index = pars.index(name.c_str());
    removeRows.insert(index);
  }

  RooArgList redpars(pars);
  for (RooLinkedListIter it = redpars.iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    string name = v->GetName();
    int index = redpars.index(name.c_str());
    if (removeRows.find(index) == removeRows.end()) {
      keepRows.push_back(index);
    } else {
      int remindex = pars.index(name.c_str());
      pars.remove(*pars.at(remindex));
    }
  }

  // Remove specified rows and columns from Hesse matrix
  int* a = &keepRows[0];
  TArrayI keepRow(keepRows.size(), a);

  for (Int_t i = 0; i < keepRow.GetSize(); i++) {
    TMatrixDColumn(hes,i) = TMatrixDColumn(hes,keepRow[i]);
    TMatrixDRow(hes,i) = TMatrixDRow(hes,keepRow[i]);
  }

  hes.ResizeTo(keepRow.GetSize(), keepRow.GetSize());
}

// ____________________________________________________________________________|__________
// Print ranking, which is stored in a (ordered) set
void PrintRanking( set< pair< double, string > > uncerts,
                   double initTotalError )
{
  for (set< pair< double, string> >::reverse_iterator itr = uncerts.rbegin(); itr != uncerts.rend(); ++itr) {
    string varName(itr->second);
    double uncert(itr->first);

    double delta = uncert/initTotalError - 1;

    cout << "PrintRanking() "
         << Form("%09f%% (%05f -> %05f): %s",
                 delta*100, initTotalError, uncert, varName.c_str()) << endl;
  }
  cout << endl;
}

// ____________________________________________________________________________|__________
// Given a value and an error, round and format them according to the PDG rules
// for significant digits
pair< double, double > PDGrounding( double value,
                                    double error,
                                    int additionalDigit )
{
  int threeDigits = GetThreeDigits(error);
  int nSignificantDigits = GetNSigDigits(threeDigits);
  nSignificantDigits += additionalDigit;

  // extraRound is meant for the special case of threeDigits > 950
  int extraRound;
  if (threeDigits >= 950) extraRound = 1;
  else extraRound = 0;

  // Convert to mantissa + exponent representation
  int expVal, expErr;
  frexp10(value, &expVal);
  frexp10(error, &expErr);

  // Format the value and error
  double formVal = FormatValue(value, expVal, expVal-expErr+nSignificantDigits, extraRound);
  double formErr = FormatValue(error, expErr, nSignificantDigits, extraRound);

  return make_pair(formVal, formErr);
}

// ____________________________________________________________________________|__________
// Get three digits
int GetThreeDigits( double error )
{
  // Extract the three most significant digits and return them as an integer
  ostringstream stream;
  stream << Form("%.2e", error);
  TString str(stream.str());
  str = ((TObjString*)(str.Tokenize("e"))->At(0))->GetString();
  str.ReplaceAll(".", "");
  str.ReplaceAll("+", "");
  str.ReplaceAll("-", "");

  int threeDigits = atoi(str.Data());

  return threeDigits;
}

// ____________________________________________________________________________|__________
// Get the number of significant digits
int GetNSigDigits( int threeDigits )
{
  // Find the number of significant digits
  assert(threeDigits < 1000);
  int nSignificantDigits;
  if (threeDigits < 101) nSignificantDigits = 2;
  else if (threeDigits < 356) nSignificantDigits = 2;
  else if (threeDigits < 950) nSignificantDigits = 1;
  else nSignificantDigits = 2;

  return nSignificantDigits;
}

// ____________________________________________________________________________|__________
// Convert a number to mantissa + exponent representation in base 10
double frexp10( double x,
                int* exp )
{
  double mantissa = .0 > x ? - x : x;
  *exp = 0;

  if (mantissa >= 10.) {
    *exp = 1;
    for (; !((mantissa /= 10.) < 10.); ++(*exp));
  } else if (!(mantissa > 1.)) {
    *exp = -1;
    for (; !((mantissa *= 10.) > 1.); --(*exp));
  }

  return mantissa;
}

// ____________________________________________________________________________|__________
// Format a value correctly and remove not needed digits
double FormatValue( double value,
                    int exponent,
                    int nDigits,
                    int extraRound )
{
  int roundAt = nDigits - 1 - exponent - extraRound;
  int nDec;
  if (exponent < nDigits) nDec = roundAt;
  else nDec = 0;

  ostringstream stream;
  stream << "%." << nDec << "f";

  double tmp = pow(10, roundAt);
  double formVal = atof(Form(stream.str().c_str(), round(value * tmp) / tmp));

  return formVal;
}

// ____________________________________________________________________________|__________
// Compare two floating point numbers, see http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
union MyFloat_t
{
  MyFloat_t(float num = 0.0f) : f(num) {}
  // Portable extraction of components.
  bool Negative() const { return (i >> 31) != 0; }
  Int_t RawMantissa() const { return i & ((1 << 23) - 1); }
  Int_t RawExponent() const { return (i >> 23) & 0xFF; }

  Int_t i;
  float f;
};

bool AlmostEqualUlpsAndAbs( float A,
                            float B,
                            float maxDiff,
                            int maxUlpsDiff )
{
  // Check if the numbers are really close -- needed  when comparing numbers near zero.
  float absDiff = fabs(A - B);
  if (absDiff <= maxDiff)
    return true;

  MyFloat_t uA(A);
  MyFloat_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative())
    return false;

  // Find the difference in ULPs.
  int ulpsDiff = abs(uA.i - uB.i);
  if (ulpsDiff <= maxUlpsDiff)
    return true;

  return false;
}
