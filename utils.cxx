// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Common helper functions

#include <chrono>

#include "TROOT.h"
#include "TTime.h"
#include "TSystem.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"

#include "RooProdPdf.h"
#include "RooArgSet.h"

#include "utils.hxx"
#include "log.hxx"

#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"

using namespace std;



TOwnedList::TOwnedList() : TList() { SetOwner(); }
TOwnedList::~TOwnedList()          { Clear(); }
void TOwnedList::Clear(Option_t *option)
{
  if (!option || strcmp(option,"nodelete")!=0)
    for (TIter it(this); TObject* obj= it();) {
      // cout << "Delete "<<obj->ClassName()<<"::"<<obj->GetName()<<(obj->IsOnHeap()?"":" (not on heap)")<<endl;
      delete obj;
    }
  TList::Clear("nodelete");
}


// _____________________________________________________________________________
// Load custom functions
extern void loadCustom() {
  string path = "custom";
  vector<string> ext = {".cxx"};

  LOG(logINFO) << "Load custom classes from " << path;

  boost::filesystem::path targetDir(path);
  boost::filesystem::directory_iterator it(targetDir), eod;

  BOOST_FOREACH(boost::filesystem::path const &fname, std::make_pair(it, eod)) {
    if (boost::filesystem::is_regular_file(fname) && std::find(ext.begin(), ext.end(), it->path().extension()) != ext.end()) {
      gSystem->CompileMacro(fname.c_str(), "k");
    }
  }
}

// _____________________________________________________________________________
// Compare two floating point numbers, see
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
// for original code
bool AlmostEqualUlpsAndAbs(float A, float B, float maxDiff, int maxUlpsDiff) {
  // Check if the numbers are really close -- needed  when comparing numbers
  // near zero.
  float absDiff = fabs(A - B);
  if (absDiff <= maxDiff) return true;

  MyFloat_t uA(A);
  MyFloat_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative()) return false;

  // Find the difference in ULPs.
  int ulpsDiff = abs(uA.i - uB.i);
  if (ulpsDiff <= maxUlpsDiff) return true;

  return false;
}

// _____________________________________________________________________________
// Print used resources
// Courtesy of Tim Adye <T.J.Adye@rl.ac.uk>.
void PrintResourcesUsed(const TTime& progStart) {
  ProcInfo_t info;
  if (gSystem->GetProcInfo(&info)<0) return;
  Long_t cput= TMath::CeilNint(info.fCpuUser);
  Long_t wall= Long64_t(gSystem->Now()-progStart+TTime(500))/Long64_t(1000);
  LOG(logINFO) << Form("resources used: cput=%02ld:%02ld:%02ld, mem=%ldkb, vmem=%ldkb, walltime=%02ld:%02ld:%02ld",
                       cput/3600, (cput/60)%60, cput%60,
                       info.fMemResident, info.fMemVirtual,
                       wall/3600, (wall/60)%60, wall%60);
}

// _____________________________________________________________________________
// Split strings according to separator
vector<string> parseString(string str, string sep) {
  vector<string> parsed;
  int pos = 0;
  bool first = true;
  if (str.size() == 0) return parsed;
  if (str.find(sep) == string::npos) {
    parsed.push_back(str);
    return parsed;
  }

  while (true) {
    int newPos = str.find(sep, pos);
    if (str.find(sep, pos) == string::npos) {
      if (!first) parsed.push_back(str.substr(pos, newPos-pos));
      break;
    }

    string sub = str.substr(pos, newPos-pos);
    parsed.push_back(sub);
    pos = newPos+1;
    first = false;
  }

  return parsed;
}

// _____________________________________________________________________________
// Convert vector to array
double* getAry(vector<double> numbers) {
  int nrPoints = numbers.size();
  double* ary = new double[nrPoints];
  for (int i=0; i<nrPoints; i++) {
    ary[i] = numbers[i];
  }
  return ary;
}

// _____________________________________________________________________________
// Convert deque to array
double* getAry(deque<double> numbers) {
  int nrPoints = numbers.size();
  double* ary = new double[nrPoints];
  for (int i=0; i<nrPoints; i++) {
    ary[i] = numbers[i];
  }
  return ary;
}

// _____________________________________________________________________________
// Split a RooProdPdf into its components
void FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components ) {
  static int counter = 0;
  counter++;

  if (counter > 50) {
    LOG(logERROR) << "FindUniqueProdComponents detected infinite loop. Please check.";
    exit(1);
  }

  RooArgList pdfList = Pdf->pdfList();
  if (pdfList.getSize() == 1) {
    LOG(logINFO) << "FindUniqueProdComponents " << pdfList.at(0)->GetName() << " is fundamental.";
    Components.add(pdfList);
  } else {
    TIterator* pdfItr = pdfList.createIterator();
    RooAbsArg* nextArg;
    while ((nextArg = (RooAbsArg*)pdfItr->Next())) {
      RooProdPdf* Pdf = (RooProdPdf*)nextArg;
      if (string(Pdf->ClassName()) != "RooProdPdf") {
        LOG(logINFO) << "FindUniqueProdComponents " << Pdf->GetName() << " is no RooProdPdf. Adding it.";
        Components.add(*Pdf);
        continue;
      }
      FindUniqueProdComponents(Pdf, Components);
    }
    delete pdfItr;
  }
  counter = 0;
}

// _____________________________________________________________________________
// Return a TGraph with the points of intersection (taken from https://root.cern.ch/phpBB3/viewtopic.php?t=12048)
TGraph* findIntersection(TGraph &a, TGraph &b) {
  TGraph *interPoint = new TGraph();
  int i = 0;

  // Loop over all points in this TGraph
  for(int a_i = 0; a_i < a.GetN()-1; ++a_i) {
    // Loop over all points in the other TGraph
    for(int b_i = 0; b_i < b.GetN()-1; ++b_i) {

      // Get the current point, and the next point for each of the objects
      double x1, y1, x2, y2 = 0;
      double ax1, ay1, ax2, ay2 = 0;
      a.GetPoint(a_i, x1, y1);
      a.GetPoint(a_i+1, x2, y2);
      b.GetPoint(b_i, ax1, ay1);
      b.GetPoint(b_i+1, ax2, ay2);

      // Calculate the intersection between two straight lines, x axis
      double x = (ax1 *(ay2 *(x1-x2)+x2 * y1 - x1 * y2 )+ ax2 * (ay1 * (-x1+x2)- x2 * y1+x1 * y2)) / (-(ay1-ay2) * (x1-x2)+(ax1-ax2)* (y1-y2));

      // Calculate the intersection between two straight lines, y axis
      double y = (ax1 * ay2 * (y1-y2)+ax2 * ay1 * (-y1+y2)+(ay1-ay2) * (x2 * y1-x1 * y2))/(-(ay1-ay2) * (x1-x2)+(ax1-ax2) * (y1-y2));

      // Find the tightest interval along the x-axis defined by the four points
      double xrange_min = max(min(x1, x2), min(ax1, ax2));
      double xrange_max = min(max(x1, x2), max(ax1, ax2));

      if ((x1 == ax1 and y1 == ay1)or (x2 == ax2 and y2 == ay2)) {
        // If points from the two lines overlap, they are trivially intersecting
        interPoint->SetPoint(i, (x1 == ax1 and y1 == ay1) ? x1 : x2, (x1 == ax1 and y1 == ay1) ? y1 : y2);
        i++;
      } else if(x > xrange_min && x < xrange_max) {
        // If the intersection between the two lines is within the tight range, add it to the list of intersections.
        interPoint->SetPoint(i,x, y);
        i++;
      }
    }
  }

  return interPoint;
}

// _____________________________________________________________________________
TGraph* makeGraph(string title, int n, double* x_ary, double* y_ary) {
  TGraph* graph = new TGraph(n, x_ary, y_ary);
  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("X");
  graph->GetYaxis()->SetTitle(title.c_str());
  return graph;
}

// _____________________________________________________________________________
TGraphAsymmErrors* makeGraphErr(string title, int n, double* x_ary, double* central, double* errlo, double* errhi) {
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(n, x_ary, central, NULL, NULL, errlo, errhi);
  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("X");
  graph->GetYaxis()->SetTitle(title.c_str());
  return graph;
}

// _____________________________________________________________________________
TGraphAsymmErrors* makeGraphErr(string title, int n, double* x_ary, double* central, double* cenlo, double* cenhi, double* errlo, double* errhi) {
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(n, x_ary, central, cenlo, cenhi, errlo, errhi);
  graph->SetTitle("");
  graph->GetXaxis()->SetTitle("X");
  graph->GetYaxis()->SetTitle(title.c_str());
  return graph;
}

// _____________________________________________________________________________
double subtract_error(double err12, double err1) {
  double q = (err12 * err12 - err1 * err1);
  if (q >= 0.0) {
    return sqrt(err12 * err12 - err1 * err1);
  } else {
    LOG(logWARNING) << "Please check: q=" << q;
    return 1e-09;
  }
}

// _____________________________________________________________________________
// Given a value and an error, round and format them according to the PDG rules
// for significant digits
pair< double, double > PDGrounding( double value, double error, int additionalDigit )
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

// _____________________________________________________________________________
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

// _____________________________________________________________________________
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

// _____________________________________________________________________________
// Convert a number to mantissa + exponent representation in base 10
double frexp10( double x, int* exp )
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

// _____________________________________________________________________________
// Format a value correctly and remove not needed digits
double FormatValue( double value, int exponent, int nDigits, int extraRound )
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
