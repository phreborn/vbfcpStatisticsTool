// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@cern.ch
// Date        : 2016-03-17
// Description : Common helper functions

#ifndef _UTILS_
#define _UTILS_

#include <chrono>

#include "TROOT.h"
#include "TTime.h"
#include "TSystem.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"

#include "RooProdPdf.h"
#include "RooArgSet.h"

#include "log.hxx"

#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"

using namespace std;

// _____________________________________________________________________________
// Defined functions
extern "C" void loadCustom();
bool AlmostEqualUlpsAndAbs(float A, float B, float maxDiff, int maxUlpsDiff);
void PrintResourcesUsed(const TTime& progStart);
vector<string> parseString(string str, string sep);
double* getAry(vector<double> numbers);
double* getAry(deque<double> numbers);
void FindUniqueProdComponents(RooProdPdf* Pdf, RooArgSet& Components);
TGraph* findIntersection(TGraph &a, TGraph &b);
TGraph* makeGraph(string title, int n, double* x_ary, double* y_ary);
TGraphAsymmErrors* makeGraphErr(string title, int n, double* x_ary, double* central, double* errlo, double* errhi);
TGraphAsymmErrors* makeGraphErr(string title, int n, double* x_ary, double* central, double* cenlo, double* cenhi, double* errlo, double* errhi);
double subtract_error(double err12, double err1);
pair< double, double > PDGrounding( double value, double error, int additionalDigit );
int GetThreeDigits( double error );
int GetNSigDigits( int threeDigits );
double frexp10( double x, int* exp );
double FormatValue( double value, int exponent, int nDigits, int extraRound );

// _____________________________________________________________________________
// Timer for benchmarks, see
// http://www.answerandquestion.net/questions/12290/c-timer-function-to-provide-time-in-nano-seconds
// for the original code
class MyTimer
{
public:
  MyTimer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const {
    return chrono::duration_cast<second_>
      (clock_::now() - beg_).count(); }

private:
  typedef chrono::high_resolution_clock clock_;
  typedef chrono::duration<double, ratio<1> > second_;
  chrono::time_point<clock_> beg_;
};

// _____________________________________________________________________________
// Compare two floating point numbers, see
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
// for original code
union MyFloat_t {
  MyFloat_t(float num = 0.0f) : f(num) {}
  // Portable extraction of components.
  bool Negative() const { return (i >> 31) != 0; }
  Int_t RawMantissa() const { return i & ((1 << 23) - 1); }
  Int_t RawExponent() const { return (i >> 23) & 0xFF; }

  Int_t i;
  float f;
};

#endif /* _UTILS_ */
