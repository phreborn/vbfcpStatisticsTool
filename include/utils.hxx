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
#include "TRegexp.h"

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
void save(string baseName, vector<string> type, TCanvas* c1);

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




struct TOwnedList : public TList {
  // A collection class for keeping TObjects for deletion.
  // TOwnedList is like TList with SetOwner(), but really deletes all objects, whether or not on heap.
  // This is a horrible hack to work round the fact that RooArgSet and RooDataSet objects have have IsOnHeap() false.
  TOwnedList();
  virtual ~TOwnedList();
  virtual void Clear (Option_t* option="");
  ClassDef(TOwnedList,0)
};




class WildcardList {
public:
  enum MatchMode { kExact, kWildcard, kRegexp };

private:
  // Use vector<void*> as a type that CINT knows about
  std::vector<void*> _regexp, _exact;

public:

  WildcardList () {}

  WildcardList (const TString& match, MatchMode mode= kWildcard, const char* delim=" \t\n\r\f")
  {
    AddMatches (match, mode, delim);
  }

  WildcardList (const char* const matches[], size_t nmatch, MatchMode mode= kWildcard)
  {
    AddMatches (matches, nmatch, mode);
  }

  Int_t AddMatches (const TString& match, MatchMode mode= kWildcard, const char* delim=" \t\n\r\f")
  {
    TObjArray* a= match.Tokenize(delim);
    if (!a) return 0;
    Int_t nmatch= a->GetEntries();
    size_t n= _exact.size();
    _exact .resize(n+nmatch);
    _regexp.resize(n+nmatch);
    Int_t i;
    for (i= 0; i<nmatch; i++) {
      TObjString* os= (TObjString*)(a->At(i));
      TString* s= new TString (os->GetName());
      _exact [n+i]= s;
      _regexp[n+i]= (mode==kExact) ? 0 : new TRegexp (*s, mode==kWildcard);
    }
    delete a;
    return nmatch;
  }

  void AddMatches (const char* const matches[], size_t nmatch, MatchMode mode= kWildcard)
  {
    size_t n= _exact.size();
    _exact .resize(n+nmatch);
    _regexp.resize(n+nmatch);
    size_t i;
    for (i= 0; i<nmatch; i++) {
      TString* s= new TString (matches[i]);
      _exact [n+i]= s;
      _regexp[n+i]= (mode==kExact) ? 0 : new TRegexp (*s, mode==kWildcard);
    }
  }

  void Add (const char* match, MatchMode mode= kWildcard)
  {
    TString* s= new TString (match);
    _exact .push_back(new TString (match));
    _regexp.push_back ((mode==kExact) ? 0 : new TRegexp (*s, mode==kWildcard));
  }

  Int_t NumTerms() const { return _exact.size(); }

  const TString* Match (const TString& str, Int_t ind=-1) const
  {
    // Matches the string against any of the saved patterns, or just one if ind is specified.
    // Returns the first matching pattern, or else 0.
    Int_t i=0, n=_exact.size();
    if        (ind>=n) {
      return 0;
    } else if (ind>=0) {
      i=ind; n=ind+1;
    }
    for (; i<n; i++) {
      if (_regexp[i]) {
        const TRegexp* r= (const TRegexp*)_regexp[i];
        if (str.Contains(*r)) return (const TString*)_exact[i];
      } else {
        const TString* s= (const TString*)_exact [i];
        if (*s==str) return s;
      }
    }
    return 0;
  }

  Int_t MatchInd (const TString& str) const
  {
    // Matches the string against any of the saved patterns.
    // Returns the first matching pattern number, or else -1.
    for (size_t i=0, n=_exact.size(); i<n; i++) {
      if (_regexp[i]) {
        if (str.Contains(*(const TRegexp*)_regexp[i])) return i;
      } else {
        if (*(const TString*)_exact[i]==str)           return i;
      }
    }
    return -1;
  }

  const TString* Pattern (Int_t i) const
  {
    // Return pattern by its index
    return (i>=0 && size_t(i)<_exact.size()) ? (const TString*)_exact[i] : 0;
  }

  ~WildcardList() {
    size_t i, n;
    for (i= 0, n= _exact .size(); i<n; i++) delete (const TString*)_exact [i];
    for (i= 0, n= _regexp.size(); i<n; i++) delete (const TRegexp*)_regexp[i];
  }

private:

  WildcardList (const WildcardList&);               // copy not implemented
  WildcardList& operator= (const WildcardList&);    // copy not implemented
};



#endif /* _UTILS_ */
