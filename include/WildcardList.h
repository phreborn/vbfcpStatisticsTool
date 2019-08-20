//=====================================================================-*-C++-*-
// File: $Id: WildcardList.h 549275 2013-05-29 18:26:11Z adye $
//==============================================================================

#ifndef WildcardList_h
#define WildcardList_h

#include <vector>

#include "TString.h"
#include "TRegexp.h"
#include "TObjArray.h"
#include "TObjString.h"

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

#endif
