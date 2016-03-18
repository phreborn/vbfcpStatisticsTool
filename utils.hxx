#ifndef _UTILS_
#define _UTILS_

#include "TROOT.h"
#include "TTime.h"
#include "TSystem.h"
#include "TMath.h"

#include "log.hxx"

using namespace std;

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
void PrintResourcesUsed(const TTime& progStart)
{
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
vector<string> parseString(string str, string sep)
{
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

#endif /* _UTILS_ */
