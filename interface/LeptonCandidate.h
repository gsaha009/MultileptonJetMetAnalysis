#ifndef __LeptonCandidate__hh
#define __LeptonCandidate__hh

#include "TLorentzVector.h"
struct LeptonCandidate {
  int lIndex;
  TLorentzVector lP4;
  TLorentzVector lFsrP4;
  TLorentzVector lRP4;
  int lCharge;
  double lEta;
  double lPt;
  double lIsolation;
  int flavour; // 1:muon, 2:electron, -1:others
};
#endif
