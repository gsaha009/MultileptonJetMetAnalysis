#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct  
{
  //  float nJet;
  float met;
  float hT;
  float hTvec;
  float lepDR;
  float lepSumPt;
  //  float lT;  
  float MetL1dPhi;
  float MetL2dPhi;
  float MetJ1dPhi;
  float MetJ2dPhi;
  //float lep1eta;
  //float lep2eta;
  float j1l1dR;
  //  float j2l1dR;
  float j1j2dR;
  float l1l2InvM;
  //  float METsqrtST;
  //float jetDPhi;
  //float TrMass1;
  //  float MTsum;
  //  float DPhiOverMTsum;
  // float l1ptOhT;
#if 0
  float lepOjetPt;
  float TrMass;
  float jetSumVecPt;  
  float lepOjetPt;

  //  float j1l2dR;
  //
  //  float j2l2dR;

  //  float HJ1dR;
  //  float HJ2dR;

  //float jetsInvM;
  float hToverMET;
  float TrMass1;
  float mt2;
#endif
} InputVariables;

class MVAnalysis {
    
public:

  MVAnalysis(const std::string& mva_algo, const std::string& xmlfile);
  virtual ~MVAnalysis() {}

  double evaluate(const std::string& tag, const InputVariables& varList);

  InputVariables varList_;
  std::unique_ptr<TMVA::Reader> reader_;
};
#endif
