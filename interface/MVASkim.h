#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  //  int nLep;
  float puevwt;
  float nJet;
  float met;
  float jet1Pt;
  float jet2Pt;
  float jet1Eta;
  float jet2Eta;
  float hT;
  float hTvec;


  float lep1Pt;
  float lep1Eta;
  float lep2Pt;
  float lep2Eta;
  float lepSumCharge;  
  float lepDR;
  float lepSumPt;
  float lT;
  float sT;
  float METsqrtST;

  float MetPhi;
  float MetL1dPhi;
  float MetL2dPhi;
  float MetJ1dPhi;
  float MetJ2dPhi;
  float TrMass1;
  float TrMass2;
  float j1l1dR;
  float j1l2dR;
  float j2l1dR;
  float j2l2dR;
  float j1j2dR;
  float l1l2InvM;
  float lepDPhi;
  float jetDPhi;
  float mR;
 
} TreeVariables;

class MVASkim {
    
public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim();

  void fill(const TreeVariables& varList);
  void close();

  TFile* _mvaFile;
  TTree* _tree;

  TreeVariables _varList;
};
#endif
