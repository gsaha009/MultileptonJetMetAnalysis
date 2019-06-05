#include <iostream>
#include <memory>
#include "MVAnalysis.h"

using std::string;
using std::cout;
using std::endl;

MVAnalysis::MVAnalysis(const string& mva_algo, const string& xmlfile) 
{
  reader_ = std::make_unique<TMVA::Reader>("!Silent");

  //  reader_->AddVariable("nJet",        &varList_.nJet);
  reader_->AddVariable("met",        &varList_.met);
  reader_->AddVariable("hT",             &varList_.hT);
  reader_->AddVariable("hTvec",             &varList_.hTvec);
  reader_->AddVariable("lepDR",       &varList_.lepDR);
  reader_->AddVariable("lepSumPt",    &varList_.lepSumPt);
  // reader_->AddVariable("lT",    &varList_.lT);
  reader_->AddVariable("MetL1dPhi",   &varList_.MetL1dPhi);
  reader_->AddVariable("MetL2dPhi",   &varList_.MetL2dPhi);
  reader_->AddVariable("MetJ1dPhi",   &varList_.MetJ1dPhi);
  reader_->AddVariable("MetJ2dPhi",   &varList_.MetJ2dPhi);
  //  reader_->AddVariable("lep1eta := abs(lep1Eta)",     &varList_.lep1eta);
  //  reader_->AddVariable("lep2eta := abs(lep2Eta)",     &varList_.lep2eta);
  reader_->AddVariable("j1l1dR", &varList_.j1l1dR);
  //  reader_->AddVariable("j2l1dR", &varList_.j2l1dR);
  reader_->AddVariable("j1j2dR", &varList_.j1j2dR);
  reader_->AddVariable("l1l2InvM",  &varList_.l1l2InvM);
  //  reader_->AddVariable("METsqrtST", &varList_.METsqrtST);
  //reader_->AddVariable("jetDPhi", &varList_.jetDPhi);
  //reader_->AddVariable("TrMass1", &varList_.TrMass1);
  //reader_->AddVariable("mt2 := TMath::Sqrt(2*lep2Pt*met*(1-TMath::Cos(MetL2dPhi)))", &varList_.mt2);
  //  reader_->AddSpectator ("nJet", &varList_.nJet);
  //  reader_->AddVariable("DPhiOverMTsum := (1/(TrMass1+TrMass2))*abs(MetJ1dPhi-MetJ2dPhi)", &varList_.DPhiOverMTsum);
  //  reader_->AddVariable("l1ptOhT := lep1Pt/hT", &varList_.l1ptOhT);

  reader_->BookMVA(mva_algo.c_str(), xmlfile);
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here
  return reader_->EvaluateMVA(mva_algo.c_str());
}
