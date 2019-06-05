#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  _mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  _tree = new TTree("RTree", "RTree");
  _tree->Branch("puevwt",        &_varList.puevwt,          "puevwt/F");
  _tree->Branch("nJet",          &_varList.nJet,            "nJet/F");
  _tree->Branch("met",           &_varList.met,             "met/F");
  _tree->Branch("jet1Pt",        &_varList.jet1Pt,          "jet1Pt/F");
  _tree->Branch("jet2Pt",        &_varList.jet2Pt,          "jet2Pt/F");
  _tree->Branch("jet1Eta",       &_varList.jet1Eta,         "jet1Eta/F");
  _tree->Branch("jet2Eta",       &_varList.jet2Eta,         "jet2Eta/F");
  _tree->Branch("hT",            &_varList.hT,              "hT/F");
  _tree->Branch("hTvec",         &_varList.hTvec,           "hTvec/F");

  _tree->Branch("lep1Pt",        &_varList.lep1Pt,          "lep1Pt/F");
  _tree->Branch("lep1Eta",       &_varList.lep1Eta,         "lep1Eta/F");
  _tree->Branch("lep2Pt",        &_varList.lep2Pt,          "lep2Pt/F");
  _tree->Branch("lep2Eta",       &_varList.lep2Eta,         "lep2Eta/F");
  _tree->Branch("lepSumCharge",  &_varList.lepSumCharge,    "lepSumCharge/F");
  _tree->Branch("lepDR",         &_varList.lepDR,           "lepDR/F");
  _tree->Branch("lepSumPt",      &_varList.lepSumPt,        "lepSumPt/F");
  _tree->Branch("lT",            &_varList.lT,              "lT/F");
  _tree->Branch("sT",            &_varList.sT,              "sT/F");
  _tree->Branch("METsqrtST",     &_varList.METsqrtST,       "METsqrtST/F");

  _tree->Branch("MetPhi",        &_varList.MetPhi,          "MetPhi/F");
  _tree->Branch("MetL1dPhi",     &_varList.MetL1dPhi,       "MetL1dPhi/F");
  _tree->Branch("MetL2dPhi",     &_varList.MetL2dPhi,       "MetL2dPhi/F");
  _tree->Branch("MetJ1dPhi",     &_varList.MetJ1dPhi,       "MetJ1dPhi/F");
  _tree->Branch("MetJ2dPhi",     &_varList.MetJ2dPhi,       "MetJ2dPhi/F");
  _tree->Branch("TrMass1",       &_varList.TrMass1,         "TrMass1/F");
  _tree->Branch("TrMass2",       &_varList.TrMass2,         "TrMass2/F");
  _tree->Branch("j1l1dR",        &_varList.j1l1dR,          "j1l1dR/F");
  _tree->Branch("j1l2dR",        &_varList.j1l2dR,          "j1l2dR/F");
  _tree->Branch("j2l1dR",        &_varList.j2l1dR,          "j2l1dR/F");
  _tree->Branch("j2l2dR",        &_varList.j2l2dR,          "j2l2dR/F");
  _tree->Branch("j1j2dR",        &_varList.j1j2dR,          "j1j2dR/F");
  _tree->Branch("l1l2InvM",      &_varList.l1l2InvM,        "l1l2InvM/F");
  _tree->Branch("lepDPhi",        &_varList.lepDPhi,          "lepDPhi/F");
  _tree->Branch("jetDPhi",        &_varList.jetDPhi,          "jetDPhi/F"); 
  _tree->Branch("mR",        &_varList.mR,          "mR/F"); 

  _mvaFile->ls();
}
MVASkim::~MVASkim() {
  if (_tree) delete _tree;  
  if (_mvaFile) delete _mvaFile;
}
void MVASkim::fill(const TreeVariables& varList) {
  memcpy(&_varList, &varList, sizeof(varList));
  _mvaFile->cd();
  _tree->Fill();
}
void MVASkim::close() {
  //_mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  //  _tree->Print();
  _tree->Write();
  _mvaFile->Save();
  _mvaFile->Write();
  //_mvaFile->Close();
}
