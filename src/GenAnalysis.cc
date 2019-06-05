#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <functional>
#include <numeric>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <utility> 
#include <typeinfo>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"

#include "AnaUtil.h"
#include "LL4JMETUtil.h"
#include "GenAnaBase.h"
#include "GenAnalysis.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::setw;
using std::vector;
using std::string;

GenAnalysis::GenAnalysis() 
{
  //using GenAnaBase::GenAnaBase;
}
GenAnalysis::~GenAnalysis() {
}

void GenAnalysis::bookHistograms(std::unique_ptr<TFile>& histf) {
  histf->cd();
  histf->mkdir("GenLevel");
  histf->cd("GenLevel");
  new TH1D("nlep", "#leptons from H++/H--", 10, -0.5, 9.5);
  new TH1D("nGenLep", "", 5, -0.5, 4.5);
  new TH1D("nq", "#quarks from H++/H-- & H+/H-", 10, -0.5, 9.5);
  new TH1D("nNeu", "#neutrinos from H++/H-- & H+/H-", 10, -0.5, 9.5);
  new TH1D("lepPt", "lep pt from H++/H--", 600, 0., 300.);
  new TH1D("qPt", "quark pt from H++/H-- & H+/H-", 600, 0., 300.);
  new TH1D("qP", "quark pt from H++/H-- & H+/H-", 800, 0., 400.);  
  new TH1D("qEta", "quark pt from H++/H-- & H+/H-", 20, -5., 5.);
  new TH1D("qPhi", "quark pt from H++/H-- & H+/H-", 16, -4., 4.);
  new TH1D("qE", "quark pt from H++/H-- & H+/H-", 800, 0., 400.);
  
  new TH1D("neutrinoPt", "neutrino pt from H++/H--", 600, 0., 300.);
  //new TH1D("Pt", "quark pt from H++/H-- & H+/H-", 400, 0., 200.);
  new TH1D("hplepPt", "lep pt from H++", 400, 0., 200.);
  new TH1D("hmlepPt", "lep pt from H--", 400, 0., 200.);
  new TH1D("neuEt", "Et", 400, 0., 200.);
  new TH1D("hpmLL", "H++ mass", 400, 0., 400.);
  new TH1D("hmmLL", "H-- mass", 400, 0., 400.);
  new TH1D("MET", "", 400, 0., 400.);
  new TH1D("neuPt", "", 400, 0., 400.);
  new TH1D("hpmcombined", "", 320, 70., 230.);
  new TH1D("wM", "", 500, 198., 202.);
  new TH1D("mhAll", "", 1000, 0., 250.);
  new TH1D("lepTrkEta", "", 100, -5.0, 5.0);
  new TH1D("lepEta", "", 100, -5.0, 5.0);
  histf->cd();
}

void GenAnalysis::clearLists (){
  leptonBox.clear();
  quarkBox.clear();
  neutrinoBox.clear();
  hP_Decay.clear();
  hM_Decay.clear();
}

void GenAnalysis::analyze(const std::unique_ptr<TFile>& histf) {
  histf->cd();
  histf->cd("GenLevel");
  genPass_ = false;
  //  vector <int> dIndices;
  /*leptonBox.clear();
  hP_Decay.clear();
  hM_Decay.clear();*/
  TLorentzVector hplepP4;
  hplepP4.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector hmlepP4;
  hmlepP4.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector hneuP4;
  hneuP4.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector hpcombinedP4;
  hpcombinedP4.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector wpP4;
  wpP4.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector wmP4;
  wmP4.SetPtEtaPhiE(0., 0., 0., 0.);
  TLorentzVector qP4;
  qP4.SetPtEtaPhiE(0., 0., 0., 0.);
  clearLists();  
  
  //std::cout<<"*******************************"<<std::endl;
  //Entering into genparticleCollection
  for (const auto& gp: *genPList_) {
    int pdgid = gp.pdgId;
    int abs_pid = std::abs(gp.pdgId);
    int status = gp.status;
  
    if (status == 1) { // T R E A T I N G  L E P T O N S
      if (abs_pid != 11 && abs_pid != 13 && abs_pid != 9000012 && abs_pid != 9000014 && abs_pid != 9000016) continue;
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      if (index < 0) continue;
      const vhtm::GenParticle& gp_mom = genPList_ -> at(index);
      int lepMomId = gp_mom.pdgId;
      //std::cout<<lepMomId<<std::endl;
      if (std::abs(lepMomId) != 24) continue;
      int index_gm = getMotherId(gp_mom, mmid);
      if (index_gm < 0) continue;
      const vhtm::GenParticle& gp_gmom = genPList_ -> at(index_gm);
      int lepgMomId = gp_gmom.pdgId;
      //std::cout<<"pid: "<<pdgid<<"\t"<<"mom: "<<lepMomId<<"\t"<<"gmom: "<<lepgMomId<<std::endl;

      if (std::abs(lepgMomId) != 9900041) continue;
      if (abs_pid == 11 || abs_pid == 13) {
	leptonBox.push_back(gp);
	AnaUtil::fillHist1D("lepPt", gp.pt);
	if (std::abs(gp.eta) < 2.5) AnaUtil::fillHist1D("lepTrkEta", gp.eta);
	else AnaUtil::fillHist1D("lepEta", gp.eta);
      }
      if (abs_pid == 9000012 || abs_pid == 9000014 || abs_pid == 9000016){
	neutrinoBox.push_back(gp);
	AnaUtil::fillHist1D("neutrinoPt", gp.pt);

      }

      if (lepgMomId == 9900041) {
	hP_Decay.push_back(gp);
	if (pdgid == -11 || pdgid == -13) {
	  wpP4 += (LL4JMETUtil::getP4(gp_mom));
	}    
      }
      else if (lepgMomId == -9900041) {
	hM_Decay.push_back(gp);
	if (pdgid == 11 || pdgid == 13) { 
	  wmP4 += (LL4JMETUtil::getP4(gp_mom));
	}
      }
    }
    //    std::cout<<"jet*****************jet"<<std::endl;    
    if (status == 71) { // T R E A T I N G  H A D R O N S
      if (abs_pid != 1 && abs_pid != 2 && abs_pid != 3 && abs_pid != 4 && abs_pid != 5 && abs_pid != 6) continue;
      int mmid = -1;
      int index = getMotherIdForQ(gp, mmid);
      if (index < 0) continue;
      const vhtm::GenParticle& gp_mom = genPList_ -> at(index);
      int qMomId = gp_mom.pdgId;
      //      std::cout<<"jpid: "<<pdgid<<"\t"<<"jmom: "<<qMomId<<"\t"<<"pt: "<<gp.pt<<std::endl;

      if (std::abs(qMomId) != 24 && std::abs(qMomId) != 23) continue;
      int index_gm = getMotherId(gp_mom, mmid);
      if (index_gm < 0) continue;
      const vhtm::GenParticle& gp_gmom = genPList_ -> at(index_gm);
      int qgMomId = gp_gmom.pdgId;
      //std::cout<<"jpid: "<<pdgid<<"\t"<<"jmom: "<<qMomId<<"\t"<<"jgmom: "<<qgMomId<<"\t"<<"pt: "<<gp.pt<<"\t"<<"eta: "<<gp.eta<<"\t"<<"phi: "<<gp.phi<<std::endl;
      if (std::abs(qgMomId) == 9900041 || std::abs(qgMomId) == 37){
	quarkBox.push_back(gp);
	AnaUtil::fillHist1D("qPt", gp.pt);
	AnaUtil::fillHist1D("qEta", gp.eta);
	AnaUtil::fillHist1D("qPhi", gp.phi);
	AnaUtil::fillHist1D("qE", gp.energy);
	AnaUtil::fillHist1D("qP", gp.p);
      }
    }
  } //Loop finished

  AnaUtil::fillHist1D("nlep", leptonBox.size());
  AnaUtil::fillHist1D("nq", quarkBox.size());
  // if (quarkBox.size() == 2) dumpEvent();
  //std::cout<<"__________________________________________________________________"<<std::endl;
  if (quarkBox.size() == 4) {
    for (auto& ob: quarkBox){
      //std::cout<<ob.pdgId<<"\t"<<ob.status<<"\t"<<ob.pt<<"\t"<<ob.p<<"\t"<<ob.energy<<std::endl;
      //      std::cout<<ob.pdgId<<"\t"<<ob.eta<<std::endl;
      qP4 += (LL4JMETUtil::getP4(ob));
    }
  }
  //  std::cout<<qP4.M()<<std::endl;
  AnaUtil::fillHist1D("mhAll", qP4.M());
  AnaUtil::fillHist1D("nNeu", neutrinoBox.size());

  bool size4 = false;
  if (hP_Decay.size() == 4) size4 = true;
  for (auto &a: hP_Decay) {
    if (size4) hpcombinedP4 += (LL4JMETUtil::getP4(a));
    if (std::fabs(a.pdgId) == 11 || std::fabs(a.pdgId) == 13){
      AnaUtil::fillHist1D("hplepPt", a.pt);
      hplepP4 += (LL4JMETUtil::getP4(a));
    }
    else if (std::fabs(a.pdgId) == 9000012 || std::fabs(a.pdgId) == 9000014 || std::fabs(a.pdgId) == 9000016){
      AnaUtil::fillHist1D("neuEt", a.pt);
      hneuP4 += (LL4JMETUtil::getP4(a));
    }
  }  
  for (auto &a: hM_Decay) {
    if (std::fabs(a.pdgId) == 11 || std::fabs(a.pdgId) == 13){
      AnaUtil::fillHist1D("hmlepPt", a.pt);
      hmlepP4 += (LL4JMETUtil::getP4(a));
    }
    else if (std::fabs(a.pdgId) == 9000012 || std::fabs(a.pdgId) == 9000014 || std::fabs(a.pdgId) == 9000016){
      AnaUtil::fillHist1D("neuEt", a.pt);
      hneuP4 += (LL4JMETUtil::getP4(a));
    }
  }  
  //hasError = false;
  AnaUtil::fillHist1D("hpmLL", hplepP4.M());
  AnaUtil::fillHist1D("hmmLL", hmlepP4.M());
  AnaUtil::fillHist1D("MET", hneuP4.Et());
  AnaUtil::fillHist1D("neuPt", hneuP4.Pt());
  AnaUtil::fillHist1D("hpmcombined", hpcombinedP4.M());
  AnaUtil::fillHist1D("wM", wpP4.M());
  genPass_ = true;
}

//**************** F I L T E R *******************//
bool GenAnalysis::filter() const {
  int ngenp = genPList_->size();
  if (!ngenp) return false;
  bool hasbJet = false;
  int nJet = 0;
  int nLep = 0;
  double charge_ = 1.;
  for (const auto& gp: *genPList_) {
    int pdgid = std::abs(gp.pdgId);
    int status = gp.status;
    double genEta = std::abs(gp.eta);
    double charge = gp.charge;
    
    if (pdgid == 5 && gp.pt > 20) {
      hasbJet = true;
      break;
    }
    int mmid = -1;
    //    if ((pdgid == 11 || pdgid == 13) && gp.pt >= 10. && status == 1 && genEta < 2.5) {
    if ((pdgid == 11 || pdgid == 13) && status == 1) {
      int index = getMotherId(gp, mmid);
      if (index < 0) continue;
      int lepMomId = std::abs((genPList_->at(index)).pdgId);
      if (lepMomId == 24) {
	nLep++;
	charge_ *= charge; 
      }
    }
#if 0
    if ((pdgid >= 1 && pdgid <= 4) && genEta < 4.7 && gp.pt > 20 && status == 71){
      vector<int> m = gp.motherIndices;
      if (m.size() != 1) continue;
      int index = getMotherIdForQ(gp, mmid);
      if (index < 0) continue;
      int qMomId = std::abs((genPList_->at(index)).pdgId);
      if (qMomId == 24 || qMomId == 23) nJet++;
    }
#endif 
  }
  //  std::cout<<"pid: "<<gp.pdgId<<"\t"<<"mom: "<<
  //  std::cout<<"#jets: "<<"\t"<<nJet<<std::endl;
  //  if (nLep > 1)  std::cout<<"charge  : "<<"\t"<<charge_<<std::endl;
  if (hasbJet) return false;
  AnaUtil::fillHist1D("nGenLep", nLep);
  bool has2Leps = (charge_ > 0 && nLep == 2) ? true : false;
  //bool has4Jets  = (nJet >= 4 ) ? true : false;
  //bool has2Leps4Jets = (has2Leps && has4Jets);
  //if (has2Leps4Jets) return true;
  if (has2Leps) return true;
  return false;
}
