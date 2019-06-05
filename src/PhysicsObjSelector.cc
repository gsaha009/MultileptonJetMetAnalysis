#include "configana.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <utility> 
#include <typeinfo>

#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TVector2.h"

#include "AnaUtil.h"
#include "LL4JMETUtil.h"
#include "PhysicsObjSelector.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;

PhysicsObjSelector::PhysicsObjSelector()
  : AnaBase()
{
}
bool PhysicsObjSelector::beginJob() {
  if (!AnaBase::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("ObjectSelection");
  return true;
}
void PhysicsObjSelector::endJob() {
  objectEfficiency();
  AnaBase::endJob();
}
bool PhysicsObjSelector::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;
  return true;
}
void PhysicsObjSelector::objectEfficiency() {
  // Object selection Efficiency
  histf()->cd();
  histf()->cd("ObjectSelection");

  // Muon
  vector<string> muLabels {
    "All",
    "pT",
    "eta",
    "dxyPV",
    "dzPV",
    "muType",
    "bestTrack",
    "sip3d",
    "GhostCleaned",
    "PF"
  };
  LL4JMETUtil::showEfficiency("muCutFlow", muLabels, "Muon Selection", "Muons");  

  // Electron
  vector<string> eleLabels {
    "All",
    "pT",
    "eta",
    "dxyPV",
    "dzPV",
    "Missing hits",
    "sip3d",
    "CrossCleaned",
    "BDT"
  };
  LL4JMETUtil::showEfficiency("eleCutFlow", eleLabels, "Electron Selection", "Electrons");  

  // Tau
  vector<string> tauLabels {
    "All",
    "dz",
    "pT",
    "eta",
    "decay mode",
    "isolation",
    "muon veto",
    "electron veto"
  };
  LL4JMETUtil::showEfficiency("tauCutFlow", tauLabels, "Tau Selection", "Taus");  

  // Jet
  vector<string> jetLabels {
    "All",
    "pT",
    "eta",
    "MVAid",
    "jetLeptonCleaning",
    "LooseJet",
    "JetIsTight"
  };
  LL4JMETUtil::showEfficiency("jetCutFlow", jetLabels, "Jet Selection", "Jets");  

}
void PhysicsObjSelector::bookHistograms() {
  histf()->cd();
  histf()->cd("ObjectSelection");
  new TH1D("muCutFlow", "Muon Cut Flow", 10, -0.5, 9.5);
  new TH1D("eleCutFlow", "Electron Cut Flow", 9, -0.5, 8.5);
  new TH1D("tauCutFlow", "Tau Cut Flow", 8, -0.5, 7.5);
  new TH1D("jetCutFlow", "Jet Cut Flow", 7, -0.5, 6.5);
  new TH1D("photonCutFlow", "Photon Cut Flow", 8, -0.5, 7.5);

  new TH1D("nJets", "Number of jets cleaned w.r.t tight leptons passing iso per event", 20, -0.5, 19.5);
  new TH1D("bDisc", "bDiscriminator", 50, 0., 2.);
  new TH1D("nbJets", "Number of b-jets per event", 5, -0.5, 4.5);

  new TH1F("tauPt", "Tau p_{T}", 200, 0., 200.);
  new TH1F("muPt", "Muon p_{T}", 200, 0., 200.);
  new TH1F("elePt", "Electron p_{T}", 200, 0., 200.);
  new TH1F("phoPt", "Photon p_{T}", 200, 0., 200.);
  new TH1F("jetPt", "Jet p_{T}", 100, 0., 500.);
  new TH1F("JetMuDR", "", 100, 0., 10.);
  new TH1F("JetEleDR", "", 100, 0., 10.);
  new TH1D("nJets_0", "", 20, -0.5, 19.5);
  new TH1D("nJets_1", "", 20, -0.5, 19.5);
  new TH1D("nJets_2", "", 20, -0.5, 19.5);
  new TH1D("nJets_3", "", 20, -0.5, 19.5);
  new TH1D("nJets_4", "", 20, -0.5, 19.5);
  new TH1D("nJets_5", "", 20, -0.5, 19.5);
}
// clear lists
void PhysicsObjSelector::clear() {
  preSIPLooseMuList_.clear();
  looseMuList_.clear();
  tightMuList_.clear();
  tightIsoMuList_.clear();

  preSIPLooseEleList_.clear();
  looseEleList_.clear();
  tightEleList_.clear();
  tightIsoEleList_.clear();

  looseJetList_.clear();
  tightJetList_.clear();
  nbJets_ = 0;

  looseMuPhotonPairList_.clear();
  tightMuPhotonPairList_.clear();
  tightIsoMuPhotonPairList_.clear();

  looseElePhotonPairList_.clear();
  tightElePhotonPairList_.clear();
  tightIsoElePhotonPairList_.clear();

  fsrPhotonList_.clear();

  tauList_.clear();
  isoTauList_.clear();

  searchedMu_ = false; 
  searchedEle_ = false; 
  searchedPhoton_ = false;
}
void PhysicsObjSelector::findObjects(double vz, double wt) {
//void PhysicsObjSelector::findObjects(double wt) {
  // order of execution is crucial!
  const vhtm::Event& evt = eventColl()->at(0);
  fGridRhoFastjetAll_ = evt.fGridRhoFastjetAll;

  // muonSelector must precede electronSelector
  muonSelector(wt);
  electronSelector(wt);
  tauSelector(vz, wt);
  //tauSelector(wt);
  
  // after electrons and muons are found, find photons
  photonSelector(wt);

  // Now find isolated leptons after removing contribution from the associated photons
  isoLeptonSelector();

  // and finally the Jets
  jetSelector(wt);
}
// Tau selection
void PhysicsObjSelector::tauSelector(double vz, double wt) {
  //void PhysicsObjSelector::tauSelector(double wt) {
  histf()->cd();
  histf()->cd("ObjectSelection");
  for (const auto& tau: *tauColl()) {
    AnaUtil::fillHist1D("tauCutFlow", 0, wt);
    
    if (std::fabs(tau.dzPV) >= AnaUtil::cutValue(tauCutMap(), "dz")) continue;   
    AnaUtil::fillHist1D("tauCutFlow", 1, wt);
    
    AnaUtil::fillHist1D("tauPt", tau.pt, wt);
    if (tau.pt <= AnaUtil::cutValue(tauCutMap(), "pt")) continue; 
    AnaUtil::fillHist1D("tauCutFlow", 2, wt);

    if (std::fabs(tau.eta) >= AnaUtil::cutValue(tauCutMap(), "eta")) continue;
    AnaUtil::fillHist1D("tauCutFlow", 3, wt);
    
    if (tau.decayModeFinding < 0.99) continue;
    AnaUtil::fillHist1D("tauCutFlow", 4, wt);
    
    if (tau.againstMuonTight3 < AnaUtil::cutValue(tauCutMap(), "muVeto")) continue;                        
    AnaUtil::fillHist1D("tauCutFlow", 5, wt);
    
    if (tau.againstElectronTightMVA <= AnaUtil::cutValue(tauCutMap(), "eleVeto")) continue;                                            
    //    if (tau.againstElectronMediumMVA < AnaUtil::cutValue(tauCutMap(), "eleVeto")) continue;                                            
    AnaUtil::fillHist1D("tauCutFlow", 6, wt);
    
    tauList_.push_back(tau);
  }
  if (tauList_.size() > 1) 
    std::sort(std::begin(tauList_), std::end(tauList_), PtComparator<vhtm::Tau>());

  for (const auto& tau: tauList_) {
    if (tau.byMediumIsolationMVArun2v1DBoldDMwLT < AnaUtil::cutValue(tauCutMap(), "isol")) continue;                     
    AnaUtil::fillHist1D("tauCutFlow", 7, wt);
    
    isoTauList_.push_back(tau);
  }
}
// Muon selection
void PhysicsObjSelector::muonSelector(double wt) {
  histf()->cd();
  histf()->cd("ObjectSelection");

  for (const auto& muon: *muonColl()) {
    AnaUtil::fillHist1D("muCutFlow", 0, wt);

    AnaUtil::fillHist1D("muPt", muon.pt, wt);
    if (muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt"))                             continue;
    AnaUtil::fillHist1D("muCutFlow", 1, wt);

    if (std::fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta"))                continue;
    AnaUtil::fillHist1D("muCutFlow", 2, wt);

    if (std::fabs(muon.dxyPV) >= AnaUtil::cutValue(muonCutMap(), "dxyPV") )           continue;
    AnaUtil::fillHist1D("muCutFlow", 3, wt);

    if (std::fabs(muon.dzPV) >= AnaUtil::cutValue(muonCutMap(), "dzPV") )             continue;
    AnaUtil::fillHist1D("muCutFlow", 4, wt);

    bool muType = muon.isGlobalMuon || (muon.isTrackerMuon && muon.matches > 0);
    if (!muType)                                                                      continue;
    AnaUtil::fillHist1D("muCutFlow", 5, wt);

    if (muon.muonBestTrackType == 2)                                                  continue;
    AnaUtil::fillHist1D("muCutFlow", 6, wt);

    if (!muon.isghostCleaned)                                                         continue;
    AnaUtil::fillHist1D("muCutFlow", 7, wt);

    preSIPLooseMuList_.push_back(muon);
  }
  if (preSIPLooseMuList_.size() > 1) 
    std::sort(std::begin(preSIPLooseMuList_), std::end(preSIPLooseMuList_), PtComparator<vhtm::Muon>());

  for (const auto& muon: preSIPLooseMuList_) {
    //std::cout<<"sip3D: "<<std::fabs(muon.dB3D/muon.edB3D)<<std::endl;
    if (std::fabs(muon.dB3D/muon.edB3D) >= AnaUtil::cutValue(muonCutMap(), "SIP3D"))  continue;
    AnaUtil::fillHist1D("muCutFlow", 8, wt);

    // loose muon 
    looseMuList_.push_back(muon);

    // prepare a vector<pair<vhtm::Muon, vector<vhtm::PackedPFCandidate>>> for each loose muon
    // attaching an empty photon vector
    std::vector<vhtm::PackedPFCandidate> phov;
    looseMuPhotonPairList_.push_back({muon, phov});

    // tight muon 
    bool isTight = muon.isPFMuon || (muon.passTrackerhighPtid && muon.pt > 200.);
    if (!isTight) continue;
    AnaUtil::fillHist1D("muCutFlow", 9, wt);
    tightMuList_.push_back(muon);
  }
  
  searchedMu_ = true;
}
// Electron selecton
void PhysicsObjSelector::electronSelector(double wt) {
  histf()->cd();
  histf()->cd("ObjectSelection");

  if (!searchedMu_) muonSelector();

  for (const auto& electron: *electronColl()) {
    AnaUtil::fillHist1D("eleCutFlow", 0, wt);
    
    AnaUtil::fillHist1D("elePt", electron.pt, wt);
    if (electron.pt <= AnaUtil::cutValue(electronCutMap(), "pt"))                                continue;
    AnaUtil::fillHist1D("eleCutFlow", 1, wt);

    if (std::fabs(electron.eta) >= AnaUtil::cutValue(electronCutMap(), "eta"))                   continue;
    AnaUtil::fillHist1D("eleCutFlow", 2, wt);

    if (std::fabs(electron.dxyPV) >= AnaUtil::cutValue(electronCutMap(), "dxyPV"))               continue;
    AnaUtil::fillHist1D("eleCutFlow", 3, wt);

    if (std::fabs(electron.dzPV) >= AnaUtil::cutValue(electronCutMap(), "dzPV"))                 continue;
    AnaUtil::fillHist1D("eleCutFlow", 4, wt);

    // we don't use the missingHits cut anymore since this variable is now part of the BDT inputs. 
    //if (electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits"))               continue;
    AnaUtil::fillHist1D("eleCutFlow", 5, wt);

    preSIPLooseEleList_.push_back(electron);
  }
  if (preSIPLooseEleList_.size() > 1)
    std::sort(std::begin(preSIPLooseEleList_), std::end(preSIPLooseEleList_), PtComparator<vhtm::Electron>());

  for (const auto& electron: preSIPLooseEleList_) {
    // Cross cleaning
    if (!crossCleaned(electron))  continue;
    AnaUtil::fillHist1D("eleCutFlow", 6, wt);
    //cout<<"sip3D: "<<std::fabs(electron.dB3D/electron.edB3D)<<endl;
    if (std::fabs(electron.dB3D/electron.edB3D) >= AnaUtil::cutValue(electronCutMap(), "SIP3D")) continue;
    AnaUtil::fillHist1D("eleCutFlow", 7, wt);
    
    // loose electrons
    looseEleList_.push_back(electron);

    // prepare a vector<pair<vhtm::Electron, vector<vhtm::PackedPFCandidate>>> for each loose electron
    // attaching an empty photon vector
    std::vector<vhtm::PackedPFCandidate> phov;
    looseElePhotonPairList_.push_back({electron, phov});
    
    // tight electrons
    if (!LL4JMETUtil::electronBDT(electron)) continue;
    AnaUtil::fillHist1D("eleCutFlow", 8, wt);
    tightEleList_.push_back(electron);
  }
  searchedEle_ = true;
}
// Jet selection
void PhysicsObjSelector::jetSelector(double wt) {
  if (!searchedPhoton_) photonSelector(wt);

  histf()->cd();
  histf()->cd("ObjectSelection");
  for (const auto& jet: *jetColl()) {
    AnaUtil::fillHist1D("jetCutFlow", 0, wt);

    AnaUtil::fillHist1D("nJets_0", jetColl()->size(), wt);

    if (jet.pt <= AnaUtil::cutValue(jetCutMap(), "pt"))                   continue;
    AnaUtil::fillHist1D("jetCutFlow", 1, wt);
    AnaUtil::fillHist1D("nJets_1", jetColl()->size(), wt);

    if (std::fabs(jet.eta) >= AnaUtil::cutValue(jetCutMap(), "eta"))      continue;
    AnaUtil::fillHist1D("jetCutFlow", 2, wt);
    AnaUtil::fillHist1D("nJets_2", jetColl()->size(), wt);

    if (!LL4JMETUtil::jetpuMVAid(jet))                                      continue; 
    AnaUtil::fillHist1D("jetCutFlow", 3, wt);
    AnaUtil::fillHist1D("nJets_3", jetColl()->size(), wt);

    if (!jetLeptonCleaning(jet, AnaUtil::cutValue(jetCutMap(), "dRlep"))) continue;
    AnaUtil::fillHist1D("jetCutFlow", 4, wt);
    AnaUtil::fillHist1D("nJets_4", jetColl()->size(), wt);

    if (!LL4JMETUtil::isLooseJet(jet))                                      continue;   
    AnaUtil::fillHist1D("jetCutFlow", 5, wt);
    AnaUtil::fillHist1D("nJets_5", jetColl()->size(), wt);
    looseJetList_.push_back(jet);
  }
  int nJets = looseJetList_.size();
  AnaUtil::fillHist1D("nJets", nJets, wt);
  if (nJets > 1)
    std::sort(std::begin(looseJetList_), std::end(looseJetList_), PtComparator<vhtm::Jet>());

  int nbj = 0;
  for (const auto& jet: looseJetList_) {
    AnaUtil::fillHist1D("bDisc", jet.pfCombinedInclusiveSecondaryVertexV2BJetTags, wt);
    if (jet.pfCombinedInclusiveSecondaryVertexV2BJetTags > AnaUtil::cutValue(jetCutMap(), "btagFactor")) nbj++;
    AnaUtil::fillHist1D("jetPt", jet.pt, wt);

    if (!LL4JMETUtil::isTightJet(jet)) continue;
    tightJetList_.push_back(jet);
    AnaUtil::fillHist1D("jetCutFlow", 6, wt);
  }
  nbJets_ = nbj;    
  AnaUtil::fillHist1D("nbJets", nbj, wt);
}
// photon from PFCandidates
void PhysicsObjSelector::photonSelector(double wt) {
  if (!searchedMu_)  muonSelector(wt);
  if (!searchedEle_) electronSelector(wt);

  // -- from twiki --
  // FSR Photon selection: for each photon, consider the closest lepton:
  // Note: Multiple photons can be attached to a lepton
  // Muons

  // Pre-selection of PF Photons
  histf()->cd();
  histf()->cd("ObjectSelection");
  if (looseEleList_.size() || looseMuList_.size()) {
    for (const auto& pfcand: *packedPFCandidateColl()) {
      if (std::abs(pfcand.pdgId) != 22) continue;
      AnaUtil::fillHist1D("photonCutFlow", 0, wt);
      
      // Preselection: pT > 2 GeV, |eta| < 2.4 
      if (pfcand.pt <= AnaUtil::cutValue(photonCutMap(), "pt")) continue; 
      AnaUtil::fillHist1D("photonCutFlow", 1, wt);

      if (std::fabs(pfcand.eta) >= AnaUtil::cutValue(photonCutMap(), "eta")) continue;
      AnaUtil::fillHist1D("photonCutFlow", 2, wt);

      // Photon PF relative isolation less than 1.8.
      // The PF isolation is computed using a cone of 0.3, a threshold of 0.2 GeV on charged 
      // hadrons with a veto cone of 0.0001, and 0.5 GeV on neutral hadrons and photons with 
      // a veto cone of 0.01, including also the contribution from PU vertices (same radius 
      // and threshold as per charged isolation) 
      double iso = LL4JMETUtil::pfiso(pfcand);
      if (iso/pfcand.pt >= AnaUtil::cutValue(photonCutMap(), "isol"))             continue;
      AnaUtil::fillHist1D("photonCutFlow", 3, wt);

      // new in 2016
      // Supercluster veto
      if (!passedSuperClusterVetobyReference(pfcand, false))                      continue;
      AnaUtil::fillHist1D("photonCutFlow", 4, wt);
      
      // Photons are associated to the closest lepton in the event among all those passing loose ID + SIP cut
      // It is constructed in such a way that we must check elindx first 
      int muindx = -1, 
	  elindx = -1;
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      if (muindx < 0 && elindx < 0)                                               continue;
      AnaUtil::fillHist1D("photonCutFlow", 5, wt);
      
      if (dRmin >= AnaUtil::cutValue(photonCutMap(), "dRmin"))                    continue;
      AnaUtil::fillHist1D("photonCutFlow", 6, wt);

      // Discard photons that do not satisfy the cuts dR(pho,l)/ET_pho^2 < 0.012, and dR(pho,l) < 0.5. 
      TLorentzVector pfcandP4(LL4JMETUtil::getP4(pfcand)); 
      double dRovEt2 = dRmin/pfcandP4.Et2();
      if (dRovEt2 >= AnaUtil::cutValue(photonCutMap(), "dRovEt2"))                continue;
      AnaUtil::fillHist1D("photonCutFlow", 7, wt);

      AnaUtil::fillHist1D("phoPt", pfcand.pt, wt);
      // If more than one photon is associated with the same lepton, only the lowest dR(pho,l)/ET_pho^2 is stored
      if (elindx > -1) {  // check electron index first
        if (looseElePhotonPairList_.at(elindx).second.empty())
	  looseElePhotonPairList_.at(elindx).second.push_back(pfcand);
        else {
          TLorentzVector prephoP4(LL4JMETUtil::getP4(looseElePhotonPairList_.at(elindx).second.at(0)));
          TLorentzVector eleP4(LL4JMETUtil::getP4(looseElePhotonPairList_.at(elindx).first));
          if (eleP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
            looseElePhotonPairList_.at(elindx).second.clear(); // reset all
            looseElePhotonPairList_.at(elindx).second.push_back(pfcand);
          }
        }
      }
      else if (muindx > -1) {
        if (looseMuPhotonPairList_.at(muindx).second.empty())
	  looseMuPhotonPairList_.at(muindx).second.push_back(pfcand);
        else {
          TLorentzVector prephoP4(LL4JMETUtil::getP4(looseMuPhotonPairList_.at(muindx).second.at(0)));
          TLorentzVector muP4(LL4JMETUtil::getP4(looseMuPhotonPairList_.at(muindx).first));
          if (muP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
            looseMuPhotonPairList_.at(muindx).second.clear(); // reset all
            looseMuPhotonPairList_.at(muindx).second.push_back(pfcand);
          }
        }
      }
    }
  }
  // collect FSR photon candidates in fsrPhotonList_
  // also fill the vectors for tight leptons
  for (const auto& elem: looseElePhotonPairList_) {
    if (!elem.second.empty())
      fsrPhotonList_.push_back(elem.second.at(0));

    // tight electron
    if (LL4JMETUtil::electronBDT(elem.first))
      tightElePhotonPairList_.push_back(elem);
  }
  for (const auto& elem: looseMuPhotonPairList_) {
    if (!elem.second.empty()) 
      fsrPhotonList_.push_back(elem.second.at(0));

    // tight muon 
    const auto& muon = elem.first;
    bool highPtId = muon.passTrackerhighPtid && muon.pt > 200.;
    if (muon.isPFMuon || highPtId)
      tightMuPhotonPairList_.push_back(elem);
  }
  searchedPhoton_ = true;
}
void PhysicsObjSelector::isoLeptonSelector(double muIso, double eleIso) {
  // Tight Muons
  for (const auto& elem: tightMuPhotonPairList_) {
    if (LL4JMETUtil::computeMuonReliso(elem.first, fsrPhotonList_, 0.01, 0.3) >= muIso) continue;
    tightIsoMuList_.push_back(elem.first);
    tightIsoMuPhotonPairList_.push_back(elem);
  }

  // Tight electrons
  for (const auto& elem: tightElePhotonPairList_) {
    if (LL4JMETUtil::computeElectronReliso(elem.first, fsrPhotonList_, getEventGridRho(), 0.08, 0.3) >= eleIso) continue;
    tightIsoEleList_.push_back(elem.first);
    tightIsoElePhotonPairList_.push_back(elem);
  }
}
void PhysicsObjSelector::leptonCrossCleaning() {
  vector<pair<vhtm::Electron, vector<vhtm::PackedPFCandidate>>> list;
  for (const auto& elem: tightElePhotonPairList_)
    if (crossCleaned(elem.first)) list.push_back(elem);

  tightElePhotonPairList_.clear();
  tightElePhotonPairList_ = list;
}
// from twiki
// Lepton cross cleaning: Remove electrons which are 
// within dR(eta,phi)<0.05 of a muon passing tight ID && SIP<4 
bool PhysicsObjSelector::crossCleaned(const vhtm::Electron& electron) const {
  bool flag = true;
  for (const auto& mu: tightMuList_) {
    if (LL4JMETUtil::getP4(electron).DeltaR(LL4JMETUtil::getP4(mu)) < 0.05) {
      flag = false;
      break;
    }
  }
  return flag;
}
bool PhysicsObjSelector::jetLeptonCleaning(const vhtm::Jet& jet, double dR) const {
  TLorentzVector jetP4(LL4JMETUtil::getP4(jet));
  // Tight muons
  for (const auto& mu: tightIsoMuList_){
    AnaUtil::fillHist1D("JetMuDR", jetP4.DeltaR(LL4JMETUtil::getP4(mu)), 1.0);
    if (jetP4.DeltaR(LL4JMETUtil::getP4(mu)) <= dR) return false;
  }
  // Tight electrons
  for (const auto& ele: tightIsoEleList_){
    AnaUtil::fillHist1D("JetEleDR", jetP4.DeltaR(LL4JMETUtil::getP4(ele)), 1.0);
    if (jetP4.DeltaR(LL4JMETUtil::getP4(ele)) <= dR) return false;
  }
  // clean against FSR also
  for (const auto& v: fsrPhotonList_) {
    // photon isolation is inbuilt
    TLorentzVector fsrP4(LL4JMETUtil::getP4(v));
    if (fsrP4.Pt() > 10 && jetP4.DeltaR(fsrP4) <= dR) return false;
  }
  return true;
}
bool PhysicsObjSelector::passedSuperClusterVeto(const vhtm::PackedPFCandidate& pfcand, bool verbose) const {
  // Supercluster veto: remove all PF photons that match with any electron passing loose ID and SIP cuts; 
  // matching is according to (|deta| < 2, |dphi| < 0.05) OR (dR < 0.15), with respect to the electron's supercluster. 
  bool passedVeto {true};
  if (verbose && looseEleList_.size()) {
    cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << endl;
    cout << setprecision(3);
  }
  for (const auto& ele: looseEleList_) {
    double deta = ele.scEta - pfcand.eta;
    double dphi = TVector2::Phi_mpi_pi(ele.scPhi - pfcand.phi);
    double dR = std::sqrt(deta * deta + dphi * dphi);
    if (verbose) {
      cout << setw(8) << pfcand.pt
	   << setw(8) << pfcand.eta
	   << setw(8) << pfcand.phi
	   << setw(8) << ele.pt
	   << setw(8) << ele.scEta
	   << setw(8) << ele.scPhi
	   << setw(8) << deta
	   << setw(8) << dphi
	   << setw(8) << dR
	   << endl;
    }
    if ( (std::fabs(deta) < 0.05 && std::fabs(dphi) < 2.0) || dR < 0.15 ) {
      passedVeto = false;
      break;
    }
  }
  return passedVeto;
}
bool PhysicsObjSelector::passedSuperClusterVetobyReference(const vhtm::PackedPFCandidate& pfcand, bool verbose) const {
  // Supercluster veto by PF reference: veto all the PF candidates used in the PF cluster, as returned by the method 
  bool passedVeto {true};
  TLorentzVector pfcandP4(LL4JMETUtil::getP4(pfcand));
  if (verbose && looseEleList_.size())
    cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << endl;
  for (const auto& ele: looseEleList_) {
    for (const auto& p4ref: ele.associatedPackedPFCandidatesP4) {
      if (AnaUtil::sameObject(p4ref, pfcandP4)) {
        passedVeto = false;
        break;
      }
    }
    if (!passedVeto) break;
    if (verbose) {
      cout << setprecision(3);
      cout << setw(8) << pfcand.pt
	   << setw(8) << pfcand.eta
	   << setw(8) << pfcand.phi
	   << setw(8) << ele.pt
	   << setw(8) << ele.scEta
	   << setw(8) << ele.scPhi
	   << endl;
    }
  }
  return passedVeto;
}
double PhysicsObjSelector::findClosestLepton(const vhtm::PackedPFCandidate& pfPho, int& muindx, int& elindx) const {
  TLorentzVector phoP4(LL4JMETUtil::getP4(pfPho));
  double dRmin {999};
  muindx = -1;
  // Consider loose muons first
  for (size_t i = 0; i < looseMuPhotonPairList_.size(); ++i) {
    const vhtm::Muon& muon = looseMuPhotonPairList_[i].first;
    TLorentzVector muP4(LL4JMETUtil::getP4(muon));
    double dR = muP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      muindx = i;
    }
  }
  // Then consider loose electron
  elindx = -1;
  for (size_t i = 0; i < looseElePhotonPairList_.size(); ++i) {
    const vhtm::Electron& ele = looseElePhotonPairList_[i].first;
    TLorentzVector eleP4(LL4JMETUtil::getP4(ele));
    double dR = eleP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      elindx = i;
    }
  }
  return dRmin;
}
/*void PhysicsObjSelector::ZZMass(ZCandidate& Za, ZCandidate& Zb, std::vector<std::pair<ZCandidate, ZCandidate> >& ZZVec) {
  // -- from twiki --
  // Require both Z candidate masses (computed including the FSR photons, if present) to be 12 < m(ll(g)) < 120 GeV
  // leptons are already isolated in the new FSR scheme
  // Z pair pass the mass cut, we have a ZZ Candidate
  bool ZZmasscond = (Za.mass > 12 && Za.mass < 120) && (Zb.mass > 12 && Zb.mass < 120);
  
  // -- from twiki --
  // define the Z1 as the one with mass closest to the nominal mZ; require mZ1 > 40 GeV. The other Z is the Z2.
  if (ZZmasscond) {
    if (Za.massDiff < Zb.massDiff) {
      if (Za.mass > 40.)
  	ZZVec.push_back({Za, Zb});
    } 
    else {
      if (Zb.mass > 40.)
	ZZVec.push_back({Zb, Za});
    }
  }
  }*/
/*void PhysicsObjSelector::addLeptonIsolation(const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& elePhotonPairList, 
					    const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& muPhotonPairList, 
					    std::vector<ZCandidate>& ZCandList) 
{
  for (auto& v: ZCandList) {
    if (v.flavour == LL4JMETUtil::ZType::mumu) {
      const vhtm::Muon& mu1 = muPhotonPairList.at(v.l1Index).first;
      v.l1Isolation = LL4JMETUtil::computeMuonReliso(mu1, fsrPhotonList_, 0.01, 0.3);

      const vhtm::Muon& mu2 = muPhotonPairList.at(v.l2Index).first;
      v.l2Isolation = LL4JMETUtil::computeMuonReliso(mu2, fsrPhotonList_, 0.01, 0.3);
    }
    else if (v.flavour == LL4JMETUtil::ZType::ee) {
      const vhtm::Electron& ele1 = elePhotonPairList.at(v.l1Index).first;
      v.l1Isolation = LL4JMETUtil::computeElectronReliso(ele1, fsrPhotonList_, getEventGridRho(), 0.08, 0.3);

      const vhtm::Electron& ele2 = elePhotonPairList.at(v.l2Index).first;
      v.l2Isolation = LL4JMETUtil::computeElectronReliso(ele2, fsrPhotonList_, getEventGridRho(), 0.08, 0.3);
    }
  }
  }
void 
PhysicsObjSelector::addLeptonIsolationForLL(const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& elePhotonPairList,
					    const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& muPhotonPairList, 
					    std::vector<ZCandidate>& ExoCandList) 
{
  for (auto& v: ExoCandList) {
    if (v.flavour == LL4JMETUtil::llType::mue) {
      const vhtm::Muon& mu = muPhotonPairList.at(v.l1Index).first;
      v.l1Isolation = LL4JMETUtil::computeMuonReliso(mu, fsrPhotonList_, 0.01, 0.3);

      const vhtm::Electron& ele = elePhotonPairList.at(v.l2Index).first;
      v.l2Isolation = LL4JMETUtil::computeElectronReliso(ele, fsrPhotonList_, getEventGridRho(), 0.08, 0.3);
    }
    else if (v.flavour == LL4JMETUtil::llType::emu) {
      const vhtm::Electron& ele = elePhotonPairList.at(v.l1Index).first;
      v.l1Isolation = LL4JMETUtil::computeElectronReliso(ele, fsrPhotonList_, getEventGridRho(), 0.08, 0.3);

      const vhtm::Muon& mu = muPhotonPairList.at(v.l2Index).first;
      v.l2Isolation = LL4JMETUtil::computeMuonReliso(mu, fsrPhotonList_, 0.01, 0.3);
    }
  }
}
void 
PhysicsObjSelector::addLeptonIsolationForLTau(const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& elePhotonPairList, 
					      const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& muPhotonPairList, 
					      std::vector<ZCandidate>& ZCandList) 
{
  for (auto& v: ZCandList) {
    if (v.flavour == LL4JMETUtil::lTauType::mutau) {  // light lepton is the first entry
      const vhtm::Muon& mu = muPhotonPairList.at(v.l1Index).first;
      v.l1Isolation = LL4JMETUtil::computeMuonReliso(mu, fsrPhotonList_, 0.01, 0.3);
      v.l2Isolation = -1.0;
    }
    else if (v.flavour == LL4JMETUtil::lTauType::etau) {
      const vhtm::Electron& ele = elePhotonPairList.at(v.l1Index).first;
      v.l1Isolation = LL4JMETUtil::computeElectronReliso(ele, fsrPhotonList_, getEventGridRho(), 0.08, 0.3);
      v.l2Isolation = -1.0;
    }
  }
  }*/
void PhysicsObjSelector::addLeptonIsolation(std::vector<LeptonCandidate>& lepCandList,
                                            const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& elePhotonPairList,
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& muPhotonPairList)
{
  for (auto& v: lepCandList) {
    if (v.flavour == 1) {
      const vhtm::Muon& mu = muPhotonPairList.at(v.lIndex).first;
      v.lIsolation = LL4JMETUtil::computeMuonReliso(mu, fsrPhotonList_, 0.01, 0.3);
      //      if (v.lIsolation < 0.35)                                                                                                                                  
    }
    else if (v.flavour == 2) {
      const vhtm::Electron& ele = elePhotonPairList.at(v.lIndex).first;
      v.lIsolation = LL4JMETUtil::computeElectronReliso(ele, fsrPhotonList_, getEventGridRho(), 0.08, 0.3);
    }
  }
}


int PhysicsObjSelector::getTrueLooseLeptons(double dRCut) const {
  int ntlep {0};
  for (const auto& mu: looseMuList_) {
    TLorentzVector muP4(LL4JMETUtil::getP4(mu));
    bool res = true;
    for (const auto& jet: looseJetList_) {
      TLorentzVector jetP4(LL4JMETUtil::getP4(jet));
      if (jetP4.DeltaR(muP4) <= dRCut) {
        res = false;
        break;
      }
    }
    if (res) ntlep++;
  }
  for (const auto& ele: looseEleList_) {
    TLorentzVector eleP4(LL4JMETUtil::getP4(ele));
    bool res = true;
    for (const auto& jet: looseJetList_) {
      TLorentzVector jetP4(LL4JMETUtil::getP4(jet));
      if (jetP4.DeltaR(eleP4) <= dRCut) {
        res = false;
        break;
      }
    }
    if (res) ntlep++;
  }
  return ntlep;
}
std::tuple<int,int,int,int> PhysicsObjSelector::findIsoLeptons(double dRCut) const {
  // Muon
  int nLooseIsoMu {0},
      nTightIsoMu {0};
  for (const auto& l: looseMuList_) {
    TLorentzVector lP4(LL4JMETUtil::getP4(l));
    bool isInsideJet = false;
    for (const auto& j: looseJetList_) {
      TLorentzVector jP4(LL4JMETUtil::getP4(j));
      double dR = lP4.DeltaR(jP4);
      if (dR <= dRCut) {
	isInsideJet = true;
	break;
      }
    }
    if (isInsideJet) continue;
    double iso = LL4JMETUtil::computeMuonReliso(l, fsrPhotonList_, 0.01, 0.3)/lP4.Pt();
    if (iso >= 0.35) continue;
    nLooseIsoMu++;
    
    // tight muon 
    bool isTight = l.isPFMuon || (l.passTrackerhighPtid && l.pt > 200.);
    if (!isTight) continue;
    nTightIsoMu++;
  }

  // Electrons
  int nLooseIsoEle {0},
      nTightIsoEle {0};
  for (const auto& l: looseEleList_) {
    TLorentzVector lP4(LL4JMETUtil::getP4(l));
    bool isInsideJet = false;
    for (const auto& j: looseJetList_) {
      TLorentzVector jP4(LL4JMETUtil::getP4(j));
      double dR = lP4.DeltaR(jP4);
      if (dR <= dRCut) {
	isInsideJet = true;
	break;
      }
    }
    if (isInsideJet) continue;
    double iso = LL4JMETUtil::computeElectronReliso(l, fsrPhotonList_, getEventGridRho(), 0.08, 0.3)/lP4.Pt();
    if (iso >= 0.35) continue;
    nLooseIsoEle++;
    
    if (!LL4JMETUtil::electronBDT(l)) continue;
    nTightIsoEle++;
  }
  return std::make_tuple(nLooseIsoMu, nTightIsoMu, nLooseIsoEle, nTightIsoEle);
}
void PhysicsObjSelector::dumpEvent(double vz, bool dumpGen, bool showEvent, ostream& os) const {
  os << std::setprecision(3);

  // Event
  if (showEvent) showEventNumber(os);

  // Muons
  if (muonColl()->size()) {
    os << " -- # Muons: " << muonColl()->size() << endl;
    os << "  indx      pT     eta     phi  charge      dxy       dz  global tracker      PF  nMatch  Type       SIP ghostCleaned Loose Tight  reliso"
       << endl;
    int indx {0};
    for (const auto& muon: *muonColl()) {
      bool muType = muon.isGlobalMuon || (muon.isTrackerMuon && muon.matches > 0);
      bool isLoose = (muon.pt > AnaUtil::cutValue(muonCutMap(), "pt")                     &&
		      std::fabs(muon.eta) < AnaUtil::cutValue(muonCutMap(), "eta")        &&
		      std::fabs(muon.dxyPV) < AnaUtil::cutValue(muonCutMap(), "dxyPV")    &&
		      std::fabs(muon.dzPV) < AnaUtil::cutValue(muonCutMap(), "dzPV")      && 
		      muType                                                              &&
		      muon.muonBestTrackType != 2                                         &&
		      muon.isghostCleaned                                                 &&
		      std::fabs(muon.dB3D/muon.edB3D) < AnaUtil::cutValue(muonCutMap(), "SIP3D"));
      bool isTight = isLoose && (muon.isPFMuon || (muon.passTrackerhighPtid && muon.pt > 200.));
      os << setw(6)  << indx++ 
	 << setw(8)  << muon.pt 
	 << setw(8)  << muon.eta
	 << setw(8)  << muon.phi
         << setw(8)  << muon.charge
	 << setw(9)  << muon.dxyPV
	 << setw(9)  << muon.dzPV
	 << setw(8)  << (muon.isGlobalMuon ? "T" : "F")
	 << setw(8)  << (muon.isTrackerMuon ? "T" : "F")
	 << setw(8)  << (muon.isPFMuon ? "T" : "F")
	 << setw(8)  << muon.matches
	 << setw(6)  << muon.muonBestTrackType
	 << setw(10) << muon.dB3D/muon.edB3D
	 << setw(13) << (muon.isghostCleaned ? "T" : "F")
	 << setw(6) << (isLoose ? "T" : "F")
	 << setw(6) << (isTight ? "T" : "F")
         << setw(8) << LL4JMETUtil::computeMuonReliso(muon, fsrPhotonList_, 0.01, 0.3)
	 << endl;
    }
  }
  // Electrons
  if (electronColl()->size()) {
    os << " -- # Electrons: " << electronColl()->size() << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz  misHit       SIP   BDT crossCleaned Loose Tight  reliso"
       << endl; 
    int indx {0};
    for (const auto& electron: *electronColl()) {
      bool isLoose = (electron.pt > AnaUtil::cutValue(electronCutMap(), "pt")                   &&
		      std::fabs(electron.eta) < AnaUtil::cutValue(electronCutMap(), "eta")      &&
		      std::fabs(electron.dxyPV) < AnaUtil::cutValue(electronCutMap(), "dxyPV")  &&
		      std::fabs(electron.dzPV) < AnaUtil::cutValue(electronCutMap(), "dzPV")    &&
		      crossCleaned(electron)                                                    &&
		      std::fabs(electron.dB3D/electron.edB3D) < AnaUtil::cutValue(electronCutMap(), "SIP3D"));
      bool isTight = isLoose && LL4JMETUtil::electronBDT(electron);

      os << setw(6)  << indx++
	 << setw(8)  << electron.pt 
	 << setw(8)  << electron.eta
	 << setw(8)  << electron.phi
         << setw(8)  << electron.charge
	 << setw(8)  << electron.dxyPV
	 << setw(8)  << electron.dzPV
	 << setw(8)  << electron.missingHits
	 << setw(10) << electron.dB3D/electron.edB3D
	 << setw(6)  << (LL4JMETUtil::electronBDT(electron) ? "T" : "F")
         << setw(13) << (crossCleaned(electron) ? "T" : "F")
	 << setw(6)  << (isLoose ? "T" : "F")
	 << setw(6)  << (isTight ? "T" : "F")
         << setw(8)  << LL4JMETUtil::computeElectronReliso(electron, fsrPhotonList_, getEventGridRho(), 0.08, 0.3)
	 << endl;
    }
  }
  // Photons
  int npho {0};
  for (const auto& pfcand: *packedPFCandidateColl()) {
    // Preselection: pT > 2 GeV, |eta| < 2.4 
    if (pfcand.pdgId != 22 || pfcand.pt < 2. || std::fabs(pfcand.eta) >= 2.4) continue;
    ++npho;
  }
  if (npho) {
    os << " -- Photons: " << endl;
    os << "  indx      pT     eta     phi  reliso  SCVeto closestLepton   dRmin        Et2 dRovEt2(%) selected"
       << endl; 
    int indx {0};
    for (const auto& pfcand: *packedPFCandidateColl()) {
      // Preselection: pT > 2 GeV, |eta| < 2.4 
      if (pfcand.pdgId != 22 || pfcand.pt < 2. || std::fabs(pfcand.eta) >= 2.4) continue;
      int muindx = -1, 
	  elindx = -1;
      ostringstream ltype;
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      if (elindx > -1) ltype << "electron(" << elindx << ")";
      else if (muindx > -1) ltype << "muon(" << muindx << ")";
      else ltype << "none";

      double iso = LL4JMETUtil::pfiso(pfcand);
      bool scVetoPassed = passedSuperClusterVetobyReference(pfcand);
      TLorentzVector pfcandP4(LL4JMETUtil::getP4(pfcand)); 
      double dRovEt2 = dRmin/pfcandP4.Et2();
      bool selected = (iso/pfcand.pt < 1.8          &&  
		       scVetoPassed 		    &&  
		       (muindx > -1 || elindx > -1) &&  
		       dRmin < 0.5 		    &&  
		       dRovEt2 < 0.012);
      os << setw(6)  << indx++
	 << setw(8)  << pfcand.pt 
	 << setw(8)  << pfcand.eta
	 << setw(8)  << pfcand.phi
	 << setw(8)  << iso/pfcand.pt
	 << setw(8)  << (scVetoPassed ? "T" : "F")
	 << setw(14) << ltype.str()
	 << setw(8)  << dRmin
         << setprecision(2)
	 << setw(11) << pfcandP4.Et2()
         << setprecision(3)
	 << setw(11) << dRovEt2*100
	 << setw(9) << (selected ? "T" : "F")
	 << endl;
    }
  }
  // Taus
  if (tauColl()->size()) {
    os << " -- # Taus: " << tauColl()->size() << endl;
    os << "  indx       pT      eta      phi charge decayMode isolation muonVeto eleVeto  zvertex     dzPV selected"
       << endl;
    int indx {0};
    for (const auto& tau: *tauColl()) {
      bool selected = (std::fabs(tau.dzPV) < AnaUtil::cutValue(tauCutMap(), "dz")               &&
		       tau.pt > AnaUtil::cutValue(tauCutMap(), "pt")                            &&
		       std::fabs(tau.eta) < AnaUtil::cutValue(tauCutMap(), "eta")               &&
		       tau.decayModeFinding > 0.99                                              &&
		       tau.againstMuonTight3 > AnaUtil::cutValue(tauCutMap(), "muVeto")         &&
		       tau.againstElectronTightMVA > AnaUtil::cutValue(tauCutMap(), "eleVeto")  &&                                            
		       tau.byMediumIsolationMVArun2v1DBoldDMwLT > AnaUtil::cutValue(tauCutMap(), "isol"));                     

      os << setw(6) << indx++
	 << setw(9) << tau.pt
	 << setw(9) << tau.eta
	 << setw(9) << tau.phi
	 << setw(7) << tau.charge
         << setprecision(1)
         << setw(10) << ((tau.decayModeFinding > 0.999) ? "T" : "F")
         << setw(10) << ((tau.byLooseIsolationMVArun2v1DBoldDMwLT > 0.5) ? "T" : "F")
         << setw(9) << ((tau.againstMuonTight3 > 0.5) ? "T" : "F")
         << setw(8) << ((tau.againstElectronLooseMVA > 0.5) ? "T" : "F")
         << setprecision(3)
	 << setw(9) << tau.zvertex - vz
         << setw(9) << tau.dzPV
	 << setw(9) << (selected ? "T" : "F")
	 << endl;
    }
  }
  // Jets
  if (jetColl()->size()) {
    os << " -- # Jets: " << jetColl()->size() << endl;
    os << "  indx       pT      eta      phi NConst   CHM      CHF     CEMF      NHF     NEMF      MUF puID   bDisc lepCleaned  looseId  tightId ljet bjet"
       << endl;
    int indx {0};
    for (const auto& jet: *jetColl()) {
      bool is_ljet = (jet.pt > AnaUtil::cutValue(jetCutMap(), "pt")                   &&
		      std::fabs(jet.eta) < AnaUtil::cutValue(jetCutMap(), "eta")      && 
		      LL4JMETUtil::jetpuMVAid(jet)                                      && 
		      jetLeptonCleaning(jet, AnaUtil::cutValue(jetCutMap(), "dRlep")) &&
		      LL4JMETUtil::isLooseJet(jet));
      bool is_bjet = is_ljet && jet.pfCombinedInclusiveSecondaryVertexV2BJetTags > AnaUtil::cutValue(jetCutMap(), "btagFactor");
      os << setw(6)  << indx++
	 << setw(9)  << jet.pt
	 << setw(9)  << jet.eta
	 << setw(9)  << jet.phi
         << setw(7)  << (jet.chargedMultiplicity + jet.neutralMultiplicity)
         << setw(6)  << jet.chargedMultiplicity
         << setw(9)  << jet.chargedHadronEnergyFraction
         << setw(9)  << jet.chargedEmEnergyFraction
         << setw(9)  << jet.neutralHadronEnergyFraction
         << setw(9)  << jet.neutralEmEnergyFraction
         << setw(9)  << jet.muonEnergyFraction
	 << setw(5)  << (LL4JMETUtil::jetpuMVAid(jet) ? "T" : "F")
	 << setw(8)  << jet.pfCombinedInclusiveSecondaryVertexV2BJetTags
	 << setw(11) << (jetLeptonCleaning(jet, AnaUtil::cutValue(jetCutMap(), "dRlep")) ? "T" : "F")
         << setw(9)  << (LL4JMETUtil::isLooseJet(jet) ? "T" : "F")
         << setw(9)  << (LL4JMETUtil::isTightJet(jet) ? "T" : "F")
	 << setw(5)  << ((is_ljet) ? "T" : "F")
	 << setw(5)  << ((is_bjet) ? "T" : "F")
	 << endl;
    }
  }
  // Selected vertices
  if (vertexColl()->size()) {
    os << " -- # Vertices: " << vertexColl()->size() << endl;
    os << "  indx      ndf     rho    chi2     dxy       z  isFake"
       << endl; 
    int indx {0};
    for (const auto& vtx: *vertexColl()) {
      double dxy = std::sqrt(pow(vtx.x, 2) + pow(vtx.y, 2));
      os << setw(6) << indx++
	 << setw(9) << vtx.ndf
	 << setw(8) << vtx.rho
	 << setw(8) << vtx.chi2
	 << setw(8) << dxy
	 << setw(8) << vtx.z
	 << setw(8) << (vtx.isfake ? "T" : "F")
	 << endl;
    }
  }
  // PF Candidates
  if (packedPFCandidateColl()->size()) {
    os << " -- # PFCandidates: " << packedPFCandidateColl()->size() << endl;
    os << "  indx   pdgId      pT     eta     phi isolation"
       << endl; 
    int indx {0};
    for (const auto& pfcand: *packedPFCandidateColl()) {
      double iso = LL4JMETUtil::pfiso(pfcand);
      os << setw(6)  << indx++
	 << setw(8)  << pfcand.pdgId 
	 << setw(8)  << pfcand.pt 
	 << setw(8)  << pfcand.eta
	 << setw(8)  << pfcand.phi
	 << setw(10) << iso/pfcand.pt
	 << endl;
    }
  }
  // MET
  os << " -- # MET: " << endl
     << "      type      met      phi    sumet"
     << endl; 
  os << setw(10) << "pfMet"
     << setw(9) << metColl()->at(0).met
     << setw(9) << metColl()->at(0).metphi
     << setw(9) << metColl()->at(0).sumet
     << endl;
  os << setw(10) << "corrMet"
     << setw(9) << corrmetColl()->at(0).met
     << setw(9) << corrmetColl()->at(0).metphi
     << setw(9) << corrmetColl()->at(0).sumet
     << endl;
  os << setw(10) << "puppiMet"
     << setw(9) << puppimetColl()->at(0).met
     << setw(9) << puppimetColl()->at(0).metphi
     << setw(9) << puppimetColl()->at(0).sumet
     << endl << endl;

  if (isMC() && dumpGen) dumpGenInfo(os); 
}
