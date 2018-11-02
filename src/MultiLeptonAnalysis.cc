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
#include <memory>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"

#include "HZZ4lUtil.h"
#include "MultiLeptonAnalysis.h"
#include "GenAnalysis.h"

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

using namespace vhtm;

// -----------
// Constructor
// -----------
MultiLeptonAnalysis::MultiLeptonAnalysis()
  : PhysicsObjSelector()
{
}
// ----------
// Destructor
// ----------
MultiLeptonAnalysis::~MultiLeptonAnalysis() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MultiLeptonAnalysis::beginJob() 
{ 
  if (!PhysicsObjSelector::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("MultiLeptonAnalysis");
  
  bookHistograms();
  if (genAna_ != nullptr) genAna_->bookHistograms(histf());

#ifdef SKIP_DUPLICATE_ALL
  eventIdStore_.clear();
#endif

  return true;
}
// ---------------
// Book histograms
// ---------------
void MultiLeptonAnalysis::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("MultiLeptonAnalysis");

  // book histograms to be filled at different stages
  new TH1D("nvtx_stage0", "Number of Good vertices (stage 0)", 60, 0, 60);
  new TH1D("nvtx_stage1", "Number of Good vertices (stage 1)", 60, 0, 60);
  new TH1D("nvtx_stage2", "Number of Good vertices (stage 2)", 60, 0, 60);
  new TH1D("nvtx_stage3", "Number of Good vertices (stage 3)", 60, 0, 60);
  new TH1D("nvtx_stage4", "Number of Good vertices (stage 4)", 60, 0, 60);
  new TH1D("nvtx_stage5", "Number of Good vertices (stage 5)", 60, 0, 60);
  new TH1D("nvtx_stage6", "Number of Good vertices (stage 6)", 60, 0, 60);
  new TH1D("nvtx_stage7", "Number of Good vertices (stage 7)", 60, 0, 60);
  new TH1D("nvtx_stage8", "Number of Good vertices (stage 8)", 60, 0, 60);
  if (isMC()) {
    new TH1D("puweight", "PU reweight factor", 100, 0, 2);
    new TH1D("evtweight", "Event weight factor (MC events)", 100, -5, 5);
  }

  //------- Object PLots -----------------------------------------------
  new TH1D("nGoodMuon", "Number of Good muons", 6, -0.5, 5.5);
  new TH1D("nGoodEle",  "Number of Good electrons", 6, -0.5, 5.5);
  new TH1D("nGoodTau",  "Number of Good Tau", 5, -0.5, 4.5);
  new TH1D("nJets",     "Number of Tight Jets", 10, -0.5, 9.5);
  new TH1D("nbJets",    "Number of b-Jets", 10, -0.5, 9.5);

  new TH1D("evtCutFlow", "Event CutFlow", 15, -0.5, 14.5);
  if (isMC()) 
    new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 15, -0.5, 14.5);

  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);

  new TH1D("nLooseLeptons",    "Total tight isolated leptons", 10, -0.5, 9.5);
  new TH1D("nTightIsoLeptons", "Total tight isolated leptons", 8, -0.5, 7.5);
  new TH1D("hLepPt_stage0",    "Highest pT lepton pT (stage0)", 300, 0, 300);
  new TH1D("hLepPt",           "Highest pT lepton pT", 300, 0, 300);
  new TH1D("met",              "Missing Transverse Energy", 150, 0, 150);
  new TH1D("corrmet",          "Missing Transverse Energy (Type-I Corrected)", 150, 0, 150);
  new TH1D("puppimet",         "Missing Transverse Energy (PUPPI)", 150, 0, 150);

  // Z and h1
  new TH1F("nZcand", "Number of selected Z candidates", 10, 0, 10);
  new TH1F("massZcand", "Mass of selected Z candidates", 100, 0., 200.);

  new TH1F("Zmass",   "ll invariant mass of the Z candidate", 200, 0., 200.);
  new TH1F("Z2mass",  "ll invariant mass of the second Z candidate", 200, 0., 200.);
  new TH1F("ZPt",     "pT of the ll system of the Z candidate", 200, 0., 400.);
  new TH1F("ZmassC0", "ll invariant mass of the Z candidate (mmem)", 200, 0., 200.);
  new TH1F("ZmassC1", "ll invariant mass of the Z candidate (eeem)", 200, 0., 200.);

  new TH1F("h1mass",   "mass of the selected Higgs candidate", 300, 0, 300);
  new TH1F("h1Pt",     "pT of the e-mu system", 200, 0., 400.);
  new TH1F("h1massC0", "Mass of the selected Higgs candidate (mmem)", 300, 0, 300);
  new TH1F("h1massC1", "Mass of the selected Higgs candidate (eeem)", 300, 0, 300);

  // angle
  new TH1F("dRh1Z",   "dR(h1, Z)", 100, 0, 2 * M_PI);
  new TH1F("dRh1ZC0", "dR(h1, Z) (mmem)", 100, 0, 2 * M_PI);
  new TH1F("dRh1ZC1", "dR(h1, Z) (eeem)", 100, 0, 2 * M_PI);

  new TH1F("dEtall",  "dEta(lep, lep) of Z", 100, -5.0, 5.0);
  new TH1F("dEtaemu", "dEta(ele, mu) of h1 candidate", 100, -5.0, 5.0);

  new TH1F("dPhill",  "dPhi(lep, lep) of Z", 100, 0, 2 * M_PI);
  new TH1F("dPhiemu", "dPhi(ele, mu) of h1 candidate", 100, 0, 2 * M_PI);

  new TH1F("dRll",  "dR(lep, lep) of Z", 100, 0, 2 * M_PI);
  new TH1F("dRemu", "dR(ele, mu) of h1 candidate", 100, 0, 2 * M_PI);

  // Lepton properties (pT, eta etc.)
  new TH1F("Zl1Pt", "pT of the first lepton from Z", 300, 0, 300);
  new TH1F("Zl2Pt", "pT of the second lepton from Z", 300, 0, 300);
  new TH1F("h1ElePt", "pT of the electron from the h1 candidate", 300, 0, 300);
  new TH1F("h1MuPt",  "pT of the muon from the h1 candidate", 300, 0, 300);

  new TH1F("Zl1Eta", "Eta of the first lepton from Z", 200, -4., 4);
  new TH1F("Zl2Eta", "Eta of the second lepton from Z", 200, -4, 4);
  new TH1F("h1EleEta", "Eta of the electron from the h1 candidate", 200, -4, 4);
  new TH1F("h1MuEta",  "Eta of the muon from the h1 candidate", 100, -4, 4);

  new TH1F("Zl1Phi", "Phi of the first lepton from Z", 200, -M_PI, M_PI);
  new TH1F("Zl2Phi", "Phi of the second lepton from Z", 200, -M_PI, M_PI);
  new TH1F("h1ElePhi", "Phi of the electron from the h1 candidate", 200, -M_PI, M_PI);
  new TH1F("h1MuPhi",  "Phi of the muon from the h1 candidate", 100, -M_PI, M_PI);

  new TH1F("h1PtRatio", "pT(ele+muon) over pT(electron)+pT(muon)", 100, 0, 1.2);
  new TH1F("h1Lep1Angle", "Angle of lepton 1 from Z with the e-mu decay plane", 100, 0, M_PI);
  new TH1F("h1Lep2Angle", "Angle of lepton 2 from Z with the e-mu decay plane", 100, 0, M_PI);

  new TH1F("nEventsPerType", "Selected events per category", 3, -1.5, 1.5);
  if (isMC()) 
    new TH1F("nEventsPerTypeWt", "Selected events per category (Weighted)", 3, -1.5, 1.5);

  new TH1F("sLTSum", "Scalar sum of Transverse Energy of the leptons", 200, 0, 800);
  new TH1F("vLTSum", "Vector sum of Transverse Energy of the leptons", 100, 0, 500);

  new TH2F("h1VsZpT", "h1 vs Z pT scatter plot", 100, 0, 400, 100, 0, 400);
  new TH2F("lep2VsLep1pT", "lepton 1 vs lepton 1 pT scatter plot", 100, 0, 300, 100, 0, 300);
  new TH2F("muVsLep1pT", "muon vs lepton 1 pT scatter plot", 100, 0, 300, 100, 0, 300);
  new TH2F("muVsLep2pT", "muon vs lepton 2 pT scatter plot", 100, 0, 300, 100, 0, 300);
  new TH2F("muVsElepT", "muon vs electron pT scatter plot", 100, 0, 300, 100, 0, 300);

  histf()->cd();
  histf()->ls();
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void MultiLeptonAnalysis::clearLists() {
  PhysicsObjSelector::clear();
  vtxList_.clear();
}
// -------------------
// The main event loop
// -------------------
void MultiLeptonAnalysis::eventLoop()
{
  int nPrint = std::max(10000, nEvents()/1000);
  
  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  
  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  int fevt = (firstEvent() > -1) ? firstEvent() : 0;
  int levt = (lastEvent() > -1) ? lastEvent() : nEvents();
  cout << ">>> Event range: [" << fevt << ", " << levt -1 << "]" << endl;
  int nEventSel = 0;
  evtWeightSum_ = 0;
  for (int ev = fevt; ev < levt; ++ev) {
    clearEvent(); // reset tree variables 
    clearLists(); // reset analysis related lists for each event
    int lflag = chain()->LoadTree(ev);
    int nbytes = getEntry(lflag);    // returns total bytes read
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));

    histf()->cd();
        
    // For data or for MC without pileup
    double puevWt = 1;
    if (isMC() && usePUWt()) {
      int npu = 0;
      puevWt = wtPileUp(npu);
    }

    // Show status of the run
    const Event& evt = eventColl()->at(0);
    unsigned long run   = evt.run;
    unsigned long event = evt.event;
    unsigned long lumis = evt.lumis;
    
    // Show status of the run
    if (currentFile != lastFile)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << chain()->GetCurrentFile()->GetName()
	   << " <<< Run# " << setw(8) << run
	   << " Lumis# " << setw(6) << lumis
	   << " Event# " << setw(12) << event << " >>> "
	   << " Events proc. " << setw(9) << ev << "(of " << setw(9) << levt-1 << ")"
	   << endl;
    lastFile = currentFile;
    
    // Show the status
    if (ev%nPrint == 0 || firstEvent() > -1)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << currentFile
	   << " <<< Run# " << setw(8) << run
	   << " Lumis# " << setw(6) << lumis
	   << " Event# " << setw(12) << event << " >>> "
	   << " Events proc. " << setw(8) << ((firstEvent() > -1) ? ev - firstEvent() : ev)
	   << endl;
    
    // Select a set of events by [run, event]
    if (useEventList_ && eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      if (eventIdMap().find(mkey.str()) != eventIdMap().end()) continue;
    }
    
    histf()->cd();
    histf()->cd("MultiLeptonAnalysis");
    if (isMC()) {
      AnaUtil::fillHist1D("puweight", puevWt);

      // only when exclusive jets datasets are used (for np = 0)
      if (0) std::cout << "== nMEPartons: " << genEventColl()->at(0).nMEPartons 
		       << ", lheNOutPartons: " << genEventColl()->at(0).lheNOutPartons 
		       << std::endl;
      if (selectPM_ && genEventColl()->at(0).nMEPartons != nMEPartons_) continue;

      double wt = genEventColl()->size() ? genEventColl()->at(0).evtWeight : 1;
      AnaUtil::fillHist1D("evtweight", wt);      
      evtWeightSum_ += wt;   // this is used for final normalization
    }

    AnaUtil::fillHist1D("evtCutFlow", 0);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, puevWt);
    
    if (genAna_ != nullptr) {
      genAna_->setEvent(genParticleColl());
      bool genOk = genAna_->filter();
      if (!genOk) continue;
    }

    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, puevWt);

    histf()->cd();
    histf()->cd("MultiLeptonAnalysis");
    
    // at least 1 good PV
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    size_t ngoodVtx = vtxList_.size();
    if (ngoodVtx < 1) continue;

#ifdef SKIP_DUPLICATE_ALL
    // Duplicate event removal
    if (!isMC()) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      std::string evs(mkey.str());
      if (eventIdStore_.find(evs) != eventIdStore_.end()) {
	if (0) cout << "DuplicateAll: " << evs << endl;
	continue;
      }
      else {
	eventIdStore_.insert({evs, 1});
      }
    }
#endif

    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, puevWt);
    
    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt); 
    
    // Is event triggered?
    if (0) dumpTriggerPaths(std::cout, true);
    if (useTrigger() && !isTriggered(true, false)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, puevWt);
    
    // Main analysis object selection
    double vz = (vtxList_.size() > 0) ? vtxList_.at(0).z : -999;
    findObjects(vz, 1); // do not use event weight
    
    if (genAna_ != nullptr && dumpGenInfo_) genAna_->dumpEvent();
    //if (getFSRPhotonList().size()) dumpEvent(vz, false, true);

    // Access selected objects 
    const auto& elePhotonPairList = getTightIsoElePhotonPairList();
    int nEle = elePhotonPairList.size();

    const auto& muPhotonPairList  = getTightIsoMuPhotonPairList();
    int nMuon = muPhotonPairList.size();

    const auto& tauList = getIsoTauList();
    int nTau = tauList.size();
    
    histf()->cd();
    histf()->cd("MultiLeptonAnalysis");
    AnaUtil::fillHist1D("nGoodMuon", nMuon, puevWt);
    AnaUtil::fillHist1D("nGoodEle", nEle, puevWt);

    // Require at least 2 leptons of the same flavor (electron/muon)
    if (nMuon < 2 && nEle < 2) continue;

    // Highest pT lepton pT > 20
    // leptons coming from the h1 should have higher pT 
    double hLepPt = (nMuon > 0) ? muPhotonPairList[0].first.pt : 0.0;
    if (nEle && elePhotonPairList[0].first.pt > hLepPt) 
      hLepPt = elePhotonPairList[0].first.pt;
    AnaUtil::fillHist1D("hLepPt_stage0", hLepPt, puevWt);
    
    if (hLepPt < AnaUtil::cutValue(evselCutMap(), "hLepPtMin")) continue;
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, puevWt);

    // Find Z candidates
    std::vector<ZCandidate> ZCandList;
    if (nMuon >= 2) ZSelector(muPhotonPairList, ZCandList);
    if (nEle >= 2)  ZSelector(elePhotonPairList, ZCandList);
    AnaUtil::fillHist1D("nZcand", ZCandList.size(), puevWt);

    if (ZCandList.empty())                     continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, puevWt);

    // Now, add lepton isolation to the Z candidates found above
    addLeptonIsolation(elePhotonPairList, muPhotonPairList, ZCandList);

    // Z candidates
    if (ZCandList.size() > 1) 
      std::sort(std::begin(ZCandList), std::end(ZCandList), dmComparator);
    for (const auto& z: ZCandList) 
      AnaUtil::fillHist1D("massZcand", z.mass, puevWt);

    // The first Z candidate
    const ZCandidate& ZCand = ZCandList[0];
    if (ZCand.mass < AnaUtil::cutValue(evselCutMap(), "ZMassLow") || 
	ZCand.mass > AnaUtil::cutValue(evselCutMap(), "ZMassHigh")) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, puevWt);
    AnaUtil::fillHist1D("nvtx_stage0", ngoodVtx, puevWt);
    
#ifdef SKIP_DUPLICATE_ZMASS
    // Duplicate event removal
    if (!isMC()) {
      if (skipDuplicate_) { 
	std::ostringstream mkey;
	mkey << run << "-" << lumis << "-" << event;
	if (eventMap_.find(mkey.str()) != eventMap_.end()) {
	  if (0) cout << "DuplicateZMass: " << mkey.str() << endl;
	  continue;
	}
      }
      evLog() << run << " " << lumis << " " << event << std::endl;
    }
#endif
    // The second Z candidate
    if (ZCandList.size() > 1) AnaUtil::fillHist1D("Z2mass", ZCandList[1].mass, puevWt);

    // Now require at least four tight isolated leptons
    AnaUtil::fillHist1D("nTightIsoLeptons", nMuon + nEle, puevWt);    
    if (nMuon + nEle < 4) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, puevWt);
    AnaUtil::fillHist1D("nvtx_stage1", ngoodVtx, puevWt);

    // Go closer to final topology
    if (nMuon < 3 && nEle < 3) continue;
    if (nMuon < 1 || nEle < 1) continue;

    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, puevWt);
    AnaUtil::fillHist1D("nvtx_stage2", ngoodVtx, puevWt);

    // Events should not have any extra loose muons/electrons (ignore those which are inside the jets)
    AnaUtil::fillHist1D("nLooseLeptons", getLooseMuList().size() + getLooseEleList().size(), puevWt);    

    int nLooseLeptons = getTrueLooseLeptons();
    if (nLooseLeptons > 4) continue;
    AnaUtil::fillHist1D("nvtx_stage3", ngoodVtx, puevWt);
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, puevWt);

    // No tau should be present    
    AnaUtil::fillHist1D("nGoodTau", nTau, puevWt);
    if (nTau > 0)                         continue;
    AnaUtil::fillHist1D("nvtx_stage4", ngoodVtx, puevWt);
    AnaUtil::fillHist1D("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 10, puevWt);

    // No b-jets
    AnaUtil::fillHist1D("nJets", nTightJets(), puevWt);
    AnaUtil::fillHist1D("nbJets", nbJets(), puevWt);
    if (nbJets() > 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 11);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 11, puevWt);
    AnaUtil::fillHist1D("nvtx_stage5", ngoodVtx, puevWt);

    if (nLooseJets() > 0 && looseJetList().at(0).pt > AnaUtil::cutValue(evselCutMap(), "maxJetPt")) {
#if 0
      if (++nEventSel < dumpEventCount_) dumpEvent(vz, false, true);
#endif
      continue;
    }
    AnaUtil::fillHist1D("evtCutFlow", 12);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12, puevWt);
    AnaUtil::fillHist1D("nvtx_stage6", ngoodVtx, puevWt);

    AnaUtil::fillHist1D("met", metColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("corrmet", corrmetColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("puppimet", puppimetColl()->at(0).met, puevWt);

    if (corrmetColl()->at(0).met > AnaUtil::cutValue(evselCutMap(), "maxMET")) continue;
    AnaUtil::fillHist1D("evtCutFlow", 13);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 13, puevWt);
    AnaUtil::fillHist1D("nvtx_stage7", ngoodVtx, puevWt);

    // Now find the ele-mu candidate
    std::vector<ZCandidate> leptonPairCandList;
    auto result = leptonPairSelector(elePhotonPairList, muPhotonPairList, ZCand, leptonPairCandList, -1.0, false);
    if (0) std::cout << "result: " << result << std::endl;
    if (leptonPairCandList.empty())                continue;
    AnaUtil::fillHist1D("evtCutFlow", 14);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 14, puevWt);
    AnaUtil::fillHist1D("nvtx_stage8", ngoodVtx, puevWt);

    // Add lepton isolation to the e-mu candidates found above
    addLeptonIsolationForLL(elePhotonPairList, muPhotonPairList, leptonPairCandList);

    // Find an algorithm to sort
    if (leptonPairCandList.size() > 1) 
      std::sort(std::begin(leptonPairCandList), std::end(leptonPairCandList), massComparator);
    const ZCandidate& h1Cand = leptonPairCandList[0];

    AnaUtil::fillHist1D("hLepPt", hLepPt, puevWt);
    AnaUtil::fillHist1D("Zmass", ZCand.mass, puevWt);
    AnaUtil::fillHist1D("ZPt", ZCand.p4.Pt(), puevWt);
    AnaUtil::fillHist1D("h1mass", h1Cand.mass, puevWt);
    AnaUtil::fillHist1D("h1Pt", h1Cand.p4.Pt(), puevWt);
    AnaUtil::fillHist1D("dRh1Z", h1Cand.p4.DeltaR(ZCand.p4), puevWt);
    AnaUtil::fillHist1D("dEtall", ZCand.dEtall, puevWt);
    AnaUtil::fillHist1D("dPhill", ZCand.dPhill, puevWt);
    AnaUtil::fillHist1D("dRll", ZCand.dRll, puevWt);
    AnaUtil::fillHist1D("dEtaemu", h1Cand.dEtall, puevWt);
    AnaUtil::fillHist1D("dPhiemu", h1Cand.dPhill, puevWt);
    AnaUtil::fillHist1D("dRemu", h1Cand.dRll, puevWt);

    AnaUtil::fillHist2D("h1VsZpT", ZCand.p4.Pt(), h1Cand.p4.Pt(), puevWt);

    // Lepton properties
    const TLorentzVector& lep1P4 = ZCand.l1RP4;
    const TLorentzVector& lep2P4 = ZCand.l2RP4;
    const TLorentzVector& eleP4  = h1Cand.l1RP4;
    const TLorentzVector& muP4   = h1Cand.l2RP4;

    AnaUtil::fillHist1D("Zl1Pt", lep1P4.Pt(), puevWt);
    AnaUtil::fillHist1D("Zl2Pt", lep2P4.Pt(), puevWt);
    AnaUtil::fillHist1D("h1ElePt", eleP4.Pt(), puevWt);
    AnaUtil::fillHist1D("h1MuPt", muP4.Pt(), puevWt);

    AnaUtil::fillHist1D("Zl1Eta", lep1P4.Eta(), puevWt);
    AnaUtil::fillHist1D("Zl2Eta", lep2P4.Eta(), puevWt);
    AnaUtil::fillHist1D("h1EleEta", eleP4.Eta(), puevWt);
    AnaUtil::fillHist1D("h1MuEta", muP4.Eta(), puevWt);
    
    AnaUtil::fillHist1D("Zl1Phi", lep1P4.Phi(), puevWt);
    AnaUtil::fillHist1D("Zl2Phi", lep2P4.Phi(), puevWt);
    AnaUtil::fillHist1D("h1ElePhi", eleP4.Phi(), puevWt);
    AnaUtil::fillHist1D("h1MuPhi", muP4.Phi(), puevWt);

    AnaUtil::fillHist2D("lep2VsLep1pT", lep1P4.Pt(), lep2P4.Pt(), puevWt);
    AnaUtil::fillHist2D("muVsLep1pT", lep1P4.Pt(), muP4.Pt(), puevWt);
    AnaUtil::fillHist2D("muVsLep2pT", lep2P4.Pt(), muP4.Pt(), puevWt);
    AnaUtil::fillHist2D("muVsElepT", eleP4.Pt(),  muP4.Pt(), puevWt);

    AnaUtil::fillHist1D("h1PtRatio", h1Cand.p4.Pt()/(eleP4.Pt()+muP4.Pt()), puevWt);

    TVector3 lep1_Alpha(lep1P4.Px(), lep1P4.Py(), lep1P4.Pz());
    TVector3 lep2_Alpha(lep2P4.Px(), lep2P4.Py(), lep2P4.Pz());
    TVector3 Z_Alpha(0, 0, 1.0);

    TVector3 d1(eleP4.Px(), eleP4.Py(), eleP4.Pz());
    TVector3 d2(muP4.Px(), muP4.Py(), muP4.Pz());
    TVector3 h1CrossProd = d1.Cross(d2);

    TVector3 lep1Z = Z_Alpha.Cross(lep1_Alpha);
    double angle = lep1Z.Angle(h1CrossProd);
    AnaUtil::fillHist1D("h1Lep1Angle", angle, puevWt);

    TVector3 lep2Z = Z_Alpha.Cross(lep2_Alpha);
    angle = lep2Z.Angle(h1CrossProd);
    AnaUtil::fillHist1D("h1Lep2Angle", angle, puevWt);

    // Scalar eT sum of leptons
    double lt = lep1P4.Et() + lep2P4.Et() + eleP4.Et() + muP4.Et(); 
    AnaUtil::fillHist1D("sLTSum", lt, puevWt);

    // Vector eT sum of leptons
    TLorentzVector vv = lep1P4 + lep2P4 + eleP4 + muP4;
    AnaUtil::fillHist1D("vLTSum", vv.Et(), puevWt);

    // Event type
    int type = static_cast<int>(setEventType(ZCand, h1Cand));
    AnaUtil::fillHist1D("nEventsPerType", type);
    if (isMC()) AnaUtil::fillHist1D("nEventsPerTypeWt", type, puevWt);

    // fill histograms for different event categories (e.g mmem, eeem etc.)
    std::ostringstream zcat;
    zcat << "ZmassC" << setw(1) << type;
    AnaUtil::fillHist1D(zcat.str(), ZCand.mass, puevWt);

    std::ostringstream hcat;
    hcat << "h1massC" << setw(1) << type;
    AnaUtil::fillHist1D(hcat.str(), h1Cand.mass, puevWt);

    std::ostringstream drcat;
    drcat << "dRh1ZC" << setw(1) << type;
    AnaUtil::fillHist1D(drcat.str(), h1Cand.p4.DeltaR(ZCand.p4), puevWt);

    selEvLog() << run << " " << lumis << " " << event << std::endl;

    // Print only the first n events; n configurable
    if (isMC() && dumpEventCount_ > -1 && ++nEventSel >= dumpEventCount_) continue;

    cout << ">>> "
         << "<nLooseMuon>: " << getLooseMuList().size()
         << ", <nMuon>: " << muPhotonPairList.size()
         << ", <nLooseElectron>: " << getLooseEleList().size()
         << ", <nElectron>: " << elePhotonPairList.size()
         << ", <nTau>: " << getTauList().size()
         << ", <nIsoTau>: " << nTau
         << ", <nLooseJet>: " << looseJetList().size()
         << ", <nTightJet>: " << tightJetList().size()
         << ", <nZCandidate>: " << ZCandList.size()
         << ", <nemuCandidate>: " << leptonPairCandList.size()
	 << endl;

    // Print Z Candidates
    dumpEvent(vz, false, true);
    for (const auto& c: ZCandList) 
      HZZ4lUtil::printZCandidate(c, "Z Candidate");

    // Print e-mu candidates
    for (const auto& c: leptonPairCandList)
      HZZ4lUtil::printZCandidate(c, "e-mu Candidate", false);
  }
}
MultiLeptonAnalysis::EventType MultiLeptonAnalysis::setEventType(const ZCandidate& ZCand, const ZCandidate& emuCand) {
  EventType type = EventType::unkwn;
  if      (ZCand.flavour == HZZ4lUtil::ZType::mumu && 
	   (emuCand.flavour == HZZ4lUtil::llType::emu || emuCand.flavour == HZZ4lUtil::llType::mue)) type = EventType::mmem;
  else if (ZCand.flavour == HZZ4lUtil::ZType::ee && 
	   (emuCand.flavour == HZZ4lUtil::llType::emu || emuCand.flavour == HZZ4lUtil::llType::mue)) type = EventType::eeem;
  
  return type;
}
// Function to calculate number of extra tight leptons passing isolation apart from leptons in the Z-h1 candidate
int MultiLeptonAnalysis::findExtraLeptons(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand) {
  int nExtLep = 0;
  for (const auto& mu: getTightIsoMuList()) {
    TLorentzVector muP4(HZZ4lUtil::getP4(mu));
    if (AnaUtil::sameObject(muP4, Z1Cand.l1P4) || 
	AnaUtil::sameObject(muP4, Z1Cand.l2P4) || 
	AnaUtil::sameObject(muP4, Z2Cand.l1P4) || 
	AnaUtil::sameObject(muP4, Z2Cand.l2P4)) continue;
    nExtLep++;
  }
  for (const auto& ele: getTightIsoEleList()) {
    TLorentzVector eleP4(HZZ4lUtil::getP4(ele));
    if (AnaUtil::sameObject(eleP4, Z1Cand.l1P4) || 
	AnaUtil::sameObject(eleP4, Z1Cand.l2P4) || 
	AnaUtil::sameObject(eleP4, Z2Cand.l1P4) || 
	AnaUtil::sameObject(eleP4, Z2Cand.l2P4)) continue;
    nExtLep++;
  }
  return nExtLep;
}
int MultiLeptonAnalysis::findEventCategory(int nleptons, 
					   const std::vector<vhtm::Jet>& jetList, 
					   int nbjets,
					   const ZCandidate& Z1Cand, 
					   const ZCandidate& Z2Cand, 
					   bool verbose) {
  int njets = jetList.size();

  TLorentzVector j1P4, j2P4;
  j1P4.SetPtEtaPhiE(0.,0.,0.,0.);
  j2P4.SetPtEtaPhiE(0.,0.,0.,0.);
  if (njets) j1P4 = HZZ4lUtil::getP4(jetList[0]);
  if (njets > 1) j2P4 = HZZ4lUtil::getP4(jetList[1]);

  double djet = -1., mjj = 0.;
  TLorentzVector final4lP4 = Z1Cand.l1P4 + Z1Cand.l2P4 + Z1Cand.fsrPhoP4 + Z2Cand.l1P4 + Z2Cand.l2P4 + Z2Cand.fsrPhoP4;
  if (njets >= 2) {
    mjj = (j1P4 + j2P4).M();
    djet = 0.18 * std::fabs(j1P4.Eta() - j2P4.Eta()) + 1.92E-04 * mjj;
  }
  histf()->cd();
  histf()->cd("MultiLeptonAnalysis");
  AnaUtil::fillHist1D("djet", djet, 1);

  int cat = 0;
  if (nleptons == 4 && njets >= 2 && nbjets <= 1 && djet > 0.5) 
    cat = 2;
  else if ( (nleptons == 4 && 
        njets >= 2 && MultiLeptonAnalysis::hasJetPair(jetList) &&
        final4lP4.Pt() > final4lP4.M() ) ||
      (nleptons == 4 && njets == 2 && nbjets == 2) )
    cat = 4;
  else if (njets <= 2 && nbjets == 0 && nleptons >= 5)
    cat = 3;
  else if ( (njets >= 3 && nbjets >= 1) || nleptons >= 5 )
    cat = 5;
  else if (njets >= 1)
    cat = 1;

  if (verbose) {
    cout << "---- Event Category" << endl;
    cout << "  nlep  njet nbjet   jet1Pt  jet1Eta   jet2Pt  jet2Eta      mjj     4lPt      4lM     djet category"
      << endl;
    cout << setw(6) << nleptons
      << setw(6) << njets
      << setw(6) << nbjets
      << setw(9) << j1P4.Pt()
      << setw(9) << (j1P4.Pt() > 0 ? j1P4.Eta() : 99)
      << setw(9) << j2P4.Pt()
      << setw(9) << (j2P4.Pt() > 0 ? j2P4.Eta() : 99)
      << setw(9) << mjj
      << setw(9) << final4lP4.Pt() 
      << setw(9) << final4lP4.M()
      << setw(9) << djet
      << setw(9) << cat
      << endl << endl;
  }
  return cat;
}
bool MultiLeptonAnalysis::hasJetPair(const std::vector<vhtm::Jet>& jetList) {
  for (size_t i = 0; i < jetList.size(); ++i) {
    const auto& j1 = jetList[i];
    TLorentzVector j1P4(HZZ4lUtil::getP4(j1));
    for (size_t j = i+1; j < jetList.size(); ++j) {
      const auto& j2 = jetList[j];
      TLorentzVector j2P4(HZZ4lUtil::getP4(j2));
      double mjj = (j1P4 + j2P4).M();      
      if (0) cout << "mjj[" << i << ", " << j << "] = " << mjj << endl;
      if ( (std::fabs(j1.eta) < 2.4 && j1.pt > 40.) && 
	   (std::fabs(j2.eta) < 2.4 && j2.pt > 40.) && 
	   (mjj > 60. && mjj < 120.) ) return true;
    }
  }
  return false;
}
void MultiLeptonAnalysis::endJob() {
  PhysicsObjSelector::endJob();

  histf()->cd();
  histf()->cd("MultiLeptonAnalysis");
  vector<string> evLabels {
    "Events processed",
    "Gen Filter",
    "Events with > 0 good vertex",
    "Events passing trigger",
    "# leptons >= 2",
    "# of Z Candidates > 0",
    "mlow < Z mass < mhigh",
    "# leptons >= 4",
    "at least 1 ele/mu",
    "# extra leptons == 0",
    "# tau == 0",
    "# b-jets == 0",
    "# loose-jets == 0",
    "MET < 60 GeV",
    "has a (Z + e-mu) candidate"
  };
  HZZ4lUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection");  

  double lumiFac = 1.0;
  if (isMC()) {
    lumiFac = lumiWt(evtWeightSum_);
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
  }

  if (isMC()) {
    HZZ4lUtil::scaleHistogram("evtCutFlowWt", lumiFac);
    HZZ4lUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");  
  }
  
  // Event category
  vector<string> evCatList {
    "unknwn",
    "mmem",
    "eeem"
  };
  HZZ4lUtil::showCount("nEventsPerType", evCatList, "Events Per type");  
  if (isMC()) {
    HZZ4lUtil::scaleHistogram("nEventsPerTypeWt", lumiFac);
    HZZ4lUtil::showCount("nEventsPerTypeWt", evCatList, "Events Per type (Weighted)", 3);  
  }

  syncDumpf_.close();
}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default
// -------------------------------------------------------------------------------
bool MultiLeptonAnalysis::readJob(const string& jobFile, int& nFiles)
{
  if (!PhysicsObjSelector::readJob(jobFile, nFiles)) return false;
  
  // Open the same file containing the datacards again to read analysis specific cards
  ifstream fin(jobFile.c_str(), std::ios::in);    
  if (!fin) {
    cerr << "==> Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }
  
  eventFilelist_.clear();  

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   
    
    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;
    
    // Split the line into words
    vector<string> tokens;
    AnaUtil::tokenize(line, tokens);
    assert(tokens.size() > 1);
    const string& key   = tokens[0];
    const string& value = tokens[1];
    if (key == "useEventList")
      useEventList_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "skipDuplicate")
      skipDuplicate_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "inputEventFile")
      eventFilelist_.push_back(value.c_str());
    else if (key == "syncDumpFile")
      dumpFilename_ = value.c_str();
    else if (key == "dumpEventMax")
      dumpEventCount_ = std::stoi(value.c_str());
    else if (key == "selectPartons")
      selectPM_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "nMEPartons")
      nMEPartons_ = std::stoi(value.c_str());
  }
  // Close the file
  fin.close();
  
  if (!dumpFilename_.empty()) {
    syncDumpf_.open(dumpFilename_, std::ios::out);
    if (!syncDumpf_) {
      cerr << "Output File: " << dumpFilename_ << " could not be opened!" << endl;
      return false;
    }
  }  
#ifdef SKIP_DUPLICATE_ZMASS
  if (skipDuplicate_ && !eventFilelist_.empty()) {
    eventMap_.clear();
    for (const auto& f: eventFilelist_) {
      cout << ">>> Reading file: " << f << endl;
      ifstream fin(f, std::ios::in);
      if (!fin) {
	cerr << "Input file: " << f << " could not be opened!" << endl;
	continue;
      }
      char buf[BUF_SIZE];
      vector<string> tokens;
      while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
	string line(buf);
	if (line.empty()) continue;   
	
	// enable '#' and '//' style comments
	if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    
	// Split the line into words
	AnaUtil::tokenize(line, tokens);
	assert(tokens.size() > 2);
	string key = tokens.at(0) + "-" + tokens.at(1) + "-" + tokens.at(2);
	eventMap_.insert({key, 1});

	tokens.clear();
      }
      // Close the file
      fin.close();
    }
    cout << ">>> Total events present: " << eventMap_.size() << endl;
  }  
#endif
  printJob();
  
  if (readGenInfo()) genAna_ = std::make_unique<GenAnalysis>();
  return true;
}
void MultiLeptonAnalysis::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  os << endl;
  os << "   useEventList: " << std::boolalpha << useEventList_ << endl
     << "  skipDuplicate: " << std::boolalpha << skipDuplicate_ << endl
     << " dumpEventCount: " << dumpEventCount_ << endl
     << "   syncDumpFile: " << dumpFilename_ << endl
     << "   dumpEventMax: " << dumpEventCount_ << endl
     << "  selectPartons: " << std::boolalpha << selectPM_ << endl
     << "     nMEPartons: " << nMEPartons_ << endl;
  if (isMC()) 
  os << "    dumpGenInfo: " << std::boolalpha << dumpGenInfo_ << endl;

  AnaUtil::showList(eventFilelist_, ">>> INFO. Input event files:", os);
}
