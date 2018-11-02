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
#include "LeptonPlusTauAnalysis.h"
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
LeptonPlusTauAnalysis::LeptonPlusTauAnalysis()
  : PhysicsObjSelector(),
  dumpGenInfo_(false),
  useEventList_(false),
  skipDuplicate_(false), 
  selectEvType_(false),
  evtype_(-1),
  dumpFilename_("syncDumpFile.txt"),
  dumpEventCount_(0)
{
  genAna_ = (readGenInfo()) ? new GenAnalysis() : nullptr;
}
// ----------
// Destructor
// ----------
LeptonPlusTauAnalysis::~LeptonPlusTauAnalysis() 
{
  if (genAna_) delete genAna_;
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool LeptonPlusTauAnalysis::beginJob() 
{ 
  // Open the output ROOT file (in AnaBase)
  PhysicsObjSelector::beginJob();

  histf()->cd();
  histf()->mkdir("LeptonPlusTauAnalysis");
  
  bookHistograms();
  if (genAna_) genAna_->bookHistograms(histf());

  return true;
}
// ---------------
// Book histograms
// ---------------
void LeptonPlusTauAnalysis::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("LeptonPlusTauAnalysis");

  new TH1D("nvtx", "Number of Good vertices", 60, 0, 60);
  if (isMC()) {
    new TH1D("puweight", "PU reweight factor", 100, 0, 2);
    new TH1D("evtweight", "Event weight factor (MC events)", 100, -5, 5);
  }

  //------- Object PLots -----------------------------------------------
  new TH1D("nGoodMuon", "Number of Good muons", 6, -0.5, 5.5);
  new TH1D("nGoodEle", "Number of Good electrons", 6, -0.5, 5.5);
  new TH1D("nGoodTau", "Number of Good Tau", 5, -0.5, 4.5);
  new TH1D("nJets", "Number of Tight Jets", 10, -0.5, 9.5);
  new TH1D("nbJets", "Number of b-Jets", 10, -0.5, 9.5);
  new TH1D("evtCutFlow", "Event CutFlow", 10, -0.5, 9.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 10, -0.5, 9.5);
  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);

  new TH1D("hLepPt", "Highest pT lepton pT", 300, 0, 300);
  new TH1D("met", "Missing Transver Energy", 150, 0, 150);
  new TH1D("corrmet", "Missing Transver Energy (Corrected)", 150, 0, 150);
  new TH1D("puppimet", "Missing Transver Energy (PUPPI)", 150, 0, 150);

  // Z and h1
  new TH1F("nZcand", "Number of selected Z candidates", 10, 0, 10);
  new TH1F("massZcand", "Mass of selected Z candidates", 100, 0., 200.);

  new TH1F("Zmass",   "ll invariant mass of the Z candidate", 200, 0., 200.);
  new TH1F("Z2mass",  "ll invariant mass of the second Z candidate", 200, 0., 200.);
  new TH1F("ZPt",     "pT of the ll system of the Z candidate", 200, 0., 400.);
  new TH1F("ZmassC0", "ll invariant mass of the Z candidate (mmmtau)", 200, 0., 200.);
  new TH1F("ZmassC1", "ll invariant mass of the Z candidate (mmetau)", 200, 0., 200.);
  new TH1F("ZmassC2", "ll invariant mass of the Z candidate (eeetau)", 200, 0., 200.);
  new TH1F("ZmassC3", "ll invariant mass of the Z candidate (eemtau)", 200, 0., 200.);

  new TH1F("h1mass",   "mass of the selected Higgs candidate", 300, 0, 300);
  new TH1F("h1Pt",     "pT of the ltau system", 200, 0., 400.);
  new TH1F("h1massC0", "Mass of the selected Higgs candidate (mmmtau)", 300, 0, 300);
  new TH1F("h1massC1", "Mass of the selected Higgs candidate (mmetau)", 300, 0, 300);
  new TH1F("h1massC2", "Mass of the selected Higgs candidate (eeetau)", 300, 0, 300);
  new TH1F("h1massC3", "Mass of the selected Higgs candidate (eemtau)", 300, 0, 300);

  // angle
  new TH1F("dRh1Z",   "dR(h1, Z)", 100, 0, 2 * M_PI);
  new TH1F("dRh1ZC0", "dR(h1, Z) (mmmtau)", 100, 0, 2 * M_PI);
  new TH1F("dRh1ZC1", "dR(h1, Z) (mmetau)", 100, 0, 2 * M_PI);
  new TH1F("dRh1ZC2", "dR(h1, Z) (eeetau)", 100, 0, 2 * M_PI);
  new TH1F("dRh1ZC3", "dR(h1, Z) (eemtau)", 100, 0, 2 * M_PI);

  new TH1F("dRll",    "dR(lep, lep) of Z", 100, 0, 2 * M_PI);
  new TH1F("dRltau",  "dR(lep, tau) of h1 candidate", 100, 0, 2 * M_PI);

  // Lepton properties (pT, eta etc.)
  new TH1F("lep1Pt", "pT of the fist lepton from Z", 300, 0, 300);
  new TH1F("lep2Pt", "pT of the second lepton from Z", 300, 0, 300);
  new TH1F("lep3Pt", "pT of the first lepton from the h1 candidate", 300, 0, 300);
  new TH1F("tauPt",  "pT of the tau lepton from the h1 candidate", 300, 0, 300);

  new TH1F("lep1Eta", "Eta of the first lepton from Z", 200, -4., 4);
  new TH1F("lep2Eta", "Eta of the second lepton from Z", 200, -4, 4);
  new TH1F("lep3Eta", "Eta of the lepton from the h1 candidate", 200, -4, 4);
  new TH1F("tauEta",  "Eta of the tau lepton from the h1 candidate", 100, -4, 4);

  new TH1F("h1PtRatio", "pT(h1+lepton) over pT(h1)+pT(lepton)", 100, 0, 2);
  new TH1F("h1Lep1Angle", "Angle of lepton 1 from Z with the ltau decay plane", 100, 0, M_PI);
  new TH1F("h1Lep2Angle", "Angle of lepton 2 from Z with the ltau decay plane", 100, 0, M_PI);

  new TH1F("nEventsPerType", "Selected events per category", 5, -1.5, 3.5);
  if (isMC()) new TH1F("nEventsPerTypeWt", "Selected events per category (Weighted)", 5, -1.5, 3.5);

  new TH1F("sLTSum", "Scalar sum of Transverse Energy of the leptons", 200, 0, 800);
  new TH1F("vLTSum", "Vector sum of Transverse Energy of the leptons", 100, 0, 500);

  new TH2F("h1VsZpT", "h1 vs Z pT scatter plot", 100, 0, 400, 100, 0, 400);
  new TH2F("lep2Vslep1pT", "lepton 1 vs lepton 1 pT scatter plot", 100, 0, 300, 100, 0, 300);
  new TH2F("tauVslep1pT", "tau vs lepton 1 pT scatter plot", 100, 0, 300, 100, 0, 300);
  new TH2F("tauVslep2pT", "tau vs lepton 2 pT scatter plot", 100, 0, 300, 100, 0, 300);
  new TH2F("tauVslep3pT", "tau vs lepton 3 pT scatter plot", 100, 0, 300, 100, 0, 300);

  histf()->cd();
  histf()->ls();
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void LeptonPlusTauAnalysis::clearLists() {
  PhysicsObjSelector::clear();
  vtxList_.clear();
  genZList_.clear();
  ZCandList_.clear();
  lepTauCandList_.clear();
  evtype_ = -1;
}
// -------------------
// The main event loop
// -------------------
void LeptonPlusTauAnalysis::eventLoop()
{
  // Initialize analysis
  if (!beginJob()) return;
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
    clearLists(); // reset analysis related lists  for each event
    int lflag = chain()->LoadTree(ev);
    int nbytes = getEntry(lflag);    // returns total bytes read
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));
    const Event& evt = eventColl()->at(0);

    histf()->cd();
        
    // For data or for MC without pileup
    puevWt_ = 1;
    if (isMC() && usePUWt()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }

    // Show status of the run
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
      if (eventIdMap().count(mkey.str()) == 0) continue;
      //if (eventIdMap().find(mkey.str()) != eventIdMap().end()) continue;
    }
    
    histf()->cd();
    histf()->cd("LeptonPlusTauAnalysis");
    if (isMC()) {
#ifdef __MC__
      double wt = genEventColl()->at(0).evtWeight;
#else
      double wt = 1;
#endif
      AnaUtil::fillHist1D("evtweight", wt);      
      evtWeightSum_ += wt;

      AnaUtil::fillHist1D("puweight", puevWt_);
    }

    AnaUtil::fillHist1D("evtCutFlow", 0);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, puevWt_);
    
    if (genAna_) {
      genAna_->setEvent(genParticleColl());
      bool genOk = genAna_->filter();
      if (!genOk) continue;
    }
    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, puevWt_);

    histf()->cd();
    histf()->cd("LeptonPlusTauAnalysis");
    
    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt_); 
    
    // is event triggered?
    if (useTrigger() && !isTriggered(true, false)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, puevWt_);
    
    // at least 1 good PV
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    size_t ngoodVtx = vtxList_.size();
    AnaUtil::fillHist1D("nvtx", ngoodVtx, puevWt_);
    
    if (ngoodVtx < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, puevWt_);
    
    // main analysis object selection
    double vz = (vtxList_.size() > 0) ? vtxList_.at(0).z : -999;
    findObjects(vz, 1); // do not use event weight
    
    if (genAna_ && dumpGenInfo_) genAna_->dumpEvent();
    if (0) dumpEvent(vz, false, true);

    // access selected objects 
    const auto& elePhotonPairList = getTightIsoElePhotonPairList();
    int nEle = elePhotonPairList.size();

    const auto& muPhotonPairList  = getTightIsoMuPhotonPairList();
    int nMuon = muPhotonPairList.size();

    const auto& tauList = getIsoTauList();
    int nTau = tauList.size();
    
    histf()->cd();
    histf()->cd("LeptonPlusTauAnalysis");
    AnaUtil::fillHist1D("nGoodMuon", nMuon, puevWt_);
    AnaUtil::fillHist1D("nGoodEle", nEle, puevWt_);

    // Require at least 2 leptons (electron/muon)
    if (nMuon < 2 && nEle < 2) continue;

    AnaUtil::fillHist1D("nGoodTau", nTau, puevWt_);
    AnaUtil::fillHist1D("nJets", getNTightJets(), puevWt_);
    AnaUtil::fillHist1D("nbJets", getNbJets(), puevWt_);

    // Highest pT lepton pT > 20
    double hLepPt = 0;
    if (nMuon) hLepPt = muPhotonPairList[0].first.pt;
    if (nEle && elePhotonPairList[0].first.pt > hLepPt) hLepPt = elePhotonPairList[0].first.pt;
    
    AnaUtil::fillHist1D("hLepPt", hLepPt, puevWt_);
    if (hLepPt < 20)           continue; // add in the job card
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, puevWt_);

    // find Z candidates
    if (nMuon >= 2) ZSelector<vhtm::Muon>(muPhotonPairList, ZCandList_);
    if (nEle >= 2)  ZSelector<vhtm::Electron>(elePhotonPairList, ZCandList_);
    AnaUtil::fillHist1D("nZcand", ZCandList_.size(), puevWt_);

    if (ZCandList_.empty())                     continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, puevWt_);

    if (eventMap_.size()) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      if (eventMap_.count(mkey.str()) > 0) {
	cout << "Duplicate: " << mkey.str() << endl;
	continue;
      }
    }

    // Now, add lepton isolation to the Z candidates found above
    addLeptonIsolation(elePhotonPairList, muPhotonPairList, ZCandList_);

    // The true Z candidate
    if (ZCandList_.size() > 1) 
      std::sort(ZCandList_.begin(), ZCandList_.end(), dmComparator);
    for (const auto& z: ZCandList_) 
      AnaUtil::fillHist1D("massZcand", z.mass, puevWt_);
    const ZCandidate& ZCand = ZCandList_[0];
    if (ZCandList_.size() > 1) AnaUtil::fillHist1D("Z2mass", ZCandList_[1].mass, puevWt_);

    // Require at least one tau    
    if (tauList.empty())                         continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, puevWt_);

    // Now require a third lepton
    if (nMuon + nEle < 3) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, puevWt_);

    // events should not have any extra loose muons/electrons (tau?)
    if (getLooseMuList().size() + getLooseEleList().size() > 3) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, puevWt_);
  
    // Now find the lepton-tau candidate
    lepTauSelector<vhtm::Muon>(muPhotonPairList, tauList, ZCand, lepTauCandList_);
    lepTauSelector<vhtm::Electron>(elePhotonPairList, tauList, ZCand, lepTauCandList_);
    
    if (lepTauCandList_.empty())                continue;
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, puevWt_);

    // add lepton isolation to the Z candidates found above
    addLeptonIsolationForLTau(elePhotonPairList, muPhotonPairList, lepTauCandList_);

    // find an algorithm to sort
    if (lepTauCandList_.size() > 1) 
      std::sort(lepTauCandList_.begin(), lepTauCandList_.end(), massComparator);
    const ZCandidate& h1Cand = lepTauCandList_[0];

    AnaUtil::fillHist1D("met", metColl()->at(0).met, puevWt_);
    AnaUtil::fillHist1D("corrmet", corrmetColl()->at(0).met, puevWt_);
    AnaUtil::fillHist1D("puppimet", puppimetColl()->at(0).met, puevWt_);

    AnaUtil::fillHist1D("Zmass", ZCand.mass, puevWt_);
    AnaUtil::fillHist1D("ZPt", ZCand.p4.Pt(), puevWt_);
    AnaUtil::fillHist1D("h1mass", h1Cand.mass, puevWt_);
    AnaUtil::fillHist1D("h1Pt", h1Cand.p4.Pt(), puevWt_);
    AnaUtil::fillHist1D("dRh1Z", h1Cand.p4.DeltaR(ZCand.p4), puevWt_);
    AnaUtil::fillHist1D("dRll", ZCand.dRll, puevWt_);
    AnaUtil::fillHist1D("dRltau", h1Cand.dRll, puevWt_);

    AnaUtil::fillHist2D("h1VsZpT", ZCand.p4.Pt(), h1Cand.p4.Pt(), puevWt_);

    // Lepton properties
    const TLorentzVector& lep1P4 = ZCand.l1RP4;
    const TLorentzVector& lep2P4 = ZCand.l2RP4;
    const TLorentzVector& lep3P4 = h1Cand.l1RP4;
    const TLorentzVector& tauP4  = h1Cand.l2RP4;

    AnaUtil::fillHist1D("lep1Pt", lep1P4.Pt(), puevWt_);
    AnaUtil::fillHist1D("lep2Pt", lep2P4.Pt(), puevWt_);
    AnaUtil::fillHist1D("lep3Pt", lep3P4.Pt(), puevWt_);
    AnaUtil::fillHist1D("tauPt",   tauP4.Pt(), puevWt_);

    AnaUtil::fillHist1D("lep1Eta", lep1P4.Eta(), puevWt_);
    AnaUtil::fillHist1D("lep2Eta", lep2P4.Eta(), puevWt_);
    AnaUtil::fillHist1D("lep3Eta", lep3P4.Eta(), puevWt_);
    AnaUtil::fillHist1D("tauEta",   tauP4.Eta(), puevWt_);
    
    AnaUtil::fillHist2D("lep2Vslep1pT", lep1P4.Pt(), lep2P4.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauVslep1pT", lep1P4.Pt(), tauP4.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauVslep2pT", lep2P4.Pt(), tauP4.Pt(), puevWt_);
    AnaUtil::fillHist2D("tauVslep3pT", lep3P4.Pt(), tauP4.Pt(), puevWt_);

    AnaUtil::fillHist1D("h1PtRatio", h1Cand.p4.Pt()/(lep3P4.Pt()+tauP4.Pt()), puevWt_);

    TVector3 lep1_Alpha(lep1P4.Px(), lep1P4.Py(), lep1P4.Pz());
    TVector3 lep2_Alpha(lep2P4.Px(), lep2P4.Py(), lep2P4.Pz());
    TVector3 Z_Alpha(0, 0, 1.0);

    TVector3 d1(lep3P4.Px(), lep3P4.Py(), lep3P4.Pz());
    TVector3 d2(tauP4.Px(), tauP4.Py(), tauP4.Pz());
    TVector3 h1CrossProd = d1.Cross(d2);

    TVector3 lep1Z = Z_Alpha.Cross(lep1_Alpha);
    double angle = lep1Z.Angle(h1CrossProd);
    AnaUtil::fillHist1D("h1Lep1Angle", angle, puevWt_);

    TVector3 lep2Z = Z_Alpha.Cross(lep2_Alpha);
    angle = lep2Z.Angle(h1CrossProd);
    AnaUtil::fillHist1D("h1Lep2Angle", angle, puevWt_);

    // Scalar eT sum of leptons
    double lt = lep1P4.Et() + lep2P4.Et() + lep3P4.Et() + tauP4.Et();
    AnaUtil::fillHist1D("sLTSum", lt, puevWt_);

    // Vector eT sum of leptons
    TLorentzVector vv = lep1P4 + lep2P4 + lep3P4 + tauP4;
    AnaUtil::fillHist1D("vLTSum", vv.Et(), puevWt_);

    // 
    int type = static_cast<int>(setEventType(ZCand, h1Cand));
    AnaUtil::fillHist1D("nEventsPerType", type);
    if (isMC()) AnaUtil::fillHist1D("nEventsPerTypeWt", type, puevWt_);

    // fill histograms for different event categories (e.m mmmtau, mmetau etc.)
    std::ostringstream zcat;
    zcat << "ZmassC" << setw(1) << type;
    AnaUtil::fillHist1D(zcat.str(), ZCand.mass, puevWt_);

    std::ostringstream hcat;
    hcat << "h1massC" << setw(1) << type;
    AnaUtil::fillHist1D(hcat.str(), h1Cand.mass, puevWt_);

    std::ostringstream drcat;
    drcat << "dRh1ZC" << setw(1) << type;
    AnaUtil::fillHist1D(drcat.str(), h1Cand.p4.DeltaR(ZCand.p4), puevWt_);

    evLog() << run << " " << lumis << " " << event << std::endl;

    // Print only the first n events; n configurable
    if (isMC() && dumpEventCount_ > -1 && ++nEventSel >= dumpEventCount_) continue;

    cout << ">>> "
         << "<nLooseMuon>: " << getLooseMuList().size()
         << ", <nMuon>: " << muPhotonPairList.size()
         << ", <nLooseElectron>: " << getLooseEleList().size()
         << ", <nElectron>: " << elePhotonPairList.size()
         << ", <nTau>: " << tauList.size()
         << ", <nIsoTau>: " << nTau
         << ", <nZCandidate>: " << ZCandList_.size()
         << ", <nltauCandidate>: " << lepTauCandList_.size()
	 << endl;

    // Print Z Candidates
    dumpEvent(vz, false, true);
    for (const auto& c: ZCandList_) 
      HZZ4lUtil::printZCandidate(c, "Z Candidate");

    // Print Mu-Tau candidates
    for (const auto& c: lepTauCandList_)
      HZZ4lUtil::printZCandidate(c, "Lepton-Tau Candidate", false);
  }
  // Analysis over
  endJob();
}
LeptonPlusTauAnalysis::EventType LeptonPlusTauAnalysis::setEventType(const ZCandidate& ZCand, const ZCandidate& lTauCand) {
  EventType type = EventType::unkwn;
  if      (ZCand.flavour == HZZ4lUtil::ZType::mumu && lTauCand.flavour == HZZ4lUtil::lTauType::mutau) type = EventType::mmmtau;
  else if (ZCand.flavour == HZZ4lUtil::ZType::mumu && lTauCand.flavour == HZZ4lUtil::lTauType::etau)  type = EventType::mmetau;
  else if (ZCand.flavour == HZZ4lUtil::ZType::ee   && lTauCand.flavour == HZZ4lUtil::lTauType::etau)  type = EventType::eeetau;
  else if (ZCand.flavour == HZZ4lUtil::ZType::ee   && lTauCand.flavour == HZZ4lUtil::lTauType::mutau) type = EventType::eemtau;
  
  return type;
}
// Function to calculate number of extra tight leptons passing isolation apart from leptons in
// ZZ candidate
int LeptonPlusTauAnalysis::findExtraLeptons(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand) {
  int nExtLep = 0;
  for (const auto& mu: getTightMuList()) {
    TLorentzVector muP4(HZZ4lUtil::getP4(mu));
    if (AnaUtil::sameObject(muP4, Z1Cand.l1P4) || 
	AnaUtil::sameObject(muP4, Z1Cand.l2P4) || 
	AnaUtil::sameObject(muP4, Z2Cand.l1P4) || 
	AnaUtil::sameObject(muP4, Z2Cand.l2P4)) continue;
    if (HZZ4lUtil::pfiso(mu)/mu.pt < 0.4) nExtLep++;
  }
  for (const auto& ele: getTightEleList()) {
    TLorentzVector eleP4(HZZ4lUtil::getP4(ele));
    if (AnaUtil::sameObject(eleP4, Z1Cand.l1P4) || 
	AnaUtil::sameObject(eleP4, Z1Cand.l2P4) || 
	AnaUtil::sameObject(eleP4, Z2Cand.l1P4) || 
	AnaUtil::sameObject(eleP4, Z2Cand.l2P4)) continue;
    if (HZZ4lUtil::pfiso(ele, getEventGridRho())/ele.pt < 0.5) nExtLep++;
  }
  return nExtLep;
}
int LeptonPlusTauAnalysis::findEventCategory(int nleptons, const std::vector<vhtm::Jet>& jetList, int nbjets,
				      const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, bool verbose) {
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
  histf()->cd("LeptonPlusTauAnalysis");
  AnaUtil::fillHist1D("djet", djet, 1);

  int cat = 0;
  if (nleptons == 4 && njets >= 2 && nbjets <= 1 && djet > 0.5) 
    cat = 2;
  else if ( (nleptons == 4 && 
        njets >= 2 && LeptonPlusTauAnalysis::hasJetPair(jetList) &&
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
bool LeptonPlusTauAnalysis::hasJetPair(const std::vector<vhtm::Jet>& jetList) {
  for (unsigned int i = 0; i < jetList.size(); ++i) {
    auto const& j1 = jetList[i];
    TLorentzVector j1P4(HZZ4lUtil::getP4(j1));
    for (unsigned int j = i+1; j < jetList.size(); ++j) {
      auto const& j2 = jetList[j];
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

void LeptonPlusTauAnalysis::endJob() {
  syncDumpf_.close();
  closeFiles();
  
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
  HZZ4lUtil::showEfficiency("muCutFlow", muLabels, "Muon Selection", "Muons");  

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
  HZZ4lUtil::showEfficiency("eleCutFlow", eleLabels, "Electron Selection", "Electrons");  

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
  HZZ4lUtil::showEfficiency("tauCutFlow", tauLabels, "Tau Selection", "Taus");  

  histf()->cd();
  histf()->cd("LeptonPlusTauAnalysis");
  vector<string> evLabels {
    "Events processed",
    "Gen Filter",
    "Events passing trigger",
    "Events with > 0 good vertex",
    "# leptons >= 2",
    "# of Z Candidates > 0",
    "# tau > 0",
    "# leptons >= 3",
    "# extra leptons == 0",
    "has a (Z + l-tau) candidate"
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
  vector<string> evCatList[] {
    "unknwn",
    "mmmtau",
    "mmetau",
    "eeetau",
    "eemtau",
  };
  HZZ4lUtil::showCount("nEventsPerType", evCatList, "Events Per type");  
  if (isMC()) {
    HZZ4lUtil::scaleHistogram("nEventsPerTypeWt", lumiFac);
    HZZ4lUtil::showCount("nEventsPerTypeWt", evCatList, "Events Per type (Weighted)", 3);  
  }

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
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
bool LeptonPlusTauAnalysis::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;
  
  eventFileList_.clear();  

  static const int BUF_SIZE = 256;
  
  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), std::ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }
  
  char buf[BUF_SIZE];
  vector<string> tokens;
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   
    
    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;
    
    // Split the line into words
    AnaUtil::tokenize(line, tokens);
    assert(tokens.size() > 1);
    string key = tokens[0];
    string value = tokens[1];
    if (key == "useEventList")
      useEventList_ = (std::stoi(value.c_str()) > 0) ? true : false;
    else if (key == "skipDuplicate")
      skipDuplicate_ = (std::stoi(value.c_str()) > 0) ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = (std::stoi(value.c_str()) > 0) ? true : false;
    else if (key == "inputEventFile")
      eventFileList_.push_back(value.c_str());
    else if (key == "syncDumpFile")
      dumpFilename_ = value.c_str();
    else if (key == "dumpEventMax")
      dumpEventCount_ = std::stoi(value.c_str());
    
    tokens.clear();
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
  if (skipDuplicate_ && !eventFileList_.empty()) {
    eventMap_.clear();
    for (auto const& f: eventFileList_) {
      cout << ">>> Reading file: " << f << endl;
      ifstream fin(f, std::ios::in);
      if (!fin) {
	cerr << "Input File: " << f << " could not be opened!" << endl;
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
	eventMap_.insert(pair<string, int>(key, 1));

	tokens.clear();
      }
      // Close the file
      fin.close();
    }
    cout << ">>> Total events present: " << eventMap_.size() << endl;
  }  
  selectEvType_ = (static_cast<int>(AnaUtil::cutValue(evselCutMap(), "selectEvType")) > 0) ? true : false;
  printJob();
  
  return true;
}
void LeptonPlusTauAnalysis::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  os << endl;
  os << "   useEventList = " << std::boolalpha << useEventList_ << endl
     << "  skipDuplicate = " << std::boolalpha <<  skipDuplicate_ << endl
     << "    dumpGenInfo = " << std::boolalpha << dumpGenInfo_ << endl;
  if (isMC()) cout << " dumpEventCount = " << dumpEventCount_ << endl;
}
