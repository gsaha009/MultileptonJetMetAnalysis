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

#include "LL4JMETUtil.h"
#include "MultiLeptonJetMetAnalysis.h"
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
MultiLeptonJetMetAnalysis::MultiLeptonJetMetAnalysis()
  : PhysicsObjSelector()
{
}
// ----------
// Destructor
// ----------
MultiLeptonJetMetAnalysis::~MultiLeptonJetMetAnalysis() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MultiLeptonJetMetAnalysis::beginJob() 
{ 
  if (!PhysicsObjSelector::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("MultiLeptonJetMetAnalysis");
  
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
void MultiLeptonJetMetAnalysis::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("MultiLeptonJetMetAnalysis");

  // book histograms to be filled at different stages


  new TH1D("nvtx", "Number of Good vertices", 60, 0, 60);
  if (isMC()) {
    new TH1D("puweight", "PU reweight factor", 100, 0, 2);
    new TH1D("evtweight", "Event weight factor (MC events)", 100, -5, 5);
  }

  //------- Object PLots -----------------------------------------------
  new TH1D("nGoodMuon", "Number of Good muons", 6, -0.5, 5.5);
  new TH1D("nGoodEle", "Number of Good electrons", 6, -0.5, 5.5);
  new TH1D("nGoodTau", "Number of Good Tau", 5, -0.5, 4.5);
  new TH1D("nTJets", "Number of Tight Jets", 10, -0.5, 9.5);
  new TH1D("nLJets", "Number of Loose Jets", 10, -0.5, 9.5);
  new TH1D("jetpT", "pT of all Tight Jets", 600, 0.0, 300.);
  new TH1D("jetpT1", "pT of 1st Tight Jets", 600, 0.0, 300.);
  new TH1D("jetpT2", "pT of 2nd Tight Jets", 600, 0.0, 300.);
  new TH1D("jetpT3", "pT of 3rd Tight Jets", 600, 0.0, 300.);
  new TH1D("jetpT4", "pT of 4th Tight Jets", 600, 0.0, 300.);
  new TH1D("jetpT5", "pT of 5th Tight Jets", 600, 0.0, 300.);
  new TH1D("nbJets", "Number of b-Jets", 10, -0.5, 9.5);
  new TH1D("nLeptonCand", "No of total mu e candidates", 10, -0.5, 9.5);
  new TH1D("nFinalLeptons", "Final leptons after pT cut and dR matching", 10, -0.5, 9.5);
  new TH1D("evtCutFlow", "Event CutFlow", 13, -0.5, 12.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 13, -0.5, 12.5);
  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);
  new TH1D("nLooseLeptons", "Total tight isolated leptons", 10, -0.5, 9.5);
  new TH1D("nTightIsoLeptons", "Total tight isolated leptons", 8, -0.5, 7.5);
  new TH1D("LeptonsPt", "pT of the tight isolated lepton candidates", 200, 0., 400.);
  new TH1D("met", "Missing Transver Energy", 150, 0, 150);
  new TH1D("corrmet", "Missing Transver Energy (Corrected)", 150, 0, 150);
  new TH1D("puppimet", "Missing Transver Energy (PUPPI)", 150, 0, 150);
  //  new TH1F("dRll",  "dR(lep, lep) of Z", 100, 0, 2 * M_PI);
  
  //  new TH2F("h1VsZpT", "h1 vs Z pT scatter plot", 100, 0, 400, 100, 0, 400);
  //histf()->cd("InFinalState");

  new TH1D("FinalLepPt", "pT of the tight isolated final lepton candidates", 200, 0., 400.);  
  new TH1D("#TightJets", "no of finally selected jets", 10, -0.5, 9.5);
  new TH1D("TightJetsPt", "Pt of finally selected jets", 200, 0., 400.);
  new TH1F("difference", "difference bet #LooseIsoEle & TightIsoEle", 9, -4.5, 4.5);  
  new TH1F("MatchedJet_pt", "", 200, 0., 400.);
  new TH1F("genq_pt", "", 200, 0., 400.);

  histf()->cd();
  histf()->ls();
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void MultiLeptonJetMetAnalysis::clearLists() {
  PhysicsObjSelector::clear();
  vtxList_.clear();
  LepCandList_.clear();
  finalLepCandList_.clear();
}
// -------------------
// The main event loop
// -------------------
void MultiLeptonJetMetAnalysis::eventLoop()
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
    histf()->cd("MultiLeptonJetMetAnalysis");
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
      genAna_ -> analyze(histf());            
      bool genOk = genAna_->filter();
      if (!genOk) continue;
    }

    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, puevWt);

    histf()->cd();
    histf()->cd("MultiLeptonJetMetAnalysis");
    
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
    
    //*********** M A I N  A N A L Y S I S  O B J E C T  S E L E C T I O N ************//
    double vz = (vtxList_.size() > 0) ? vtxList_.at(0).z : -999;
    findObjects(vz, 1); // do not use event weight
    





    //if (genAna_ != nullptr && dumpGenInfo_) genAna_->dumpEvent();
    //if (getFSRPhotonList().size()) dumpEvent(vz, false, true);

    // Access selected objects 
    const auto& elePhotonPairList = getTightIsoElePhotonPairList();
    int nEle = elePhotonPairList.size();

    //const auto& TightIsoEleList = getTightIsoEleList();
    //int nTightIsoEle = TightIsoEleList.size();

    //const auto& TightIsoMuList = getTightIsoMuList();
    //int nTightIsoMu = TightIsoMuList.size();

    const auto& muPhotonPairList  = getTightIsoMuPhotonPairList();
    int nMuon = muPhotonPairList.size();

    const auto& tauList = getIsoTauList();
    //int nTau = tauList.size();
    
    const auto& TightJets = getTightJetList();
    int nTightJets = TightJets.size();
    //    std::cout<<"#tight jets: "<<nTightJets<<std::endl;
    //const auto& LooseEleList = getLooseEleList();
    //const auto& LooseMuList = getLooseMuList();

    const auto& JetList = getLooseJetList();
    //    std::cout<<"#loose jets: "<<JetList.size()<<std::endl;
    //********************** A N A L Y S I S **********************//
    
    histf()->cd();
    histf()->cd("MultiLeptonJetMetAnalysis");
    AnaUtil::fillHist1D("nGoodMuon", nMuon, puevWt);
    AnaUtil::fillHist1D("nGoodEle", nEle, puevWt);
    AnaUtil::fillHist1D("met", metColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("corrmet", corrmetColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("puppimet", puppimetColl()->at(0).met, puevWt);

    //P A C K I N G  L E P T O N S
    if (nMuon > 0) packLeptons<vhtm::Muon>(muPhotonPairList, LepCandList_); 
    if (nEle > 0)  packLeptons<vhtm::Electron>(elePhotonPairList, LepCandList_);
    AnaUtil::fillHist1D("nLeptonCand", LepCandList_.size(), puevWt);
    
    // Printing dR value between TightIsoLeptons and TightJets
    /*
    for (auto& l: LepCandList_){
      TLorentzVector LepP4 = l.lP4;
      for (auto& j: TightJets){
	TLorentzVector JetP4(LL4JMETUtil::getP4(j));
	//	if (LepP4.DeltaR(JetP4) <= 0.4)
	std::cout<<"dR value between TightIsoLeptons and TightJets"<<"\t"<<LepP4.DeltaR(JetP4)<<std::endl;
      }
    }
    */
    
    addLeptonIsolation(LepCandList_, elePhotonPairList, muPhotonPairList); //add lepton isolation to lepcandidates    

    //N O  b-J E T S
    if (nbJets() > 0) continue;  
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, puevWt);

    AnaUtil::fillHist1D("nTJets", nTightJets, puevWt);
    AnaUtil::fillHist1D("nLJets", JetList.size(), puevWt);
    int i=0;
    for (auto& j: TightJets){
      AnaUtil::fillHist1D("jetpT", j.pt, puevWt);
      i++;
      if (i == 1) AnaUtil::fillHist1D("jetpT1", j.pt, puevWt);
      if (i == 2) AnaUtil::fillHist1D("jetpT2", j.pt, puevWt);
      if (i == 3) AnaUtil::fillHist1D("jetpT3", j.pt, puevWt);
      if (i == 4) AnaUtil::fillHist1D("jetpT4", j.pt, puevWt);
      if (i == 5) AnaUtil::fillHist1D("jetpT5", j.pt, puevWt);
    }
    
    //O N L Y  4  J E T S  W I T H  PT > 30 GEV (S E E  J O B  C A R D)
    if (nTightJets < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, puevWt);

    if (nTightJets < 3) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, puevWt);

    if (nTightJets < 4) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, puevWt);
    //    std::cout<<"Evt No: "<<ev<<"\t"<<"#jets >=4 after bjet veto: "<<"\t"<<TightJets.size()<<std::endl;

    //N O  O F  L E P T O N S = 2 (PT > 10 GEV)
    if (LepCandList_.size() < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, puevWt);
    //    std::cout<<"Evt No: "<<ev<<"\t"<<"#Leptons befr charge and dR matching: "<<"\t"<<LepCandList_.size()<<std::endl;

    //L E A D I N G  J E T  PT >= 60 GEV
    if (TightJets.at(0).pt < 90) continue;
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, puevWt);

    //S U B  L E A D I N G  J E T  PT >= 40 GEV
    if (TightJets.at(1).pt < 60) continue;
    AnaUtil::fillHist1D("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 10, puevWt);
    //    std::cout<<"Evt No: "<<ev<<"\t"<<"#jets >=4 after bjet veto and pT cuts on 2 leading jets: "<<"\t"<<TightJets.size()<<std::endl;

    //M E T > 30 GEV
    if (puppimetColl()->at(0).met < 30) continue;
    AnaUtil::fillHist1D("evtCutFlow", 11);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 11, puevWt);

    //S E A R C H I N G  F O R  T W O  T I G H T  I S O L A T E D  L E P T O N S  O F  S A M E  C H A R G E
      //   std::vector <LeptonCandidate> finalLepCandList_;
    for (unsigned int i = 0; i < LepCandList_.size(); ++i){
      const auto& lepi = LepCandList_.at(i);
      double lepicharge = lepi.lCharge;
      TLorentzVector lepiP4 = lepi.lP4;
      for (unsigned int j = i+1; j < LepCandList_.size(); ++j){
        const auto& lepj = LepCandList_.at(j);
        double lepjcharge = lepj.lCharge;
        TLorentzVector lepjP4 = lepj.lP4;
        if ((lepicharge*lepjcharge) > 0 && (lepiP4.DeltaR(lepjP4) > 0.4  && lepiP4.DeltaR(lepjP4) < 1.5)) {
          finalLepCandList_.push_back(lepi);
          finalLepCandList_.push_back(lepj); //??????????????????????????????
          // lepCandPairList_.push_back({lepi, lepj});
	}
      }
    }
    AnaUtil::fillHist1D("nFinalLeptons", finalLepCandList_.size(), puevWt);
    //new_addition
    //    if (getLooseMuList().size() + getLooseEleList().size() > 2) continue;   // events should not have any extra loose muons/electrons
    if (!tauList.empty()) continue;
    //new_addition
    if (finalLepCandList_.size() != 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 12);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12, puevWt);

    evLog() << run << " " << lumis << " " << event << std::endl;

    // Print only the first n events; n configurable
    if (isMC() && dumpEventCount_ > -1 && ++nEventSel >= dumpEventCount_) continue;

    cout << ">>> "
         << "<nLooseMuon>: " << getLooseMuList().size()
         << ", <nMuon>: " << muPhotonPairList.size()
         << ", <nLooseElectron>: " << getLooseEleList().size()
         << ", <nElectron>: " << elePhotonPairList.size()
         << ", <nTau>: " << tauList.size()
         << ", <nTightJets> " << TightJets.size()
	 << endl;

    //dumpEvent(vz, false, true);
  }
  // Analysis over
  endJob();
  //    if (corrmetColl()->at(0).met > //AnaUtil::cutValue(evselCutMap(), "maxMET")) continue;
}

void MultiLeptonJetMetAnalysis::endJob() {
  PhysicsObjSelector::endJob();

  histf()->cd();
  histf()->cd("MultiLeptonJetMetAnalysis");
  vector<string> evLabels {
    "Events processed",
    "Gen Filter",
    "Events with > 0 good vertex",
    "Events passing trigger",
    "Events have no b-Jet (C1)",
    "Events have >= 2 jets of pT > 30 Gev",
    "Events have >= 3 jets of pT > 30 Gev",
    "Events have >= 4 jets of pT > 30 Gev (C2)",
    "Events have >= 2 IsoLeptons with pT > 10 Gev (C3)",
    "Leading Jet pT >= 60 Gev (C4)",
    "Sub leading jet pT >= 40 Gev (C5)",
    "Events have MET > 30 Gev (C6)",
    "Events have same charged di leptons within 0.4<dR<1.5 (C7)"
  };
  LL4JMETUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection");  

  double lumiFac = 1.0;
  if (isMC()) {
    lumiFac = lumiWt(evtWeightSum_);
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
  }

  if (isMC()) {
    LL4JMETUtil::scaleHistogram("evtCutFlowWt", lumiFac);
    LL4JMETUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");  
  }
  
  // Event category
  vector<string> evCatList {
    "unknwn",
    "mmem",
    "eeem"
  };
  LL4JMETUtil::showCount("nEventsPerType", evCatList, "Events Per type");  
  if (isMC()) {
    LL4JMETUtil::scaleHistogram("nEventsPerTypeWt", lumiFac);
    LL4JMETUtil::showCount("nEventsPerTypeWt", evCatList, "Events Per type (Weighted)", 3);  
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
bool MultiLeptonJetMetAnalysis::readJob(const string& jobFile, int& nFiles)
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
    std::cout << line << std::endl;
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
	//AnaUtil::tokenize(line, tokens);
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
void MultiLeptonJetMetAnalysis::printJob(ostream& os) const
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

  //AnaUtil::showList(eventFilelist_, ">>> INFO. Input event files:", os);
}
