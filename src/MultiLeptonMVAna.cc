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
#include "TCanvas.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"

#include "configana.h"
#include "LL4JMETUtil.h"
#include "MultiLeptonMVAna.h"


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
MultiLeptonMVAna::MultiLeptonMVAna()
  : PhysicsObjSelector()
   {}

// ----------
// Destructor
// ----------
MultiLeptonMVAna::~MultiLeptonMVAna() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MultiLeptonMVAna::beginJob() 
{ 
  if (!PhysicsObjSelector::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("TMVAnalysis");
  
  bookHistograms();

  if (_createMVATree) {
    skimObj_ = std::make_unique <MVASkim> (_mvaInputFile);
    if (!skimObj_) return false;
  }

  else if (_readMVA) {
    _mvaObj = std::make_unique<MVAnalysis>(_MVAnetwork, _MVAxmlFile);
    if (!_mvaObj) return false;
  }

#ifdef  SKIP_DUPLICATE_ALL
  eventIdStore_.clear();
#endif

  return true;
}

// ---------------
// Book histograms
// ---------------
void MultiLeptonMVAna::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("TMVAnalysis");

  // book histograms to be filled at different stages


  new TH1D("nvtx_0", "Number of Good vertices", 60, 0, 60);
  if (isMC()) {
    new TH1D("puweight", "PU reweight factor", 1000, -20., 20.);
    new TH1D("pu_evtWt", "PU nd Evt Wt", 1000, -20., 20.);
    new TH1D("evtweight", "Event weight factor (MC events)", 20, -10., 10.);
  }
  new TH1D("evtCutFlow", "Event CutFlow", 25, -0.5, 24.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 25, -0.5, 24.5);
  new TH1D("met", "Missing Transver Energy",200, 0, 200);
  new TH1D("corrmet", "Missing Transver Energy (CORR)", 600, 0., 400.);
  new TH1D("puppimet", "Missing Transver Energy (PUPPI)", 200, 0, 200);


  //stage_0

  new TH1F("nLepCand_0" ,"", 10, -0.5, 9.5);
  new TH1F("nTightJet_0" ,"", 10, -0.5, 9.5);
  new TH1F("nLooseJet_0" ,"", 10, -0.5, 9.5);
  new TH1D("corrmet_0", "Missing Transver Energy (CORR)", 600, 0., 400.);
  new TH1F("nbJets_0" ,"", 10, -0.5, 9.5);
  
  
  //stage_1: nLep = 2

  new TH1D("lep1pt_1" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_1" ,"", 300, 0., 300.0);
  new TH1D("lepDR_1" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_1" ,"", 500, 0.0, 5.0);
  new TH1D("lSumPt_1" ,"", 500, 0., 500.0);
  new TH1D("lT_1" ,"", 1000, 0., 1000.0);
  new TH1D("lep1MetDPhi_1", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_1", "", 500, 0., 5.);
  new TH1D("l1TrMass_1", "Trans mass of lep1 & met", 500, 0., 600.);
  new TH1D("l2TrMass_1", "Trans mass of lep2 & met", 500, 0., 600.);
  new TH1D("nTightJet_1" ,"", 10, -0.5, 9.5);
  new TH1D("corrmet_1", "Missing Transver Energy", 600, 0., 400.);
  new TH1D("l1l2InvM_1", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("l1l2SameChrInvM_1", "InvMass of 2SameChargedLeptons", 400, 0., 400.);
  new TH1D("nbJets_1" ,"", 10, -0.5, 9.5);
  new TH1D("HppTrMass_1", "", 600, 0., 800.);

  //stage_2: nTightJets >= 2
  new TH1D("nTightJet_2" ,"", 10, -0.5, 9.5);
  new TH1D("jet1pt_2" ,"", 300, 0., 300.0);
  new TH1D("jet2pt_2" ,"", 300, 0., 300.0);
  new TH1D("jdR_2", "",100,0.0, 5.0);
  new TH1D("jdPhi_2", "",500,0.0, 5.0);
  new TH1D("hT_2" ,"", 1000, 0., 1000.0);
  new TH1D("sT_2" ,"", 1000, 0., 1000.0);
  new TH1D("JetLepDPhi_2", "",500,0.0, 5.0);
  new TH1D("corrmet_2", "Missing Transver Energy", 600, 0., 400.);
  new TH1D("l1l2InvM_2", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_2", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_2", "", 500, 0., 5.);
  new TH1D("nbJets_2" ,"", 10, -0.5, 9.5);
  new TH1D("MR_2" ,"", 600, 0.0, 800.0);
  new TH1D("HppTrMass_2", "", 600, 0., 800.);

  //stage_3: jet1Pt > 60
  new TH1D("nTightJet_3" ,"", 10, -0.5, 9.5);
  new TH1D("hT_3" ,"", 1000, 0., 1000.0);
  new TH1D("sT_3" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_3" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_3" ,"", 300, 0., 300.0);
  new TH1D("lepDR_3" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_3" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_3", "",100,0.0, 5.0);
  new TH1D("jdPhi_3", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_3", "",500,0.0, 5.0);
  new TH1D("lSumPt_3" ,"", 500, 0., 500.0);
  new TH1D("lT_3" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_3", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_3", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_3", "", 500, 0., 5.);
  new TH1D("nbJets_3" ,"", 10, -0.5, 9.5);



  //stage_4: jet2Pt > 40

  new TH1D("nTightJet_4" ,"", 10, -0.5, 9.5);
  new TH1D("hT_4" ,"", 1000, 0., 1000.0);
  new TH1D("sT_4" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_4" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_4" ,"", 300, 0., 300.0);
  new TH1D("lepDR_4" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_4" ,"", 500, 0.0, 5.0);
  new TH1D("jdPhi_4", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_4", "",500,0.0, 5.0);
  new TH1D("lSumPt_4" ,"", 500, 0., 500.0);
  new TH1D("lT_4" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_4", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_4", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_4", "", 500, 0., 5.);
  new TH1D("nbJets_4" ,"", 10, -0.5, 9.5);
  new TH1D("lepSumChr_4", "", 5, -2.5, 2.5);


  //stage_5: met > 30

  new TH1D("nTightJet_5" ,"", 10, -0.5, 9.5);
  new TH1D("hT_5" ,"", 1000, 0., 1000.0);
  new TH1D("sT_5" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_5" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_5" ,"", 300, 0., 300.0);
  new TH1D("lepDR_5" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_5" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_5", "",100,0.0, 5.0);
  new TH1D("jdPhi_5", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_5", "",500,0.0, 5.0);
  new TH1D("lSumPt_5" ,"", 500, 0., 500.0);
  new TH1D("lT_5" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_5", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_5", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_5", "", 500, 0., 5.);
  new TH1D("nbJets_5" ,"", 10, -0.5, 9.5);
  new TH1D("HppTrMass_5", "", 600, 0., 800.);

  //stage_6 : j1j2dr < 1.6
  //  new TH1D("nbJets_13" ,"", 10, -0.5, 9.5);
  new TH1D("nTightJet_6" ,"", 10, -0.5, 9.5);
  new TH1D("hT_6" ,"", 1000, 0., 1000.0);
  new TH1D("sT_6" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_6" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_6" ,"", 300, 0., 300.0);
  //  new TH1D("lepDR_11" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_6" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_6", "",100,0.0, 5.0);
  new TH1D("jdPhi_6", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_6", "",500,0.0, 5.0);
  new TH1D("lSumPt_6" ,"", 500, 0., 500.0);
  new TH1D("lT_6" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_6", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_6", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_6", "", 500, 0., 5.);




  //stage_7
  //  new TH1D("nbJets_13" ,"", 10, -0.5, 9.5);
  new TH1D("nTightJet_7" ,"", 10, -0.5, 9.5);
  new TH1D("hT_7" ,"", 1000, 0., 1000.0);
  new TH1D("sT_7" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_7" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_7" ,"", 300, 0., 300.0);
  //  new TH1D("lepDR_11" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_7" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_7", "",100,0.0, 5.0);
  new TH1D("jdPhi_7", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_7", "",500,0.0, 5.0);
  new TH1D("lSumPt_7" ,"", 500, 0., 500.0);
  new TH1D("lT_7" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_7", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_7", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_7", "", 500, 0., 5.);
  new TH1D("HppTrMass_7", "", 600, 0., 800.);




  //stage_8
  //  new TH1D("nbJets_13" ,"", 10, -0.5, 9.5);
  new TH1D("nTightJet_8" ,"", 10, -0.5, 9.5);
  new TH1D("hT_8" ,"", 1000, 0., 1000.0);
  new TH1D("sT_8" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_8" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_8" ,"", 300, 0., 300.0);
  //  new TH1D("lepDR_11" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_8" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_8", "",100,0.0, 5.0);
  new TH1D("jdPhi_8", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_8", "",500,0.0, 5.0);
  new TH1D("lSumPt_8" ,"", 500, 0., 500.0);
  new TH1D("lT_8" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_8", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_8", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_8", "", 500, 0., 5.);
  new TH1D("HppTrMass_8", "", 600, 0., 800.);


  //stage_9
  //Pre Selection Done
  new TH1D("nTightJet_PS" ,"", 10, -0.5, 9.5);
  new TH1D("hT_PS" ,"", 1000, 0., 1000.0);
  new TH1D("sT_PS" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_PS" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_PS" ,"", 300, 0., 300.0);
  //  new TH1D("lepDR_11" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_PS" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_PS", "",100,0.0, 5.0);
  new TH1D("jdPhi_PS", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_PS", "",500,0.0, 5.0);
  new TH1D("lSumPt_PS" ,"", 500, 0., 500.0);
  new TH1D("lT_PS" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_PS", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_PS", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_PS", "", 500, 0., 5.);



  //MVA1

  new TH1D("nvtx_MVA1", "Number of Good vertices", 60, 0, 60);
  new TH1D("nTightJet_MVA1" ,"", 10, -0.5, 9.5);
  new TH1D("hT_MVA1" ,"", 1000, 0., 1000.0);
  new TH1D("sT_MVA1" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_MVA1" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_MVA1" ,"", 300, 0., 300.0);
  new TH1D("lep1eta_MVA1" ,"", 100, -3.0, 3.0);
  new TH1D("lep2eta_MVA1" ,"", 100, -3.0, 3.0);

  new TH1D("jet1pt_MVA1" ,"", 400, 0., 600.0);
  new TH1D("jet2pt_MVA1" ,"", 400, 0., 600.0);
  new TH1D("jet1eta_MVA1" ,"", 100, -5.0, 5.0);
  new TH1D("jet2eta_MVA1" ,"", 100, -5.0, 5.0);

  new TH1D("lepDR_MVA1" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_MVA1" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_MVA1", "",100,0.0, 5.0);
  new TH1D("jdPhi_MVA1", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_MVA1", "",500,0.0, 5.0);
  new TH1D("lSumPt_MVA1" ,"", 500, 0., 500.0);
  new TH1D("lT_MVA1" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_MVA1", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_MVA1", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_MVA1", "", 500, 0., 5.);
  new TH1D("HppTrMass_MVA1", "", 600, 0., 800.);



  new TH1D("nvtx_MVA2", "Number of Good vertices", 60, 0, 60);
  new TH1D("nTightJet_MVA2" ,"", 10, -0.5, 9.5);
  new TH1D("hT_MVA2" ,"", 1000, 0., 1000.0);
  new TH1D("sT_MVA2" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_MVA2" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_MVA2" ,"", 300, 0., 300.0);
  new TH1D("lep1eta_MVA2" ,"", 100, -3.0, 3.0);
  new TH1D("lep2eta_MVA2" ,"", 100, -3.0, 3.0);

  new TH1D("jet1pt_MVA2" ,"", 400, 0., 600.0);
  new TH1D("jet2pt_MVA2" ,"", 400, 0., 600.0);
  new TH1D("jet1eta_MVA2" ,"", 100, -5.0, 5.0);
  new TH1D("jet2eta_MVA2" ,"", 100, -5.0, 5.0);

  new TH1D("lepDR_MVA2" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_MVA2" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_MVA2", "",100,0.0, 5.0);
  new TH1D("jdPhi_MVA2", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_MVA2", "",500,0.0, 5.0);
  new TH1D("lSumPt_MVA2" ,"", 500, 0., 500.0);
  new TH1D("lT_MVA2" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_MVA2", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_MVA2", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_MVA2", "", 500, 0., 5.);
  new TH1D("HppTrMass_MVA2", "", 600, 0., 800.);

  if (isMC())  new TH1D("diffFlvYield", "", 3, -0.5, 2.5);

  histf()->cd();
  histf()->ls();
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------

void MultiLeptonMVAna::clearLists() {
  PhysicsObjSelector::clear();
  vtxList_.clear();
  LepCandList_.clear();
  JetPairSelected.clear();
  JetSelected.clear();
  SignalJetsSet.clear();
  SignalJets.clear();
}

// -------------------
// The main event loop
// -------------------

void MultiLeptonMVAna::eventLoop()
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
  //  TH1F *wtDiff = new TH1F ("wtDiff", "", 20, 0., 20.);  
 
  for (int ev = fevt; ev < levt; ++ev) {
    clearEvent(); // reset tree variables 
    clearLists(); // reset analysis related lists for each event
    int lflag = chain()->LoadTree(ev);
    int nbytes = getEntry(lflag);    // returns total bytes read
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));

    // For data or for MC without pileup
    double puevWt = 1; //for Data
#if 0
    if (isMC() && usePUWt()) {
      int npu = 0;
      puevWt = wtPileUp(npu);
    }
#endif
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
    histf()->cd("TMVAnalysis");
    
    double wt = 1;
    if (isMC()) {
      AnaUtil::fillHist1D("puweight", puevWt);
      // only when exclusive jets datasets are used (for np = 0)
      if (0) std::cout << "== nMEPartons: " << genEventColl()->at(0).nMEPartons 
		       << ", lheNOutPartons: " << genEventColl()->at(0).lheNOutPartons 
		       << std::endl;
      if (selectPM_ && genEventColl()->at(0).nMEPartons != nMEPartons_) continue;

      //double wt = 1.0;
      //if (!isSignal()) wt = genEventColl()->size() ? genEventColl()->at(0).evtWeight : 1;
      wt = genEventColl()->size() ? genEventColl()->at(0).evtWeight : 1;
      //      AnaUtil::fillHist1D("evtweight", wt);      
      if (wt > 0.) wt = 1.0;
      else if (wt < 0.) wt = -1.0;
      else if (wt == 0.) wt = 0.0;
      AnaUtil::fillHist1D("evtweight", wt);      
      if (usePUWt()) {
	int npu = 0;
	puevWt  = wtPileUp(npu, false);
	AnaUtil::fillHist1D("puweight", puevWt);
	puevWt *= wt;
	AnaUtil::fillHist1D("pu_evtWt", puevWt);
      }
      evtWeightSum_ += wt;   // this is used for final normalization
    }

    /*  
#ifdef    SKIP_DUPLICATE_ALL
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

    if (!isMC()) {
#ifdef SKIP_DUPLICATE_IN_MEMORY
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      std::string evs {mkey.str()};
      if (eventIdStore_.find(evs) != eventIdStore_.end()) {
        if (0) cout << "DuplicateAll: " << evs << endl;
	continue;
      }
      else {
	eventIdStore_.insert({evs, 1});
      }
#endif
#ifdef SKIP_DUPLICATE_FROM_FILE
      if (skipDuplicate_) {
	std::ostringstream mkey;
        mkey << run << "-" << lumis << "-" << event;
	std::string evs {mkey.str()};
        if (eventIdStore_.find(evs) != eventIdStore_.end()) {
          if (0) cout << "Duplicate Event: " << evs << endl;
          continue;
        }
      }
      evLog() << run << " " << lumis << " " << event << std::endl;
#endif
    }
    */

    histf()->cd(); //required
    histf()->cd("TMVAnalysis");
    AnaUtil::fillHist1D("evtCutFlow", 0);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, puevWt);
    
    // at least 1 good PV
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    size_t ngoodVtx = vtxList_.size();
    AnaUtil::fillHist1D("nvtx_0", ngoodVtx, puevWt);
    if (ngoodVtx < 1) continue;

    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, puevWt);


    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt); 
    
    // Is event triggered?
    if (0) dumpTriggerPaths(std::cout, true);
    if (useTrigger() && !isTriggered(true, false)) continue;
    
    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, puevWt);
    
    //*********** M A I N  A N A L Y S I S  O B J E C T  S E L E C T I O N ************//
    double vz = (vtxList_.size() > 0) ? vtxList_.at(0).z : -999;
    findObjects(vz, 1); // do not use event weight
    

    //if (genAna_ != nullptr && dumpGenInfo_) genAna_->dumpEvent();
    //if (getFSRPhotonList().size()) dumpEvent(vz, false, true);

    // Access selected objects 
    const auto& elePhotonPairList = getTightIsoElePhotonPairList();
    int nEle = elePhotonPairList.size();

    const auto& muPhotonPairList  = getTightIsoMuPhotonPairList();
    int nMuon = muPhotonPairList.size();

    const auto& tauList = getIsoTauList();
    int nTau = tauList.size();
    
    const auto& TightJets = getTightJetList();
    //    int nTightJets = TightJets.size();
    
    const auto& LooseJets = getLooseJetList();
    
    //********************** A N A L Y S I S **********************//
    
    histf()->cd();//Required
    histf()->cd("TMVAnalysis");

    AnaUtil::fillHist1D("met", metColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("corrmet", corrmetColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("puppimet", puppimetColl()->at(0).met, puevWt);

    //P A C K I N G  L E P T O N S in LepCandList_
    if (nMuon > 0) packLeptons<vhtm::Muon>(muPhotonPairList, LepCandList_); 
    if (nEle > 0)  packLeptons<vhtm::Electron>(elePhotonPairList, LepCandList_);
    AnaUtil::fillHist1D("nLeptonCand", LepCandList_.size(), puevWt);

    addLeptonIsolation(LepCandList_, elePhotonPairList, muPhotonPairList); //add lepton isolation to lepcandidates    
    std::sort(std::begin(LepCandList_), std::end(LepCandList_), PtComparatorLep<LeptonCandidate>()); //sorting lepton candidates  
    histf()->cd();
    histf()->cd("TMVAnalysis");


    //_________________________________________PreSelection Begins_______________________________________________//

    //0...) Start
    const MET& mt = corrmetColl()->at(0);
    double nLeps = LepCandList_.size();
    double nJets = TightJets.size();
    double nLJets = LooseJets.size();
    AnaUtil::fillHist1D("nLepCand_0",  nLeps, puevWt);
    AnaUtil::fillHist1D("nTightJet_0", nJets, puevWt);
    AnaUtil::fillHist1D("nLooseJet_0", nLJets, puevWt);
    AnaUtil::fillHist1D("corrmet_0",  mt.met, puevWt);
    AnaUtil::fillHist1D("nbJets_0", nbJets(), puevWt);

 
    //1...) Cut (nLep = 2)
    if (LepCandList_.size() != 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, puevWt);

    TLorentzVector l1p4 = LepCandList_[0].lRP4;
    TLorentzVector l2p4 = LepCandList_[1].lRP4;

    double Mass_2l = (l1p4 + l2p4).M();
    AnaUtil::fillHist1D("l1l2InvM_1", Mass_2l, puevWt);

    double lep1c = LepCandList_[0].lCharge;
    double lep2c = LepCandList_[1].lCharge;
    if (lep1c*lep2c > 0){
      AnaUtil::fillHist1D("l1l2SameChrInvM_1", Mass_2l, puevWt);
    }

    double lep1pt = l1p4.Pt();
    double lep2pt = l2p4.Pt();
    AnaUtil::fillHist1D("lep1pt_1", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_1", lep2pt, puevWt);

    double lep1eta = l1p4.Eta();
    double lep2eta = l2p4.Eta();

    double lep1phi = l1p4.Phi();
    double lep2phi = l2p4.Phi();
 
    double dR = l1p4.DeltaR(l2p4);
    double lepDPhi = std::fabs(TVector2::Phi_mpi_pi(l1p4.Phi()-mt.metphi));
    AnaUtil::fillHist1D("lepDR_1", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_1", lepDPhi, puevWt);

    double lSumPt = lep1pt + lep2pt;
    AnaUtil::fillHist1D("lSumPt_1", lSumPt, puevWt);
    double lT = lSumPt + mt.met;
    AnaUtil::fillHist1D("lT_1", lT, puevWt);

    double mtPhi = mt.metphi;

    double l1metDPhi = std::fabs(TVector2::Phi_mpi_pi(l1p4.Phi()-mt.metphi));
    double l2metDPhi = std::fabs(TVector2::Phi_mpi_pi(l2p4.Phi()-mt.metphi));
    AnaUtil::fillHist1D("lep1MetDPhi_1", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_1", l2metDPhi, puevWt);

    //double transMass = (l1p4+l2p4).Mt();
    double mt_l1met = TMath::Sqrt(2*lep1pt*(mt.met)*(1-TMath::Cos(l1metDPhi)));
    double mt_l2met = TMath::Sqrt(2*lep2pt*(mt.met)*(1-TMath::Cos(l2metDPhi)));

    AnaUtil::fillHist1D("l1TrMass_1", mt_l1met, puevWt);
    AnaUtil::fillHist1D("l2TrMass_1", mt_l2met, puevWt);

    AnaUtil::fillHist1D("nTightJet_1", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("corrmet_1",  mt.met, puevWt);
    AnaUtil::fillHist1D("nbJets_1", nbJets(), puevWt);

    double HppTrMass = TMath::Sqrt(2*((l1p4+l2p4).Pt())*(mt.met)*(1-TMath::Cos(TVector2::Phi_mpi_pi((l1p4+l2p4).Phi()-mt.metphi))));
    AnaUtil::fillHist1D ("HppTrMass_1", HppTrMass, puevWt); 

    //2...) Cut (1 + nTightJet >= 2)
    if (TightJets.size() < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, puevWt);
    AnaUtil::fillHist1D("nTightJet_2", TightJets.size(), puevWt);


    AnaUtil::fillHist1D("jet1pt_2", TightJets[0].pt, puevWt);
    AnaUtil::fillHist1D("jet2pt_2", TightJets[1].pt, puevWt);

    TLorentzVector j1P4 = LL4JMETUtil::getP4(TightJets[0]);
    TLorentzVector j2P4 = LL4JMETUtil::getP4(TightJets[1]);
    double j1j2DR = j1P4.DeltaR(j2P4);
    double j1j2DPhi = std::fabs(TVector2::Phi_mpi_pi(j1P4.Phi()-j2P4.Phi()));
    AnaUtil::fillHist1D("jdR_2", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_2", j1j2DPhi, puevWt);

    double hT = 0.0;
    TLorentzVector jetp4;
    jetp4.SetPtEtaPhiE(0.,0.,0.,0.);
    for (auto &j: TightJets){
      TLorentzVector jp4 = LL4JMETUtil::getP4(j);
      hT += jp4.Pt();
      jetp4 += jp4;
    }

    double hTvec = jetp4.Pt();
    AnaUtil::fillHist1D("hT_2", hT, puevWt);
    double j1l1dr = j1P4.DeltaR(l1p4);
    double j1l2dr = j1P4.DeltaR(l2p4);
    double j2l1dr = j2P4.DeltaR(l1p4);
    double j2l2dr = j2P4.DeltaR(l2p4);

    double j1metDPhi = std::fabs(TVector2::Phi_mpi_pi(j1P4.Phi()-mt.metphi));
    double j2metDPhi = std::fabs(TVector2::Phi_mpi_pi(j2P4.Phi()-mt.metphi));

    double sT = hT + lep1pt + lep2pt + mt.met;
    AnaUtil::fillHist1D("sT_2", sT, puevWt);
    double metSqrtST = mt.met/TMath::Sqrt(sT);

    TLorentzVector j1j2p4 = j1P4 + j2P4;
    TLorentzVector l1l2p4 = l1p4 + l2p4;
    double jlDPhi = std::fabs(TVector2::Phi_mpi_pi(j1j2p4.Phi()-l1l2p4.Phi()));
    AnaUtil::fillHist1D("JetLepDPhi_2", jlDPhi, puevWt);
    AnaUtil::fillHist1D("corrmet_2",  mt.met, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_2", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_2", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_2", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_2", nbJets(), puevWt);

    double MR = TMath::Sqrt((j1P4.Pt()*j1P4.Pt()) + (j2P4.Pt()*j2P4.Pt()) + 2*j1P4.P()*j2P4.P()*TMath::Cos(std::abs(j1P4.Eta()-j2P4.Eta())) - 2*j1P4.Pz()*j2P4.Pz());
    AnaUtil::fillHist1D("MR_2", MR,  puevWt);
    AnaUtil::fillHist1D ("HppTrMass_2", HppTrMass, puevWt); 

    //3...) Cut (2 + jet1Pt > 60)

    //    if (j1P4.Pt() < 60) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, puevWt);

    AnaUtil::fillHist1D("nTightJet_3", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_3", hT, puevWt);
    AnaUtil::fillHist1D("sT_3", sT, puevWt);
    AnaUtil::fillHist1D("jet2pt_3", TightJets[1].pt, puevWt);
    AnaUtil::fillHist1D("corrmet_3",  mt.met, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_3", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_3", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_3", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_3", nbJets(), puevWt);


    //4...) Cut (3 + jet2Pt > 40)
    //  if (j2P4.Pt() < 40) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, puevWt);
    AnaUtil::fillHist1D("nTightJet_4", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_4", hT, puevWt);
    AnaUtil::fillHist1D("sT_4", sT, puevWt);
    AnaUtil::fillHist1D("jdR_4", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_4", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_4", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_4", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_4", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_4", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_4", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_4", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_4", lT, puevWt);
    AnaUtil::fillHist1D("corrmet_4",  mt.met, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_4", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_4", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_4", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_4", nbJets(), puevWt);



    //5...) Cut (4 + met > 30)
    if (mt.met < 30) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, puevWt);
    AnaUtil::fillHist1D("nTightJet_5", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_5", hT, puevWt);
    AnaUtil::fillHist1D("sT_5", sT, puevWt);
    AnaUtil::fillHist1D("jdR_5", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_5", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_5", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_5", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_5", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_5", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_5", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_5", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_5", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_5", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_5", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_5", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_5", nbJets(), puevWt);
    AnaUtil::fillHist1D ("HppTrMass_5", HppTrMass, puevWt); 

    //6...) Cut (3 + j1j2DR > 1.6)
    //  if (j1j2DR > 1.6) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, puevWt);
    AnaUtil::fillHist1D("nTightJet_6", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_6", hT, puevWt);
    AnaUtil::fillHist1D("sT_6", sT, puevWt);
    //    AnaUtil::fillHist1D("jdR_6", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_6", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_6", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_6", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_6", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_6", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_6", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_6", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_6", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_6", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_6", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_6", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_6", nbJets(), puevWt);
    AnaUtil::fillHist1D("lepSumChr_6", (lep1c+lep2c), puevWt);


    //7...) Cut (4 + sameChargedLeptons)
    if (lep1c*lep2c < 0.0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, puevWt);
    AnaUtil::fillHist1D("nTightJet_7", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_7", hT, puevWt);
    AnaUtil::fillHist1D("sT_7", sT, puevWt);
    AnaUtil::fillHist1D("jdR_7", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_7", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_7", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_7", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_7", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_7", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_7", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_7", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_7", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_7", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_7", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_7", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_7", nbJets(), puevWt);
    AnaUtil::fillHist1D ("HppTrMass_7", HppTrMass, puevWt); 



    AnaUtil::fillHist1D("nTightJet_PS", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_PS", hT, puevWt);
    AnaUtil::fillHist1D("sT_PS", sT, puevWt);
    AnaUtil::fillHist1D("jdR_PS", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_PS", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_PS", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_PS", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_PS", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_PS", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_PS", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_PS", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_PS", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_PS", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_PS", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_PS", l2metDPhi, puevWt);


    //    dumpEvent(vz, false, true); //dump in detail


    //_________________________________________End of PreSelection_______________________________________________//    
      
    histf()->cd();

    //_____________________________________________MVA Skimming__________________________________________________//    

    if (skimObj_) {
      TreeVariables varList;
      varList.puevwt       = puevWt;
      varList.nJet         = TightJets.size();
      varList.met          = mt.met;
      varList.jet1Pt       = j1P4.Pt();
      varList.jet2Pt       = j2P4.Pt();
      varList.jet1Eta      = j1P4.Eta();
      varList.jet2Eta      = j2P4.Eta();
      varList.hT           = hT;
      varList.hTvec        = hTvec;

      varList.lep1Pt       = lep1pt;
      varList.lep1Eta      = lep1eta;
      varList.lep2Pt       = lep2pt;     
      varList.lep2Eta      = lep2eta;
      varList.lepSumCharge = (lep1c + lep2c);
      varList.lepDR        = dR;
      varList.lepSumPt     = lSumPt;
      varList.lT           = lT;
      varList.sT           = sT;
      varList.METsqrtST    = metSqrtST;

      varList.MetPhi       = mtPhi;
      varList.MetL1dPhi    = l1metDPhi;
      varList.MetL2dPhi    = l2metDPhi;
      varList.MetJ1dPhi    = j1metDPhi;
      varList.MetJ2dPhi    = j2metDPhi;
      varList.TrMass1      = mt_l1met;
      varList.TrMass2      = mt_l2met;
      varList.j1l1dR       = j1l1dr;
      varList.j1l2dR       = j1l2dr;
      varList.j2l1dR       = j2l1dr;
      varList.j2l2dR       = j2l2dr;
      varList.j1j2dR       = j1j2DR;
      varList.l1l2InvM     = Mass_2l;
      varList.lepDPhi      = lepDPhi;
      varList.jetDPhi      = j1j2DPhi;
      varList.mR           = MR;
 
      skimObj_->fill(varList);
    }

    //_________________________________________MVA Evaluation_______________________________________________//    
      
    histf()->cd();     
    //#if 0
    double mvaOut = -999.9;

    if (_readMVA) {
      InputVariables varlist;

      //varlist.nJet         = TightJets.size();
      varlist.met          = mt.met;
      varlist.hT           = hT;
      varlist.hTvec        = hTvec;
      varlist.lepDR        = dR;            
      varlist.lepSumPt     = lSumPt;
      //varlist.lT           = lT;
      varlist.MetL1dPhi    = l1metDPhi;
      varlist.MetL2dPhi    = l2metDPhi;
      varlist.MetJ1dPhi    = j1metDPhi;
      varlist.MetJ2dPhi    = j2metDPhi;
      varlist.j1l1dR       = j1l1dr;
      varlist.j1j2dR       = j1j2DR;
      varlist.l1l2InvM     = Mass_2l;


      //varlist.jet1Pt       = j1P4.Pt();
      //varlist.jet2Pt       = j2P4.Pt();
      //varlist.jet1Eta      = j1P4.Eta();
      //varlist.jet2Eta      = j2P4.Eta();
            
      
      //      varlist.lep1eta      = std::abs(lep1eta);
      //      varlist.lep2eta      = std::abs(lep2eta);
      //varlist.METsqrtST    = metSqrtST;     
      //      varlist.j2l1dR       = j2l1dr;
      //varlist.sT           = sT;


      //varlist.jetDPhi      = j1j2DPhi;
      //varlist.TrMass1      = mt_l1met;
      //
      //      varlist.TrMass2      = mt_l2met;
      //varlist.j1l1dR       = j1l1dR_;
      //varlist.j1l2dR       = j1l2dR_;
      //
      //varlist.j2l2dR       = j2l2dR_;

      mvaOut = _mvaObj->evaluate(_MVAnetwork, varlist);
    }

    histf()->cd();
    histf()->cd("TMVAnalysis");
    
    AnaUtil::fillHist1D ("mvaOutput", mvaOut, puevWt);    

    if (mvaOut < 0.0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 10);

    if (mvaOut < 0.025) continue;
    AnaUtil::fillHist1D("evtCutFlow", 11);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 11);

    if (mvaOut < 0.05) continue;
    AnaUtil::fillHist1D("evtCutFlow", 12);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12);

    if (mvaOut < 0.075) continue;
    AnaUtil::fillHist1D("evtCutFlow", 13);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 13);

    if (mvaOut < 0.1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 14);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 14);

    if (mvaOut < 0.125) continue;
    AnaUtil::fillHist1D("evtCutFlow", 15);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 15);

    if (mvaOut < 0.15) continue;
    AnaUtil::fillHist1D("evtCutFlow", 16);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 16);

    if (mvaOut < 0.16) continue;
    AnaUtil::fillHist1D("evtCutFlow", 17);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 17);

    if (mvaOut < 0.17) continue;
    AnaUtil::fillHist1D("evtCutFlow", 18);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 18);

    if (mvaOut < 0.175) continue;
    AnaUtil::fillHist1D("evtCutFlow", 19);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 19);

    if (mvaOut < 0.18) continue;
    AnaUtil::fillHist1D("evtCutFlow", 20);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 20);

    //    if (mvaOut < 0.2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 21);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 21);


    AnaUtil::fillHist1D ("HppTrMass_8", HppTrMass, puevWt);
    AnaUtil::fillHist1D("nTightJet_8", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_8", hT, puevWt);
    AnaUtil::fillHist1D("sT_8", sT, puevWt);
    AnaUtil::fillHist1D("jdR_8", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_8", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_8", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_8", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_8", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_8", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_8", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_8", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_8", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_8", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_8", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_8", l2metDPhi, puevWt);


    if (nbJets() > 0) continue; //Signal has no b-jet
    AnaUtil::fillHist1D("evtCutFlow", 22);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 22, puevWt);

    //9...) Cut (6 + tauVeto)
    if (nTau > 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 23);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 23, puevWt);

    AnaUtil::fillHist1D("nvtx_MVA1", ngoodVtx, puevWt);
    AnaUtil::fillHist1D("nTightJet_MVA1", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_MVA1", hT, puevWt);
    AnaUtil::fillHist1D("sT_MVA1", sT, puevWt);
    AnaUtil::fillHist1D("jdR_MVA1", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_MVA1", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_MVA1", jlDPhi, puevWt);
    AnaUtil::fillHist1D("jet1pt_MVA1", j1P4.Pt(), puevWt);//
    AnaUtil::fillHist1D("jet2pt_MVA1", j2P4.Pt(), puevWt);//
    AnaUtil::fillHist1D("jet1eta_MVA1", j1P4.Eta(), puevWt);//
    AnaUtil::fillHist1D("jet2eta_MVA1", j2P4.Eta(), puevWt);//
    AnaUtil::fillHist1D("lep1pt_MVA1", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_MVA1", lep2pt, puevWt);
    AnaUtil::fillHist1D("lep1eta_MVA1", lep1eta, puevWt);
    AnaUtil::fillHist1D("lep2eta_MVA1", lep2eta, puevWt);

    AnaUtil::fillHist1D("lepDR_MVA1", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_MVA1", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_MVA1", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_MVA1", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_MVA1", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_MVA1", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_MVA1", l2metDPhi, puevWt);
    AnaUtil::fillHist1D ("HppTrMass_MVA1", HppTrMass, puevWt);


    if (TightJets.size() < 4) continue;
    AnaUtil::fillHist1D("evtCutFlow", 24);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 24, puevWt);

    AnaUtil::fillHist1D("nvtx_MVA2", ngoodVtx, puevWt);
    AnaUtil::fillHist1D("nTightJet_MVA2", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_MVA2", hT, puevWt);
    AnaUtil::fillHist1D("sT_MVA2", sT, puevWt);
    AnaUtil::fillHist1D("jdR_MVA2", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_MVA2", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_MVA2", jlDPhi, puevWt);
    AnaUtil::fillHist1D("jet1pt_MVA2", j1P4.Pt(), puevWt);//
    AnaUtil::fillHist1D("jet2pt_MVA2", j2P4.Pt(), puevWt);//
    AnaUtil::fillHist1D("jet1eta_MVA2", j1P4.Eta(), puevWt);//
    AnaUtil::fillHist1D("jet2eta_MVA2", j2P4.Eta(), puevWt);//
    AnaUtil::fillHist1D("lep1pt_MVA2", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_MVA2", lep2pt, puevWt);
    AnaUtil::fillHist1D("lep1eta_MVA2", lep1eta, puevWt);
    AnaUtil::fillHist1D("lep2eta_MVA2", lep2eta, puevWt);

    AnaUtil::fillHist1D("lepDR_MVA2", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_MVA2", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_MVA2", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_MVA2", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_MVA2", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_MVA2", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_MVA2", l2metDPhi, puevWt);
    AnaUtil::fillHist1D ("HppTrMass_MVA2", HppTrMass, puevWt);

    if (LepCandList_[0].flavour == 1 && LepCandList_[1].flavour == 1) AnaUtil::fillHist1D("diffFlvYield", 0.0, puevWt);
    else if ((LepCandList_[0].flavour == 1 && LepCandList_[1].flavour == 2)||(LepCandList_[0].flavour == 2 && LepCandList_[1].flavour == 1)) AnaUtil::fillHist1D("diffFlvYield", 1.0, puevWt);
    else if (LepCandList_[0].flavour == 2 && LepCandList_[1].flavour == 2) AnaUtil::fillHist1D("diffFlvYield", 2.0, puevWt);


    // selEvLog() << run << " " << lumis << " " << event << std::endl;
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
}

bool MultiLeptonMVAna::SearchMinTwoUniqueJetPairs(const std::vector<std::pair<vhtm::Jet, vhtm::Jet>>& JetPair) {
  for (size_t i = 0; i < JetPair.size(); ++i){
    auto& pair1 = JetPair.at(i);
    TLorentzVector j1p4 = LL4JMETUtil::getP4(pair1.first);
    TLorentzVector j2p4 = LL4JMETUtil::getP4(pair1.second);
    for (size_t j = i+1; j < JetPair.size(); ++j){
      auto& pair2 = JetPair.at(j);
      TLorentzVector j3p4 = LL4JMETUtil::getP4(pair2.first);
      TLorentzVector j4p4 = LL4JMETUtil::getP4(pair2.second);
      if (j1p4 != j3p4 && j1p4 != j4p4 && j2p4 != j3p4 && j2p4 != j4p4)	return true;
    }
  }
  return false;
}

bool MultiLeptonMVAna::hasZcandidate(const std::vector<LeptonCandidate>& lepColl, double puWt){
  bool hasZToLL {false};
  bool hasZMass {false};
  for (size_t i = 0; i < lepColl.size(); ++i){
    auto& lep1 = lepColl.at(i);
    TLorentzVector TL1 = lep1.lRP4;
    double l1c = lep1.lCharge;
    for (size_t j = i+1; j < lepColl.size(); ++j){
      auto& lep2 = lepColl.at(j);
      TLorentzVector TL2 = lep2.lRP4;
      double l2c = lep2.lCharge;
      double lepInvM = (TL1+TL2).M();
      
      if(lep1.flavour == lep2.flavour){
	if (l1c * l2c < 0.0) {
	  AnaUtil::fillHist1D("oppChLepMass", lepInvM, puWt);
	  if (lepInvM > 20 && lepInvM < 160)  hasZToLL = true;
  	}
	else if (l1c * l2c > 0.0) {
	  AnaUtil::fillHist1D("sameChLepMass", lepInvM, puWt);
	  if (lepInvM > 85. && lepInvM < 110.) hasZMass = true;
	}
      }
      if (hasZToLL||hasZMass) return true;
    }
  }
  return false;
}

void MultiLeptonMVAna::endJob() {
  PhysicsObjSelector::endJob();
  
  histf()->cd();
  histf()->cd("TMVAnalysis");
  vector<string> evLabels {
    "0) Events processed: ",
      "1) Have good Vtx: ",
      "2) triggered: ",
      "3) nLep = 2: ",
      "4) nTightJet >= 2:",
      "5) jet1pt > 60(NA): ",
      "6) jet2pt > 40(NA): ",
      "7) met > 30: ",
      "8) j1j2dR < 1.8(NA)", 
      "9) same chr leptons:",
      "10)bdt >= 0.000: ",
      "11)bdt >= 0.025: ",
      "12)bdt >= 0.050: ",
      "13)bdt >= 0.075: ",
      "14)bdt >= 0.100: ",
      "15)bdt >= 0.125: ",
      "16)bdt >= 0.150: ",
      "17)bdt >= 0.160 : ",
      "18)bdt >= 0.170: ",
      "19)bdt >= 0.175: ",
      "20)bdt >= 0.180: ",
      "21)bdt >= 0.200(NA): ",
      "22)no bjet: ",
      "23)no tau:",
      "24)nJet >= 4:"
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
    //cout<<"LumiFac: "<<lumiFac<<endl; ////////////////////////////////////////
    //LL4JMETUtil::scaleHistogram("evtCutFlowWt", lumiFac);
    LL4JMETUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");  
  }
}

void MultiLeptonMVAna::closeFiles() {
  AnaBase::closeFiles();
  // Take care of local stuff first                                                
  //  if (_mvaObj != nullptr) _mvaObj->close();
  if (skimObj_ != nullptr) skimObj_->close();
  if (syncDumpf_.is_open()) syncDumpf_.close();
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
bool MultiLeptonMVAna::readJob(const string& jobFile, int& nFiles)
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
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    //else if (key == "mvaInputFile")
    // _mvaInputFile = value;
    else if (key == "MVAnetwork")
      _MVAnetwork = value;
    else if (key == "MVAxmlFile")
      _MVAxmlFile = value;
    else if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
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
  /*
#ifdef SKIP_DUPLICATE_FROM_FILE
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
  */
  printJob();
  
  //if (readGenInfo()) genAna_ = std::make_unique<GenAnalysis>();
  //return true;
}
void MultiLeptonMVAna::printJob(ostream& os) const
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
  //  if (isMC()) 
  //  os << "    dumpGenInfo: " << std::boolalpha << dumpGenInfo_ << endl;

  //AnaUtil::showList(eventFilelist_, ">>> INFO. Input event files:", os);
}
