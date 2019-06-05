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
#include "MultiLeptonCUTAna.h"
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
MultiLeptonCUTAna::MultiLeptonCUTAna()
  : PhysicsObjSelector()
   {}

// ----------
// Destructor
// ----------
MultiLeptonCUTAna::~MultiLeptonCUTAna() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MultiLeptonCUTAna::beginJob() 
{ 
  if (!PhysicsObjSelector::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("CUTAnalysis");
  
  bookHistograms();

  if (readGenInfo() && genAna_ == nullptr) {
    genAna_ = std::make_unique<GenAnalysis>();
    genAna_->bookHistograms(histf());
  }

#ifdef  SKIP_DUPLICATE_ALL
  eventIdStore_.clear();
#endif

  return true;
}

// ---------------
// Book histograms
// ---------------
void MultiLeptonCUTAna::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("CUTAnalysis");

  // book histograms to be filled at different stages


  new TH1D("nvtx_0", "Number of Good vertices", 60, 0, 60);
  new TH1D("nvtx_1", "Number of Good vertices", 60, 0, 60);
  new TH1D("nvtx_2", "Number of Good vertices", 60, 0, 60);
  new TH1D("nvtx_3", "Number of Good vertices", 60, 0, 60);
  new TH1D("nvtx_4", "Number of Good vertices", 60, 0, 60);
  new TH1D("nvtx_5", "Number of Good vertices", 60, 0, 60);
  new TH1D("nvtx_6", "Number of Good vertices", 60, 0, 60);
  if (isMC()) {
    new TH1D("puweight", "PU reweight factor", 1000, -20., 20.);
    new TH1D("pu_evtWt", "PU nd Evt Wt", 1000, -20., 20.);
    new TH1D("evtweight", "Event weight factor (MC events)", 20, -10., 10.);
  }

  //------- Object PLots -----------------------------------------------
  
  new TH1D("evtCutFlow", "Event CutFlow", 19, -0.5, 18.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 19, -0.5, 18.5);
  new TH1D("met", "Missing Transver Energy",200, 0, 200);
  new TH1D("corrmet", "Missing Transver Energy", 600, 0., 400.);
  new TH1D("puppimet", "Missing Transver Energy", 200, 0, 200);
  new TH1D("yield" ,"", 1, -0.5, 0.5);
  if (isMC()) new TH1D("yieldWt", "", 1, -0.5, 0.5);
  new TH1D("GenMass" ,"", 800, 0., 1000.);
  new TH1D("HMassT" ,"", 800, 0., 1000.);
  new TH1D("HMassL" ,"", 800, 0., 1000.);
  new TH1D("GenMatchedMass" ,"", 800, 0., 1000.);
  new TH1D("HMassSig" ,"", 800, 0., 1000.);
  new TH1D("HTrMass" ,"", 800, 0., 1000.);

  new TH1D("nGenLep" ,"", 5, -0.5, 4.5);
  new TH1D("nGenNeu" ,"", 5, -0.5, 4.5);
  new TH1D("nGenJets" ,"", 10, -0.5, 9.5);
  new TH1D("GenLepPt" ,"", 400, 0., 400.);
  new TH1D("GenLepEta" ,"", 500, -5.0, 5.0);
  new TH1D("GenLepMass" ,"", 250, 0.0, 500.0);
  new TH1D("GenHMass" ,"", 250, 0.0, 500.0);
  new TH1D("GenLepSumCharge" ,"", 5, -2.5, 2.5);
  new TH1D("GenJetPt" ,"", 400, 0., 400.);
  new TH1D("GenJetEta" ,"", 500, -5.0, 5.0);
  new TH1D("missedQpt" ,"", 400, 0., 200.);
  new TH1D("missedQeta" ,"", 500, -5.0, 5.0);

  new TH1D("HPPpt" ,"", 1000, 0., 1000.0);
  new TH1D("HPPenergy" ,"", 1000, 0., 1000.0);
  new TH1D("HMMpt" ,"", 1000, 0., 1000.0);
  new TH1D("HMMenergy" ,"", 1000, 0., 1000.0);
  new TH1D("HPpt" ,"", 1000, 0., 1000.0);
  new TH1D("HPenergy" ,"", 1000, 0., 1000.0);
  new TH1D("HMpt" ,"", 1000, 0., 1000.0);
  new TH1D("HMenergy" ,"", 1000, 0., 1000.0);

  new TH1D("hpphmmDR" ,"", 200, 0., 5.0);
  new TH1D("hpphmDR" , "", 200, 0., 5.0);
  new TH1D("hmmhpDR" , "", 200, 0., 5.0);
  new TH1D("hpphmmDPhi" ,"", 200, 0., 4.0);
  new TH1D("hpphmDPhi" , "", 200, 0., 4.0);
  new TH1D("hmmhpDPhi" , "", 200, 0., 4.0);


  //After Trigger Cut
  new TH1D("nLepCand_0" ,"", 5, -0.5, 4.5);
  new TH1D("nTightJet_0" ,"", 10, -0.5, 9.5);
  new TH1D("nLooseJet_0" ,"", 10, -0.5, 9.5);
  new TH1D("corrmet_0", "Missing Transver Energy", 600, 0., 400.);
  new TH1D("nbJets_0" ,"", 10, -0.5, 9.5);

  //1. nLep >= 2
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

  //2. (1 + nJet>=2)
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
  new TH1D("j1l1DPhi_2", "", 400, 0., 4.);
  new TH1D("j1l1DR_2", "", 500, 0., 5.);



  //3. (2 + jet1Pt > 60)
  new TH1D("nTightJet_3" ,"", 10, -0.5, 9.5);
  new TH1D("hT_3" ,"", 1000, 0., 1000.0);
  new TH1D("sT_3" ,"", 1000, 0., 1000.0);
  new TH1D("jet2pt_3" ,"", 300, 0., 300.0);
  new TH1D("corrmet_3", "Missing Transver Energy", 600, 0., 400.);
  new TH1D("l1l2InvM_3", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_3", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_3", "", 500, 0., 5.);
  new TH1D("nbJets_3" ,"", 10, -0.5, 9.5);
  new TH1D("j1l1DPhi_3", "", 400, 0., 4.);
  new TH1D("j1l1DR_3", "", 500, 0., 5.);


  //4. (3 + jet2Pt > 40)
  new TH1D("nTightJet_4" ,"", 10, -0.5, 9.5);
  new TH1D("hT_4" ,"", 1000, 0., 1000.0);
  new TH1D("sT_4" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_4" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_4" ,"", 300, 0., 300.0);
  new TH1D("lepDR_4" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_4" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_4", "",100,0.0, 5.0);
  new TH1D("jdPhi_4", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_4", "",500,0.0, 5.0);
  new TH1D("lSumPt_4" ,"", 500, 0., 500.0);
  new TH1D("lT_4" ,"", 1000, 0., 1000.0);
  new TH1D("corrmet_4", "Missing Transver Energy", 600, 0., 400.);
  new TH1D("l1l2InvM_4", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_4", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_4", "", 500, 0., 5.);
  new TH1D("nbJets_4" ,"", 10, -0.5, 9.5);
  new TH1D("j1l1DPhi_4", "", 400, 0., 4.);
  new TH1D("j1l1DR_4", "", 500, 0., 5.);


  //5. (4 + met > 60)

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
  new TH1D("j1l1DPhi_5", "", 400, 0., 4.);
  new TH1D("j1l1DR_5", "", 500, 0., 5.);



  //6. (5 + hT > 160)

  new TH1D("nTightJet_6" ,"", 10, -0.5, 9.5);
  //  new TH1D("hT_6" ,"", 1000, 0., 1000.0);
  new TH1D("sT_6" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_6" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_6" ,"", 300, 0., 300.0);
  new TH1D("lepDR_6" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_6" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_6", "",100,0.0, 5.0);
  new TH1D("jdPhi_6", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_6", "",500,0.0, 5.0);
  new TH1D("lSumPt_6" ,"", 500, 0., 500.0);
  new TH1D("lT_6" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_6", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_6", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_6", "", 500, 0., 5.);
  new TH1D("nbJets_6" ,"", 10, -0.5, 9.5);
  new TH1D("j1l1DPhi_6", "", 400, 0., 4.);
  new TH1D("j1l1DR_6", "", 500, 0., 5.);


  //7. (6 + lepSumPt > 60)

  new TH1D("nTightJet_7" ,"", 10, -0.5, 9.5);
  new TH1D("hT_7" ,"", 1000, 0., 1000.0);
  new TH1D("sT_7" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_7" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_7" ,"", 300, 0., 300.0);
  new TH1D("lepDR_7" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_7" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_7", "",100,0.0, 5.0);
  new TH1D("jdPhi_7", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_7", "",500,0.0, 5.0);
  //  new TH1D("lSumPt_6" ,"", 500, 0., 500.0);
  new TH1D("lT_7" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_7", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_7", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_7", "", 500, 0., 5.);
  new TH1D("nbJets_7" ,"", 10, -0.5, 9.5);
  new TH1D("j1l1DPhi_7", "", 400, 0., 4.);
  new TH1D("j1l1DR_7", "", 500, 0., 5.);


  //8. (7 + j1j2dR < 1.5)

  new TH1D("nTightJet_8" ,"", 10, -0.5, 9.5);
  new TH1D("hT_8" ,"", 1000, 0., 1000.0);
  new TH1D("sT_8" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_8" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_8" ,"", 300, 0., 300.0);
  new TH1D("lepDR_8" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_8" ,"", 500, 0.0, 5.0);
  //  new TH1D("jdR_7", "",100,0.0, 5.0);
  new TH1D("jdPhi_8", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_8", "",500,0.0, 5.0);
  new TH1D("lSumPt_8" ,"", 500, 0., 500.0);
  new TH1D("lT_8" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_8", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_8", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_8", "", 500, 0., 5.);
  new TH1D("nbJets_8" ,"", 10, -0.5, 9.5);
  new TH1D("lepSumChr_8", "", 5, -2.5, 2.5);
  new TH1D("j1l1DPhi_8", "", 400, 0., 4.);
  new TH1D("j1l1DR_8", "", 500, 0., 5.);


  //9. (8 + SameChrLeps)

  new TH1D("nTightJet_9" ,"", 10, -0.5, 9.5);
  new TH1D("hT_9" ,"", 1000, 0., 1000.0);
  new TH1D("sT_9" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_9" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_9" ,"", 300, 0., 300.0);
  new TH1D("lepDR_9" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_9" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_9", "",100,0.0, 5.0);
  new TH1D("jdPhi_9", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_9", "",500,0.0, 5.0);
  new TH1D("lSumPt_9" ,"", 500, 0., 500.0);
  new TH1D("lT_9" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_9", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_9", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_9", "", 500, 0., 5.);
  new TH1D("nbJets_9" ,"", 10, -0.5, 9.5);
  new TH1D("j1l1DPhi_9", "", 400, 0., 4.);
  new TH1D("j1l1DR_9", "", 500, 0., 5.);


  //10. (9 + hasnoZToLL)

  new TH1D("nTightJet_10" ,"", 10, -0.5, 9.5);
  new TH1D("hT_10" ,"", 1000, 0., 1000.0);
  new TH1D("sT_10" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_10" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_10" ,"", 300, 0., 300.0);
  new TH1D("lepDR_10" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_10" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_10", "",100,0.0, 5.0);
  new TH1D("jdPhi_10", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_10", "",500,0.0, 5.0);
  new TH1D("lSumPt_10" ,"", 500, 0., 500.0);
  new TH1D("lT_10" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_10", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_10", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_10", "", 500, 0., 5.);
  new TH1D("nbJets_10" ,"", 10, -0.5, 9.5);
  new TH1D("j1l1DPhi_10", "", 400, 0., 4.);
  new TH1D("j1l1DR_10", "", 500, 0., 5.);


  //11. (10 + 0.4<lepDR<1.5)
  new TH1D("nTightJet_11" ,"", 10, -0.5, 9.5);
  new TH1D("hT_11" ,"", 1000, 0., 1000.0);
  new TH1D("sT_11" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_11" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_11" ,"", 300, 0., 300.0);
  //  new TH1D("lepDR_11" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_11" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_11", "",100,0.0, 5.0);
  new TH1D("jdPhi_11", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_11", "",500,0.0, 5.0);
  new TH1D("lSumPt_11" ,"", 500, 0., 500.0);
  new TH1D("lT_11" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_11", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_11", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_11", "", 500, 0., 5.);
  new TH1D("nbJets_11" ,"", 10, -0.5, 9.5);
  new TH1D("j1l1DPhi_11", "", 400, 0., 4.);
  new TH1D("j1l1DR_11", "", 500, 0., 5.);


  //12. (11 + lep1MetDPhi < 0.7)
  new TH1D("nbJets_12" ,"", 10, -0.5, 9.5);
  new TH1D("nTightJet_12" ,"", 10, -0.5, 9.5);
  new TH1D("hT_12" ,"", 1000, 0., 1000.0);
  new TH1D("sT_12" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_12" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_12" ,"", 300, 0., 300.0);
  new TH1D("lepDR_12" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_12" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_12", "",100,0.0, 5.0);
  new TH1D("jdPhi_12", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_12", "",500,0.0, 5.0);
  new TH1D("lSumPt_12" ,"", 500, 0., 500.0);
  new TH1D("lT_12" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_12", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_12", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_12", "", 500, 0., 5.);
  new TH1D("j1l1DPhi_12", "", 400, 0., 4.);
  new TH1D("j1l1DR_12", "", 500, 0., 5.);


  //13. (12 + bVeto)
  //  new TH1D("nbJets_13" ,"", 10, -0.5, 9.5);
  new TH1D("nTightJet_13" ,"", 10, -0.5, 9.5);
  new TH1D("hT_13" ,"", 1000, 0., 1000.0);
  new TH1D("sT_13" ,"", 1000, 0., 1000.0);

  new TH1D("lep1pt_13" ,"", 300, 0., 300.0);
  new TH1D("lep2pt_13" ,"", 300, 0., 300.0);
  //  new TH1D("lepDR_11" ,"", 200, 0., 5.0);
  new TH1D("lepDPhi_13" ,"", 500, 0.0, 5.0);
  new TH1D("jdR_13", "",100,0.0, 5.0);
  new TH1D("jdPhi_13", "",500,0.0, 5.0);
  new TH1D("JetLepDPhi_13", "",500,0.0, 5.0);
  new TH1D("lSumPt_13" ,"", 500, 0., 500.0);
  new TH1D("lT_13" ,"", 1000, 0., 1000.0);

  new TH1D("l1l2InvM_13", "InvMass of 2leptons", 400, 0., 400.);
  new TH1D("lep1MetDPhi_13", "", 500, 0., 5.);
  new TH1D("lep2MetDPhi_13", "", 500, 0., 5.);
  new TH1D("j1l1DPhi_13", "", 400, 0., 4.);
  new TH1D("j1l1DR_13", "", 500, 0., 5.);


  //14. (13 + tauVeto)

  if (isMC())  new TH1D("diffFlvYield", "", 3, -0.5, 2.5);










  new TH1D("sigmaPToverPT_muon", "", 500, 0., 5.0);
  new TH1D("sigmaPToverPT_ele" , "", 500, 0., 5.0);
  new TH1D("isGsfCtfScPixCharge", "", 2, -0.5, 1.5);
  new TH1D("isGsfScPixCharge", "", 2, -0.5, 1.5);


  histf()->cd();
  histf()->ls();
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------

void MultiLeptonCUTAna::clearLists() {
  PhysicsObjSelector::clear();
  vtxList_.clear();
  LepCandList_.clear();
  JetPairSelected.clear();
  JetSelected.clear();
  SignalJetsSet.clear();
  SignalJets.clear();
  GenJets.clear();
  GenLeptons.clear();
  GenMet.clear();s
}

// -------------------
// The main event loop
// -------------------

void MultiLeptonCUTAna::eventLoop()
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
  //  TH1D *wtDiff = new TH1D ("wtDiff", "", 20, 0., 20.);  
 
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
    histf()->cd("CUTAnalysis");
    
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
    histf()->cd("CUTAnalysis");
    AnaUtil::fillHist1D("evtCutFlow", 0);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, puevWt); //Events processed





    /////////////////////////////////////////////////////////GenFilter&Study////////////////////////////////////////////////////////
    //    cout<<"Ev:::::::::::::::::::::::::::::::::::::::::::"<<ev<<endl;
#if 0
    if (genAna_ != nullptr) {
      genAna_->setEvent(genParticleColl());
      //      genAna_ -> analyze(histf());            
      bool genOk = genAna_->filter();
      if (!genOk) continue;
    }
#endif
    /*  
    if (isMC() && readGenInfo()){
      //      AnaUtil::fillHist1D("nGenJets_0", genJetColl()->size(), puevWt);
      AnaUtil::fillHist1D("GenJet1Pt_0", (genJetColl()->at(0)).pt, puevWt);
      AnaUtil::fillHist1D("GenJet2Pt_0", (genJetColl()->at(1)).pt, puevWt);
      //      if (genAna_ != nullptr && dumpGenInfo_) genAna_->dumpEvent();      
      //if (genAna_ != nullptr && dumpGenInfo_) genAna_->AnaBase::dumpGenInfo();      
      //dumpGenInfo();
      int njet = 0;      
      for (const auto& genj: *genJetColl()){
	cout<<"jetPt: "<<genj.pt<<"\t"<<"jet eta: "<<genj.eta<<"\t"<<"jet phi: "<<genj.phi<<"\t"<<"jetEnergy: "<<genj.energy<<std::endl;
	if (!(genj.pt > 25 && std::abs(genj.eta) <= 2.5)) continue;
        njet++;	
	//cout<<"1st"<<endl;
	TLorentzVector jp4 = LL4JMETUtil::getP4(genj);
	int mmid = -1;
	TLorentzVector gp4;
	gp4.SetPtEtaPhiE(0.,0.,0.,0.);
	bool isnotcleaned = false;
	//bool isfromV      = false;
	for (const auto& genp: *genParticleColl()){
	  int pdgid = std::abs(genp.pdgId);
	  int status = genp.status;
	  double genEta = std::abs(genp.eta);
	  if (!((pdgid == 11 || pdgid == 13) && genp.pt >= 10. && status == 1 && genEta < 2.5)) continue;
	  int index = getMotherId(genp, mmid);
	  if (index < 0) continue;
	  int lepMomId = std::abs((genParticleColl()->at(index)).pdgId);
	  if (lepMomId != 24) continue;
	  gp4 = LL4JMETUtil::getP4(genp);
	  AnaUtil::fillHist1D("genJlepDR", jp4.DeltaR(gp4), puevWt);
	  if (jp4.DeltaR(gp4) < 0.4) {
	    isnotcleaned = true;
	    break;
	  }
	}
	if (isnotcleaned) continue;
	GenJets.push_back(genj);
      }
      AnaUtil::fillHist1D("nGenJets_0", njet, puevWt);
      AnaUtil::fillHist1D("nGenJets", GenJets.size(), puevWt);
      if (GenJets.size() >= 2){
	AnaUtil::fillHist1D("GenJet1Pt", GenJets[0].pt, puevWt);
	AnaUtil::fillHist1D("GenJet2Pt", GenJets[1].pt, puevWt);
      }
      if (GenJets.size() < 4) continue;
      int ind = 0;
      TLorentzVector gjp4;
      gjp4.SetPtEtaPhiE(0.,0.,0.,0.);
      for (auto& gj: GenJets){
	ind++;
	gjp4 += LL4JMETUtil::getP4(gj);
	if (ind == 4) break;
      }
      AnaUtil::fillHist1D("GenMass", gjp4.M(), puevWt);
    }
*/


    //////////////////////////////////////////GenLepton//////////////////////////////////////////////
    if (isMC() & readGenInfo()){
      int nLep = 0;
      TLorentzVector HPPp4;
      HPPp4.SetPtEtaPhiE(0.,0.,0.,0.);
      TLorentzVector HMMp4;
      HMMp4.SetPtEtaPhiE(0.,0.,0.,0.);
      TLorentzVector HPp4;
      HPp4.SetPtEtaPhiE(0.,0.,0.,0.);
      TLorentzVector HMp4;
      HMp4.SetPtEtaPhiE(0.,0.,0.,0.);
      bool hasHPP = false;
      bool hasHMM = false;
      bool hasHP  = false;
      bool hasHM  = false;    
      for (const auto& genp: *genParticleColl()){
	int pdgid = std::abs(genp.pdgId);
	int status = genp.status;
	double genEta = std::abs(genp.eta);
	int mmid = -1;
	int _mmid = -1;
	
	if (genp.pdgId == 9900041 && status == 22){
	  AnaUtil::fillHist1D ("HPPpt", genp.pt, 1.0);
	  AnaUtil::fillHist1D ("HPPenergy", genp.energy, 1.0);
	  HPPp4 = LL4JMETUtil::getP4(genp);
	  hasHPP = true;
	}
        else if (genp.pdgId == -9900041 && status == 22){
	  AnaUtil::fillHist1D ("HMMpt", genp.pt, 1.0);
	  AnaUtil::fillHist1D ("HMMenergy", genp.energy, 1.0);
	  HMMp4 = LL4JMETUtil::getP4(genp);
	  hasHMM = true;
	}
        else if (genp.pdgId == 37 && status == 22){
	  AnaUtil::fillHist1D ("HPpt", genp.pt, 1.0);
	  AnaUtil::fillHist1D ("HPenergy", genp.energy, 1.0);
	  HPp4 = LL4JMETUtil::getP4(genp);
	  hasHP = true;
	}
	else if (genp.pdgId == -37 && status == 22){
	  AnaUtil::fillHist1D ("HMpt", genp.pt, 1.0);
	  AnaUtil::fillHist1D ("HMenergy", genp.energy, 1.0);
	  HMp4 = LL4JMETUtil::getP4(genp);
	  hasHM = true;
	}
	if (hasHPP && hasHMM) {
	  AnaUtil::fillHist1D ("hpphmmDR", HPPp4.DeltaR(HMMp4), 1.0);
	  AnaUtil::fillHist1D ("hpphmmDPhi", std::fabs(TVector2::Phi_mpi_pi(HPPp4.Phi()-HMMp4.Phi())), 1.0);
	}
	else if (hasHPP && hasHM) {
	  AnaUtil::fillHist1D ("hpphmDR",  HPPp4.DeltaR(HMp4),  1.0);
	  AnaUtil::fillHist1D ("hpphmDPhi",  std::fabs(TVector2::Phi_mpi_pi(HPPp4.Phi()-HMp4.Phi())),  1.0);
	}
	else if (hasHMM && hasHP) {
	  AnaUtil::fillHist1D ("hmmhpDR",  HMMp4.DeltaR(HPp4),  1.0);
	  AnaUtil::fillHist1D ("hmmhpDPhi",  std::fabs(TVector2::Phi_mpi_pi(HMMp4.Phi()-HPp4.Phi())),  1.0);
	}

	if ((pdgid == 9000014 || pdgid == 9000012 || pdgid == 9000016) && status == 1){
	  int index = getMotherId(genp, _mmid);
	  auto &mom = genParticleColl()->at(index);
	  int lepMomId = std::abs(mom.pdgId);
	  if (lepMomId != 24) continue;
	  index = getMotherId(mom, _mmid);
	  auto &gmom = genParticleColl()->at(index); 
	  int lepGMomId = std::abs(gmom.pdgId); 
	  if (lepGMomId == 9900041) {
	    HPPneu.push_back(genp);
	  }
	  else if (lepGMomId == -9900041) {
	    HMMneu.push_back(genp);
	  }
	  else if (lepGMomId == 37) {
	    HPneu.push_back(genp);
	  }
	  else if (lepGMomId == -37) {
	    HMneu.push_back(genp);
	  }

	}
	if (!((pdgid == 11 || pdgid == 13) && status == 1)) continue;
	int index = getMotherId(genp, mmid);
	if (index < 0) continue;
	auto &mom = genParticleColl()->at(index);
	int lepMomId = std::abs(mom.pdgId);
	if (lepMomId != 24) continue;
	index = getMotherId(mom, mmid);
	if (index < 0) continue;
	auto &gmom = genParticleColl()->at(index);
	int lepGMomId = std::abs(gmom.pdgId);
	if (lepGMomId != 9900041) continue;
	AnaUtil::fillHist1D("GenLepPt", genp.pt, 1.0);
	AnaUtil::fillHist1D("GenLepEta", genp.eta, 1.0);
	if (!(genEta <= 2.5 && genp.pt > 10)) continue;
	GenLeptons.push_back(genp);
      }
      std::sort(std::begin(GenLeptons), std::end(GenLeptons), PtComparator<GenParticle>()); //sorting lepton candidates        
      AnaUtil::fillHist1D("nGenLep", GenLeptons.size(), 1.0);
      AnaUtil::fillHist1D("nGenNeu", GenMet.size(), 1.0);
      if (GenLeptons.size() != 2) continue;
      TLorentzVector gp4;
      gp4.SetPtEtaPhiE(0.,0.,0.,0.);
      TLorentzVector mp4;
      mp4.SetPtEtaPhiE(0.,0.,0.,0.);
      double _charge = 0.0;
      for (auto & gen: GenLeptons){
	gp4 += LL4JMETUtil::getP4(gen); 
	_charge += gen.charge;
      }
      for (auto & neu: GenMet){
	mp4 += LL4JMETUtil::getP4(neu); 
      }
      AnaUtil::fillHist1D("GenLepMass", gp4.M(), 1.0);
      AnaUtil::fillHist1D("GenHMass", (gp4+mp4).M(), 1.0);
      AnaUtil::fillHist1D("GenLepSumCharge", _charge, 1.0);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////

    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, puevWt); //Gen Filter (for lepton only)
    

    //////////////////////////////////////////////GenJets////////////////////////////////////////////
    if (isMC() & readGenInfo()){
      int njet = 0;      
      for (const auto& genj: *genJetColl()){
	if (!(genj.pt > 20 && std::abs(genj.eta) < 4.7)) continue;
	njet++;	
	//cout<<"1st"<<endl;
	TLorentzVector jp4 = LL4JMETUtil::getP4(genj);
	int mmid = -1;
	TLorentzVector gp4;
	gp4.SetPtEtaPhiE(0.,0.,0.,0.);
	for (const auto& genp: *genParticleColl()){
	  int pdgid = std::abs(genp.pdgId);
	  int status = genp.status;
	  double genEta = std::abs(genp.eta);
	  //  if (!(((pdgid >= 1 && pdgid <= 4) || pdgid == 21) && status == 71)) continue;
	  if (!((pdgid >= 1 && pdgid <= 4) && status == 23)) continue;
	  //cout<<"2nd"<<endl;
	  vector<int> m = genp.motherIndices;
	  if (m.size() != 1) continue;
	  //	  int index = getMotherIdForQ(genp, mmid);
	  int index = m.at(0);
	  if (index < 0) continue;
	  int qMomId = std::abs((genParticleColl()->at(index)).pdgId);
	  if (!(qMomId == 24 || qMomId == 23)) continue;
	  //cout<<"3rd"<<endl;
	  TLorentzVector gp4 = LL4JMETUtil::getP4(genp);
	  if (jp4.DeltaR(gp4) < 0.4){
	    GenJets.push_back(genj);
	    //	    AnaUtil::fillHist1D("genPt", genj.pt, 1.0);
	    //isfromV = true;
	    break;
	  }
	  else {
	    AnaUtil::fillHist1D("missedQpt", genp.pt, 1.0);
	    AnaUtil::fillHist1D("missedQeta", genp.eta, 1.0);
	  }
	}
      }
      AnaUtil::fillHist1D("nGenJets_0", njet, puevWt);
      AnaUtil::fillHist1D("nGenJets", GenJets.size(), puevWt);
      for (auto &j: GenJets){
	AnaUtil::fillHist1D("GenJetPt", j.pt, puevWt);
	AnaUtil::fillHist1D("GenJetEta", j.eta, puevWt);
      }
	
      if (GenJets.size() >= 2){
	AnaUtil::fillHist1D("GenJet1Pt", GenJets[0].pt, puevWt);
	AnaUtil::fillHist1D("GenJet2Pt", GenJets[1].pt, puevWt);
      }
      if (GenJets.size() < 4) continue;
    }
    /////////////////////////////////////////////////////////////////////////////////////////
      
    
    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, puevWt); //Gen Filter (for gen jets)
    

    // at least 1 good PV
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    size_t ngoodVtx = vtxList_.size();
    AnaUtil::fillHist1D("nvtx_0", ngoodVtx, puevWt);
    if (ngoodVtx < 1) continue;

    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, puevWt); //good vtx


    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt); 
    
    // Is event triggered?
    if (0) dumpTriggerPaths(std::cout, true);
    if (useTrigger() && !isTriggered(true, false)) continue;
    
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, puevWt); //is triggered
    
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
    
    const auto& LooseJets = getLooseJetList();
    
    //********************** A N A L Y S I S **********************//
    
    histf()->cd();//Required
    histf()->cd("CUTAnalysis");

    AnaUtil::fillHist1D("met", metColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("corrmet", corrmetColl()->at(0).met, puevWt);
    AnaUtil::fillHist1D("puppimet", puppimetColl()->at(0).met, puevWt);

    //P A C K I N G  L E P T O N S in LepCandList_
    if (nMuon > 0) packLeptons<vhtm::Muon>(muPhotonPairList, LepCandList_); 
    if (nEle > 0)  packLeptons<vhtm::Electron>(elePhotonPairList, LepCandList_);
    AnaUtil::fillHist1D("nLeptonCand", LepCandList_.size(), puevWt);

    addLeptonIsolation(LepCandList_, elePhotonPairList, muPhotonPairList); //add lepton isolation to lepcandidates    
    std::sort(std::begin(LepCandList_), std::end(LepCandList_), PtComparatorLep<LeptonCandidate>()); //sorting lepton candidates  
    //std::sort(std::begin(LooseJets), std::end(LooseJets), PtComparator<vhtm::Jet>());
    histf()->cd();
    histf()->cd("CUTAnalysis");

    //    cout<<"Event no: "<<ev<<endl;

    //_________________________________________PreSelection Begins_______________________________________________//




    //test/////////////////
    /*
    for (auto &mu: muPhotonPairList){
      AnaUtil::fillHist1D("sigmaPToverPT_muon", (mu.first).sigmaPToverPT, puevWt);
    }
    for (auto &el: elePhotonPairList){
      AnaUtil::fillHist1D("sigmaPToverPT_ele",  (el.first).sigmaPToverPT, puevWt);
      AnaUtil::fillHist1D("isGsfCtfScPixCharge",(el.first).isGsfCtfScPixChargeConsistent, puevWt);
      AnaUtil::fillHist1D("isGsfScPixCharge",   (el.first).isGsfScPixChargeConsistent, puevWt);
      }*/

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

 
    //1...) Cut (nLep >= 2)
    if (LepCandList_.size() != 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, puevWt);

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

    double transMass = (l1p4+l2p4).Mt();
    double mt_l1met = TMath::Sqrt(2*lep1pt*(mt.met)*(1-TMath::Cos(l1metDPhi)));
    double mt_l2met = TMath::Sqrt(2*lep2pt*(mt.met)*(1-TMath::Cos(l2metDPhi)));
    AnaUtil::fillHist1D("l1TrMass_1", mt_l1met, puevWt);
    AnaUtil::fillHist1D("l2TrMass_1", mt_l2met, puevWt);

    AnaUtil::fillHist1D("nTightJet_1", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("corrmet_1",  mt.met, puevWt);
    AnaUtil::fillHist1D("nbJets_1", nbJets(), puevWt);
   


    ///////////////////////////////////////////////GenLevelTest////////////////////////////////////////////////////////
    if (isMC() && readGenInfo()){
      int ind = 0;
      int ind_ = 0;
      TLorentzVector gjp4;
      gjp4.SetPtEtaPhiE(0.,0.,0.,0.);
      for (auto& gj: GenJets){
	ind++;
	gjp4 += LL4JMETUtil::getP4(gj);
	if (ind == 4) break;
      }
      AnaUtil::fillHist1D("GenMass", gjp4.M(), 1.0);
    
      if (LooseJets.size() >= 4){
	TLorentzVector j_p4;
	j_p4.SetPtEtaPhiE(0.,0.,0.,0.);
	for (size_t i = 0; i < LooseJets.size(); ++i){
	  auto &j = LooseJets[i];
	  TLorentzVector jp4 = LL4JMETUtil::getP4(j);
	  for (auto& g: GenJets){
	    TLorentzVector gp4 = LL4JMETUtil::getP4(g);
	    if (jp4.DeltaR(gp4) < 0.4){
	      j_p4 += jp4;
	      ind_++;
	      break;
	    }
	    if (ind_ == 4) AnaUtil::fillHist1D("GenMatchedMass", j_p4.M(), puevWt);
	  }
	}
      }
    }//isMC loop finished

    TLorentzVector ljp4;
    ljp4.SetPtEtaPhiE(0.,0.,0.,0.);
    if (LooseJets.size() >= 4){
      for (size_t i = 0; i < 4; ++i){
	auto &j = LooseJets[i];
	ljp4 += LL4JMETUtil::getP4(j);
      }
      AnaUtil::fillHist1D("HMassL", ljp4.M(), puevWt);
    }
    for (size_t i1 = 0; i1 < LooseJets.size(); ++i1){	
      const Jet& jet1 = LooseJets.at(i1);
      TLorentzVector jet1p4 = LL4JMETUtil::getP4(jet1);
      for (size_t i2 = i1 + 1; i2 < LooseJets.size(); ++i2) {
	const Jet& jet2 = LooseJets.at(i2);
	TLorentzVector jet2p4 = LL4JMETUtil::getP4(jet2);
	double wMass = (jet1p4+jet2p4).M();
	//AnaUtil::fillHist1D("JetPairMass", wMass, puevWt);
	//	if (wMass > 65. && wMass < 95.) {
	if (wMass > 65. && wMass < 95.) {
	  //	  cout<<"Pt: "<<jet1.pt<<"\t"<<"Eta: "<<jet1.eta<<"\t"<<"Phi: "<<jet1.phi<<"\t"<<"Energy: "<<jet1.energy<<endl;
	  //  cout<<"Pt: "<<jet2.pt<<"\t"<<"Eta: "<<jet2.eta<<"\t"<<"Phi: "<<jet2.phi<<"\t"<<"Energy: "<<jet2.energy<<endl;
	  //JetPairSelected.push_back(std::make_pair (jet1, jet2));
	  SignalJetsSet.insert(jet1);
	  SignalJetsSet.insert(jet2);
	}
      }
    }
    for (auto it : SignalJetsSet) {
      //      cout<<"Pt: "<<it.pt<<"\t"<<"Eta: "<<it.eta<<"\t"<<"Phi: "<<it.phi<<"\t"<<"Energy: "<<it.energy<<endl;
      SignalJets.push_back(it);
    }
    TLorentzVector sjp4;
    sjp4.SetPtEtaPhiE(0.,0.,0.,0.);
    if(SignalJets.size() >= 4) {
      for (size_t i = 0; i < 4; ++i){
	auto &j = SignalJets[i];
	sjp4 += LL4JMETUtil::getP4(j);
      }
      AnaUtil::fillHist1D("HMassSig", ljp4.M(), puevWt);
    }

    double H_TrMass = TMath::Sqrt(2*(l1p4+l2p4).Pt()*(mt.met)*(1-TMath::Cos(TVector2::Phi_mpi_pi((l1p4+l2p4).Phi()-mt.metphi))));
    AnaUtil::fillHist1D("HTrMass", H_TrMass, puevWt);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    //2...) Cut (1 + nTightJet >= 2)
    if (TightJets.size() < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, puevWt);
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
    for (auto &j: TightJets){
      TLorentzVector jp4 = LL4JMETUtil::getP4(j);
      hT += jp4.Pt();
    }
    AnaUtil::fillHist1D("hT_2", hT, puevWt);

    double sT = hT + lep1pt + lep2pt + mt.met;
    AnaUtil::fillHist1D("sT_2", sT, puevWt);


    TLorentzVector j1j2p4 = j1P4 + j2P4;
    TLorentzVector l1l2p4 = l1p4 + l2p4;
    double jlDPhi   = std::fabs(TVector2::Phi_mpi_pi(j1j2p4.Phi()-l1l2p4.Phi()));
    double j1l1DPhi = std::fabs(TVector2::Phi_mpi_pi(j1P4.Phi()-l1p4.Phi()));
    double j1l1DR   = j1P4.DeltaR(l1p4);
    AnaUtil::fillHist1D("JetLepDPhi_2", jlDPhi, puevWt);
    AnaUtil::fillHist1D("corrmet_2",  mt.met, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_2", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_2", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_2", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_2", nbJets(), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_2", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_2", j1l1DR, puevWt);

    //3...) Cut (2 + jet1Pt > 60)

    if (j1P4.Pt() < 60) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, puevWt);

    AnaUtil::fillHist1D("nTightJet_3", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_3", hT, puevWt);
    AnaUtil::fillHist1D("sT_3", sT, puevWt);
    AnaUtil::fillHist1D("jet2pt_3", TightJets[1].pt, puevWt);
    AnaUtil::fillHist1D("corrmet_3",  mt.met, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_3", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_3", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_3", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_3", nbJets(), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_3", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_3", j1l1DR, puevWt);


    //4...) Cut (3 + jet2Pt > 40)
    if (j2P4.Pt() < 40) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, puevWt);
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
    AnaUtil::fillHist1D("j1l1DPhi_4", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_4", j1l1DR, puevWt);



    //5...) Cut (4 + met > 60)
    if (mt.met < 60) continue;
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, puevWt);
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
    AnaUtil::fillHist1D("j1l1DPhi_5", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_5", j1l1DR, puevWt);




    //6...) Cut (5 + hT > 160)
    if (hT < 160) continue;
    AnaUtil::fillHist1D("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 10, puevWt);
    AnaUtil::fillHist1D("nTightJet_6", TightJets.size(), puevWt);
    //    AnaUtil::fillHist1D("hT_6", hT, puevWt);
    AnaUtil::fillHist1D("sT_6", sT, puevWt);
    AnaUtil::fillHist1D("jdR_6", j1j2DR, puevWt);
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
    AnaUtil::fillHist1D("j1l1DPhi_6", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_6", j1l1DR, puevWt);


    //7...) Cut(6 + LepSumPt > 80)
    if (lSumPt < 80) continue;
    AnaUtil::fillHist1D("evtCutFlow", 11);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 11, puevWt);
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
    AnaUtil::fillHist1D("lT_7", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_7", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_7", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_7", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_7", nbJets(), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_7", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_7", j1l1DR, puevWt);


    //8...) Cut (7 + j1j2DR > 1.5)
    if (j1j2DR > 1.5) continue;
    AnaUtil::fillHist1D("evtCutFlow", 12);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12, puevWt);
    AnaUtil::fillHist1D("nTightJet_8", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_8", hT, puevWt);
    AnaUtil::fillHist1D("sT_8", sT, puevWt);
    //    AnaUtil::fillHist1D("jdR_8", j1j2DR, puevWt);
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
    AnaUtil::fillHist1D("nbJets_8", nbJets(), puevWt);
    AnaUtil::fillHist1D("lepSumChr_8", (lep1c+lep2c), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_8", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_8", j1l1DR, puevWt);


    //9...) Cut (8 + sameChargedLeptons)
    if (lep1c*lep2c < 0.0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 13);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 13, puevWt);
    AnaUtil::fillHist1D("nTightJet_9", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_9", hT, puevWt);
    AnaUtil::fillHist1D("sT_9", sT, puevWt);
    AnaUtil::fillHist1D("jdR_9", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_9", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_9", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_9", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_9", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_9", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_9", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_9", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_9", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_9", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_9", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_9", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_9", nbJets(), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_9", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_9", j1l1DR, puevWt);


    //10...) Cut (9 + noZToLL)
    if (hasZcandidate(LepCandList_, puevWt)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 14);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 14, puevWt);
    AnaUtil::fillHist1D("nTightJet_10", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_10", hT, puevWt);
    AnaUtil::fillHist1D("sT_10", sT, puevWt);
    AnaUtil::fillHist1D("jdR_10", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_10", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_10", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_10", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_10", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_10", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_10", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_10", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_10", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_10", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_10", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_10", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_10", nbJets(), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_10", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_10", j1l1DR, puevWt);



    //11...) Cut (10 + 0.4 <lepDR < 1.5)
    if ( l1p4.DeltaR(l2p4) < 0.4 || l1p4.DeltaR(l2p4) > 1.5) continue; //lepDR
    AnaUtil::fillHist1D("evtCutFlow", 15);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 15, puevWt);
    AnaUtil::fillHist1D("nTightJet_11", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_11", hT, puevWt);
    AnaUtil::fillHist1D("sT_11", sT, puevWt);
    AnaUtil::fillHist1D("jdR_11", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_11", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_11", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_11", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_11", lep2pt, puevWt);
    //    AnaUtil::fillHist1D("lepDR_11", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_11", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_11", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_11", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_11", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_11", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_11", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_11", nbJets(), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_11", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_11", j1l1DR, puevWt);
    


    //12...) Cut (11 + lep1MetDPhi < 0.5)

    if (l1metDPhi > 0.7) continue; 
    AnaUtil::fillHist1D("evtCutFlow", 16);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 16, puevWt);
    AnaUtil::fillHist1D("nTightJet_12", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_12", hT, puevWt);
    AnaUtil::fillHist1D("sT_12", sT, puevWt);
    AnaUtil::fillHist1D("jdR_12", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_12", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_12", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_12", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_12", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_12", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_12", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_12", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_12", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_12", Mass_2l, puevWt);
    //AnaUtil::fillHist1D("lep1MetDPhi_12", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_12", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("nbJets_12", nbJets(), puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_12", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_12", j1l1DR, puevWt);



    //13...) Cut (12 + bVeto)

    if (nbJets() > 0) continue; //Signal has no b-jet
    AnaUtil::fillHist1D("evtCutFlow", 17);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 17, puevWt);
    AnaUtil::fillHist1D("nTightJet_13", TightJets.size(), puevWt);
    AnaUtil::fillHist1D("hT_13", hT, puevWt);
    AnaUtil::fillHist1D("sT_13", sT, puevWt);
    AnaUtil::fillHist1D("jdR_13", j1j2DR, puevWt);
    AnaUtil::fillHist1D("jdPhi_13", j1j2DPhi, puevWt);
    AnaUtil::fillHist1D("JetLepDPhi_13", jlDPhi, puevWt);
    AnaUtil::fillHist1D("lep1pt_13", lep1pt, puevWt);
    AnaUtil::fillHist1D("lep2pt_13", lep2pt, puevWt);
    AnaUtil::fillHist1D("lepDR_13", l1p4.DeltaR(l2p4), puevWt);
    AnaUtil::fillHist1D("lepDPhi_13", lepDPhi, puevWt);
    AnaUtil::fillHist1D("lSumPt_13", lSumPt, puevWt);
    AnaUtil::fillHist1D("lT_13", lT, puevWt);

    AnaUtil::fillHist1D("l1l2InvM_13", Mass_2l, puevWt);
    AnaUtil::fillHist1D("lep1MetDPhi_13", l1metDPhi, puevWt);
    AnaUtil::fillHist1D("lep2MetDPhi_13", l2metDPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DPhi_13", j1l1DPhi, puevWt);
    AnaUtil::fillHist1D("j1l1DR_13", j1l1DR, puevWt);



    //14...) Cut (13 + tauVeto)
    if (nTau > 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 18);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 18, puevWt);



    //dumpEvent(vz, false, true); //dump in detail

    if (isMC()){
      if (LepCandList_[0].flavour == 1 && LepCandList_[1].flavour == 1) AnaUtil::fillHist1D("diffFlvYield", 0, puevWt);
      else if ((LepCandList_[0].flavour == 1 && LepCandList_[1].flavour == 2)||(LepCandList_[0].flavour == 2 && LepCandList_[1].flavour == 1)) AnaUtil::fillHist1D("diffFlvYield", 1, puevWt);
      else if (LepCandList_[0].flavour == 2 && LepCandList_[1].flavour == 2) AnaUtil::fillHist1D("diffFlvYield", 2, puevWt);
    }
    


    //_________________________________________End of PreSelection_______________________________________________//    
      
    histf()->cd();
    histf()->cd("CUTAnalysis");
    


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
  //TCanvas *c = new TCanvas("c", "c", 600, 400);
  //wtDiff -> Draw();
  //c->Update();
  //c->SaveAs("WtDiff.png");
  // Analysis over
}

bool MultiLeptonCUTAna::SearchMinTwoUniqueJetPairs(const std::vector<std::pair<vhtm::Jet, vhtm::Jet>>& JetPair) {
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

bool MultiLeptonCUTAna::hasZcandidate(const std::vector<LeptonCandidate>& lepColl, double puWt){
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
	  if (lepInvM > 85. && lepInvM < 100.) hasZMass = true;
	}
      }
      if (hasZToLL||hasZMass) return true;
    }
  }
  return false;
}

void MultiLeptonCUTAna::endJob() {
  PhysicsObjSelector::endJob();
  
  histf()->cd();
  histf()->cd("CUTAnalysis");
  vector<string> evLabels {
    "0) Events processed: ",
      "1) GenFilter_lep: ",
      "2) GenFilter_jet: ",
      "3) Have good Vtx: ",
      "4) triggered: ",
      "5) nTightIsoLep >= 2:",
      "6) nTightJet >= 2: ",
      "7) jet1pt >= 60:",
      "8) jet2pt >= 40:",
      "9) met >= 60 GeV: ",
      "10)hT >= 160: ",
      "11)lepSumPt >= 80: ",
      "12)j1j2dR <= 1.5:",
      "13)SameChrLeptons: ",
      "14)No ZToLL: ",
      "15)0.4 <= lepDR <= 1.5: ",
      "16)abs of l1metphi <= 0.5: ",
      "17)no bJet: ",
      "18)no Tau: "
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
 
    //comment it out if run in a single job, otherwise lumiscaling would be wrong
   
    LL4JMETUtil::scaleHistogram("evtCutFlowWt", lumiFac);
    LL4JMETUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");  
  }
}

void MultiLeptonCUTAna::closeFiles() {
  AnaBase::closeFiles();
  // Take care of local stuff first                                                
  //  if (_mvaObj != nullptr) _mvaObj->close();
  // if (skimObj_ != nullptr) skimObj_->close();
  //if (syncDumpf_.is_open()) syncDumpf_.close();
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
bool MultiLeptonCUTAna::readJob(const string& jobFile, int& nFiles)
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
void MultiLeptonCUTAna::printJob(ostream& os) const
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
