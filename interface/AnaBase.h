#ifndef __AnaBase__hh
#define __AnaBase__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TVector3.h"

#include "PhysicsObjects.h"
#include "AnaUtil.h"

using uint = unsigned int;
using ulong = unsigned long;

// REQUIRED, most probably to solve circular dependence problem!!!
class TChain;
class TFile;
using std::ofstream;
using std::ifstream;

template <class T>
class PtComparator {
public:
  bool operator()(const T &a, const T &b) const {
    return a.pt > b.pt;
  }
};

template <class T>
class PtComparatorTL {
public:
  bool operator()(const T &a, const T &b) const {
    return a.Pt() > b.Pt();
  }
};

template <class T>
class PtComparatorLep {
public:
  bool operator()(const T &a, const T &b) const {
    return a.lPt > b.lPt;
  }
};

template <class T>
class MassComparator {
public:
  bool operator()(const T &a, const T &b) const {
    TLorentzVector l1,l2;
    l1.SetPtEtaPhiE(a.pt, a.eta, a.phi, a.energy);
    l2.SetPtEtaPhiE(b.pt, b.eta, b.phi, b.energy);
    return l1.M() > l2.M();
  }
};

typedef struct  
{
  bool verbose;
  bool usesbit;
  bool printselected;
} Options;

class AnaBase {
    
public:

  AnaBase();
  virtual ~AnaBase();
    
  virtual void eventLoop() = 0;  // the main analysis
  virtual bool beginJob();
  virtual void endJob();
  virtual bool readJob(const std::string& jobFile, int& nFiles);
  virtual void printJob(std::ostream& os=std::cout) const;

  virtual bool openFiles();
  virtual void closeFiles(); 
  virtual void closeHistFile(); 

  int setInputFile(const std::string& fname);
  bool branchFound(const std::string& b);
  int getEntries() const;
  void setData(int val);  
  int getRunNumber() const;
  void showEventNumber(std::ostream& os=std::cout, bool addNewline=true) const {
    const vhtm::Event& evt = eventList_->at(0);
    os << ">>> Run " << evt.run
       << " Lumis " << evt.lumis 
       << " Event " << evt.event;
    if (addNewline) os << std::endl;
  }
  bool isTriggered(bool check_prescale=true, bool verbose=false) const;
  void dumpTriggerPaths(std::ostream& os=std::cout, bool check_prescale=true) const;
  void dumpTriggerObjectInfo(const std::vector<vhtm::TriggerObject>& list, std::ostream& os=std::cout) const;
  double wtPileUp(int& nPU, bool verbose=false) const;
  bool readPileUpHist(bool verbose=false);
  bool matchTriggerPath(const std::vector<std::string>& v, const std::string& path) const;
  double matchTriggerObject(const std::vector<vhtm::TriggerObject>& trigObjList, 
                            const TLorentzVector& obj, 
                            const std::string& trigpath, 
                            int trig_skip, 
                            double maxPtDiff, 
                            int& trig_indx) const;
  
  void clearEvent();
  void enableBranches();
   int getEntry(int lflag) const;
  void setAddresses(); 

  void dumpGenInfo(std::ostream& os=std::cout) const;
  int getMotherId(const vhtm::GenParticle& gp, int& mmid) const;
  int getMotherIdForQ(const vhtm::GenParticle& gp, int& mmid) const;

  void findVtxInfo(std::vector<vhtm::Vertex>& list, Options& op, std::ostream& os=std::cout);
  void findTriggerObjectInfo(std::vector<vhtm::TriggerObject>& list);
  TVector3 findLeptonVtx(int index, bool& isGoodVtx);

  const std::vector<vhtm::Event>* eventColl() const {return eventList_;}
  const std::vector<vhtm::Vertex>* vertexColl() const {return vertexList_;}
  const std::vector<vhtm::GenEvent>* genEventColl() const {return genEventList_;}
  const std::vector<vhtm::Tau>* tauColl() const {return tauList_;}
  const std::vector<vhtm::Electron>* electronColl() const {return electronList_;}
  const std::vector<vhtm::Muon>* muonColl() const {return muonList_;}
  const std::vector<vhtm::Photon>* photonColl() const {return photonList_;}
  const std::vector<vhtm::PackedPFCandidate>* packedPFCandidateColl() const {return packedPFCandidateList_;}
  const std::vector<vhtm::Jet>* jetColl() const {return jetList_;}
  const std::vector<vhtm::MET>* metColl() const {return metList_;}
  const std::vector<vhtm::MET>* corrmetColl() const {return corrmetList_;}
  const std::vector<vhtm::MET>* puppimetColl() const {return puppimetList_;}
  const std::vector<vhtm::GenParticle>* genParticleColl() const {return genParticleList_;}
  const std::vector<vhtm::GenJet>* genJetColl() const {return genJetList_;}
  const std::vector<vhtm::GenMET>* genMetColl() const {return genMetList_;}
  const std::vector<vhtm::TriggerObject>* triggerObjColl() const {return triggerObjList_;}

  const std::vector<int>* l1physbits() const {return l1physbits_;}
  const std::vector<int>* l1techbits() const {return l1techbits_;}
  const std::vector<std::string>* hltpaths() const {return hltpaths_;}
  const std::vector<int>* hltresults() const {return hltresults_;}
  const std::vector<int>* hltprescales() const {return hltprescales_;}

  int nvertex() const {return vertexList_->size();}
  int nelectron() const {return electronList_->size();}
  int nmuon() const {return muonList_->size();}
  int nphoton() const {return photonList_->size();}
  int npackedPFCandidate() const {return packedPFCandidateList_->size();}
  int ntau() const {return tauList_->size();}
  int njet() const {return jetList_->size();}
  int nmet() const {return metList_->size();}
  int ncorrmet() const {return corrmetList_->size();}
  int npuppimet() const {return puppimetList_->size();}
  int ngenparticle() const {return genParticleList_->size();}
  int ntriggerobj() const {return triggerObjList_->size();}
  int ngenjet() const {return genJetList_->size();}
  int ngenmet() const {return genMetList_->size();}

  const std::unique_ptr<TChain>& chain() const {return chain_;}
  std::unique_ptr<TChain>& chain() {return chain_;}
  const std::unique_ptr<TFile>& histf() const {return histf_;}
  std::unique_ptr<TFile>& histf() {return histf_;}

  int nEvents() const {return nEvents_;}
  int firstEvent() const {return firstEvt_;}
  int lastEvent() const {return lastEvt_;}

  ofstream& fLog() {return fLog_;}
  ofstream& evLog() {return evLog_;}
  ofstream& selEvLog() {return selEvLog_;}

  const ofstream& fLog() const {return fLog_;}
  const ofstream& evLog() const {return evLog_;}
  const ofstream& selEvLog() const {return selEvLog_;}

  bool readGenInfo() const {return readGenInfo_;}
  bool isMC() const {return isMC_;}
  bool isSignal() const {return isSignal_;}
  int logOption() const {return logOption_;}
  bool useTrigger() const {return useTrigger_;}
  bool useLumiWt() const {return useLumiWt_;}
  double lumiWt(double evtWeightSum=-1, bool verbose=false) const;
  bool usePUWt() const {return usePUWt_;}
  const std::vector<std::string>& trigPathList() const {return trigPathList_;}
  bool useTrueNInt() const {return useTrueNInt_;}

  const std::map<std::string, double>& lumiWtMap() const {return AnaUtil::cutMap(hmap_, "lumiWtList");}
  const std::map<std::string, double>& vtxCutMap() const {return AnaUtil::cutMap(hmap_, "vtxCutList");}
  const std::map<std::string, double>& muonCutMap() const {return AnaUtil::cutMap(hmap_, "muonCutList");}
  const std::map<std::string, double>& photonCutMap() const {return AnaUtil::cutMap(hmap_, "photonCutList");}
  const std::map<std::string, double>& packedPFCandidateCutMap() const {return AnaUtil::cutMap(hmap_, "packedPFCandidateCutList");}
  const std::map<std::string, double>& electronCutMap() const {return AnaUtil::cutMap(hmap_, "electronCutList");}
  const std::map<std::string, double>& tauCutMap() const {return AnaUtil::cutMap(hmap_, "tauCutList");}
  const std::map<std::string, double>& jetCutMap() const {return AnaUtil::cutMap(hmap_, "jetCutList");}
  const std::map<std::string, double>& evselCutMap() const {return AnaUtil::cutMap(hmap_, "evselCutList");}

  const std::unordered_map<std::string, int>& eventIdMap() const {return eventIdMap_;}
  int bunchCrossing() const {return bunchCrossing_;}
  
private:
  std::unique_ptr<TChain> chain_;      // chain contains a list of root files containing the same tree
  std::unique_ptr<TFile> histf_;       // The output file with histograms

  // The tree branches
  std::vector<vhtm::Event>* eventList_;
  std::vector<vhtm::Vertex>* vertexList_;
  std::vector<vhtm::GenEvent>* genEventList_;
  std::vector<vhtm::Tau>* tauList_;
  std::vector<vhtm::Electron>* electronList_;
  std::vector<vhtm::Muon>* muonList_;
  std::vector<vhtm::Photon>* photonList_;
  std::vector<vhtm::PackedPFCandidate>* packedPFCandidateList_;
  std::vector<vhtm::Jet>* jetList_;
  std::vector<vhtm::MET>* metList_;
  std::vector<vhtm::MET>* corrmetList_;
  std::vector<vhtm::MET>* puppimetList_;
  std::vector<vhtm::GenParticle>* genParticleList_;
  std::vector<vhtm::GenJet>* genJetList_;
  std::vector<vhtm::GenMET>* genMetList_;
  std::vector<vhtm::TriggerObject>* triggerObjList_;

  std::vector<int>* l1physbits_;
  std::vector<int>* l1techbits_;
  std::vector<std::string>* hltpaths_;
  std::vector<int>* hltresults_;
  std::vector<int>* hltprescales_;

  std::vector<std::string> brList_;
  std::vector<double> puWtList_;

  int nEvents_;

  ofstream fLog_;   
  ofstream evLog_;   
  ofstream selEvLog_;   

  bool isMC_ {false};
  bool isSignal_ {false};
  bool readGenInfo_ {false};
  bool readTrigObject_ {false};
  bool readPFObject_ {false};
  std::vector<std::string> fileList_;
  int logOption_ {0};
  bool useTrigger_ {false};
  bool useLumiWt_ {false};
  bool usePUWt_ {false};
  std::string puHistFile_ {"./reweightFunctionFall11.root"};
  std::string puHistogram_ {"pileup"};
  bool useTrueNInt_ {true};
  int bunchCrossing_ {25};

  std::string histFile_ {"default.root"};
  std::string logFile_ {"default.out"};
  std::string evFile_ {"events.out"};
  std::string selEvFile_ {"selected_events.out"};
  int maxEvt_ {0};
  int firstEvt_ {-1};
  int lastEvt_ {-1};

  std::vector<std::string> trigPathList_;
  std::map<std::string, std::map<std::string, double>> hmap_;
  std::unordered_map<std::string, int> eventIdMap_;
};
#endif
