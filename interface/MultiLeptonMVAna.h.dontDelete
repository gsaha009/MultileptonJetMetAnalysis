#ifndef __MultiLeptonMVAna__hh
#define __MultiLeptonMVAna__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <bitset>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "PhysicsObjSelector.h"
#include "LeptonCandidate.h"
#include "AnaUtil.h"
#include "MVASkim.h"
#include "MVAnalysis.h"

struct MakeUnique{
  bool operator() (const vhtm::Jet& a, const vhtm::Jet& b) {
    //bool hasSameDR {false};
    //bool hasSamePt {false};
    TLorentzVector a4, b4;
    a4.SetPtEtaPhiE(a.pt, a.eta, a.phi, a.energy);
    b4.SetPtEtaPhiE(b.pt, b.eta, b.phi, b.energy);
    //return (!(std::fabs(a4.Pt() - b4.Pt()) < 1.0e-08 && a4.DeltaR(b4) < 1.0e-06));
    //return (a4.Pt() < b4.Pt() && a4.DeltaR(b4) < 1.0e-06);
    //    if (std::fabs(a4.Pt() - b4.Pt()) > 1.0e-05) return true;
    //return false;
    return a4.Pt() > b4.Pt();
  }
};


class MultiLeptonMVAna: public PhysicsObjSelector {
    
public:
  enum class EventType {
    unkwn = -1, mmem = 0, eeem
  };
  MultiLeptonMVAna();
  virtual ~MultiLeptonMVAna();
  
  virtual void eventLoop() final;  // the main analysis 
  virtual bool beginJob() override;
  virtual void endJob() override;
  virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  virtual void printJob(std::ostream& os=std::cout) const override;
  virtual void bookHistograms() override;
  virtual void closeFiles() override;
    
  void clearLists();

  bool SearchMinTwoUniqueJetPairs (const std::vector<std::pair<vhtm::Jet, vhtm::Jet>>& JetPair) ;
  bool hasZcandidate (const std::vector<LeptonCandidate>& lepColl, double puWt);
  bool _createMVATree {false};
  bool _readMVA {false};
  std::string _mvaInputFile {""};  
  std::string _MVAnetwork {""};
  std::string _MVAxmlFile {""};
  std::unique_ptr<MVAnalysis> _mvaObj {nullptr};
  std::unique_ptr<MVASkim> skimObj_ {nullptr};
  
  
private:
  std::vector<vhtm::Vertex> vtxList_;
  std::vector<LeptonCandidate> LepCandList_;
  std::vector<vhtm::Jet>JetSelected;
  std::set<vhtm::Jet, MakeUnique> SignalJetsSet;  
  //  std::set<struct MakeUnique> SignalJetsSet;  
  std::vector<vhtm::Jet>SignalJets;
  std::vector<std::pair<vhtm::Jet, vhtm::Jet>>JetPairSelected;

  double evtWeightSum_ {0};
  bool dumpGenInfo_ {false};
  bool useEventList_ {false};
  bool skipDuplicate_ {false};
  bool selectPM_ {false};
  int nMEPartons_ {-1};
  ofstream syncDumpf_;
  std::string dumpFilename_ {"syncDumpFile.txt"};
  int dumpEventCount_ {0};
  std::vector<std::string> eventFilelist_;

  std::unordered_map<std::string, int> eventIdStore_;
};
#endif
