#ifndef __MultiLeptonJetMetAnalysis__hh
#define __MultiLeptonJetMetAnalysis__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
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
//#include "ZCandidate.h"
#include "AnaUtil.h"

class GenAnaBase;

class MultiLeptonJetMetAnalysis: public PhysicsObjSelector {
    
public:
  enum class EventType {
    unkwn = -1, mmem = 0, eeem
  };
  MultiLeptonJetMetAnalysis();
  virtual ~MultiLeptonJetMetAnalysis();
  
  virtual void eventLoop() final;  // the main analysis 
  virtual bool beginJob() override;
  virtual void endJob() override;
  virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  virtual void printJob(std::ostream& os=std::cout) const override;
  virtual void bookHistograms() override;
    
  void clearLists();

  static bool hasJetPair(const std::vector<vhtm::Jet>& jetList);
  
private:
  std::vector<vhtm::Vertex> vtxList_;
  std::vector<LeptonCandidate> LepCandList_;
  std::vector<LeptonCandidate> finalLepCandList_;
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

  std::unique_ptr<GenAnaBase> genAna_;

#ifdef SKIP_DUPLICATE_ZMASS
  std::unordered_map<std::string, int> eventMap_;
#endif
#ifdef SKIP_DUPLICATE_ALL
  std::unordered_map<std::string, int> eventIdStore_;
#endif
};
#endif
