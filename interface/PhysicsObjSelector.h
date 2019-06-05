#ifndef __PhysicsObjSelector__hh
#define __PhysicsObjSelector__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "AnaBase.h"
#include "LL4JMETUtil.h"
//#include "ZCandidate.h"
#include "LeptonCandidate.h"

class PhysicsObjSelector: public AnaBase {
 public:
  PhysicsObjSelector();
  virtual ~PhysicsObjSelector() {
  }
  virtual void eventLoop() = 0;  // the main analysis
  virtual bool beginJob() override;
  virtual void endJob() override;
  virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  virtual void bookHistograms();
  void objectEfficiency();

  const std::vector<vhtm::PackedPFCandidate>& getFSRPhotonList() const {return fsrPhotonList_;}
  const std::vector<vhtm::Muon>& getPreSIPLooseMuList() const {return preSIPLooseMuList_;}
  const std::vector<vhtm::Electron>& getPreSIPLooseEleList() const {return preSIPLooseEleList_;}

  const std::vector<vhtm::Tau>& getTauList() const {return tauList_;}
  const std::vector<vhtm::Tau>& getIsoTauList() const {return isoTauList_;}

  const std::vector<vhtm::Muon>& getLooseMuList() const {return looseMuList_;}
  const std::vector<vhtm::Muon>& getTightMuList() const {return tightMuList_;}
  const std::vector<vhtm::Muon>& getTightIsoMuList() const {return tightIsoMuList_;}

  const std::vector<vhtm::Electron>& getLooseEleList() const {return looseEleList_;}
  const std::vector<vhtm::Electron>& getTightEleList() const {return tightEleList_;}
  const std::vector<vhtm::Electron>& getTightIsoEleList() const {return tightIsoEleList_;}

  const std::vector<vhtm::Jet>& getLeptonCleanedLooseJetList() const {return looseJetList_;}
  const std::vector<vhtm::Jet>& getLeptonCleanedTightJetList() const {return tightJetList_;}

  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& 
    getLooseElePhotonPairList() const {return looseElePhotonPairList_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& 
    getTightElePhotonPairList() const {return tightElePhotonPairList_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& 
    getTightIsoElePhotonPairList() const {return tightIsoElePhotonPairList_;}

  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& 
    getLooseMuPhotonPairList() const {return looseMuPhotonPairList_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& 
    getTightMuPhotonPairList() const {return tightMuPhotonPairList_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& 
    getTightIsoMuPhotonPairList() const {return tightIsoMuPhotonPairList_;}

  void setTightElePhotonPairList(const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& list) 
  {tightElePhotonPairList_ = list;}

  int nLooseJets() const {return looseJetList_.size();}
  int nTightJets() const {return tightJetList_.size();}
  int nbJets() const {return nbJets_;} 
  const std::vector<vhtm::Jet>& getLooseJetList() const {return looseJetList_;}
  const std::vector<vhtm::Jet>& getTightJetList() const {return tightJetList_;}  
  int getTrueLooseLeptons(double dRCut=0.4) const;
  std::tuple<int,int,int,int> findIsoLeptons(double dRCut=0.4) const;
  
  void findObjects(double dz,double wt=1);
  //  void findObjects(double wt=1);
  void muonSelector(double wt=1);
  void tauSelector(double vz,double wt=1);
  //void tauSelector(double wt=1);
  void electronSelector(double wt=1);
  void photonSelector(double wt=1);
  void isoLeptonSelector(double muIso=0.35, double eleIso=0.35);

  void jetSelector(double wt=1);
  bool jetLeptonCleaning(const vhtm::Jet& jet, double dR) const;
  void setEventGridRho(double evRho) {fGridRhoFastjetAll_ = evRho;}
  double getEventGridRho() const {return fGridRhoFastjetAll_;}
  void clear();

  bool passedSuperClusterVeto(const vhtm::PackedPFCandidate& pfcand, bool verbose=false) const;
  bool passedSuperClusterVetobyReference(const vhtm::PackedPFCandidate& pfcand, bool verbose=false) const;
  double findClosestLepton(const vhtm::PackedPFCandidate& photon, int& muindx, int& elindx) const;
  void leptonCrossCleaning();
  bool crossCleaned(const vhtm::Electron& electron) const;
  void addLeptonIsolation(std::vector<LeptonCandidate>& lepCandList,
                          const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& elePhotonPairList,
                          const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>>& muPhotonPairList);

  void dumpEvent(double vz, bool dumpGen=false, bool showEvent=false, std::ostream& os=std::cout) const;

  // Put all the lepton candidates into a vector                                                                                                                                                            
  template <typename T>
    void packLeptons(const std::vector<std::pair<T, std::vector<vhtm::PackedPFCandidate>>>& lepPhotonPairList, std::vector<LeptonCandidate>& candList) {
    for (unsigned int i = 0; i < lepPhotonPairList.size(); ++i) {
      const auto& ip = lepPhotonPairList[i];
      TLorentzVector lepP4(LL4JMETUtil::getP4(ip.first));
      TLorentzVector lepFsrP4;
      //bool lepHasFsr = LL4JMETUtil::fsrPhotonP4(ip.second, lepFsrP4);
      LeptonCandidate lc;
      lc.lIndex = i;
      lc.lCharge = ip.first.charge;
      lc.lPt = ip.first.pt;
      lc.lEta = ip.first.eta;
      if (typeid(ip.first) == typeid(vhtm::Muon))          lc.flavour = 1;
      else if (typeid(ip.first) == typeid(vhtm::Electron)) lc.flavour = 2;
      else                                                 lc.flavour =-1;
      lc.lP4 = lepP4;
      lc.lFsrP4 = lepFsrP4;
      lc.lRP4 = lepP4 + lepFsrP4;
      candList.push_back(lc);
    }
  }
 
 private:
  bool dumpEvent_;
  std::vector<vhtm::Tau> tauList_, 
                         isoTauList_;
  std::vector<vhtm::Muon> preSIPLooseMuList_, 
                          looseMuList_, 
                          tightMuList_, 
                          tightIsoMuList_;
  std::vector<vhtm::Electron> preSIPLooseEleList_, 
                              looseEleList_, 
                              tightEleList_, 
                              tightIsoEleList_;
  std::vector<vhtm::PackedPFCandidate> fsrPhotonList_;
  std::vector<vhtm::Jet> looseJetList_, tightJetList_;

  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>> looseElePhotonPairList_; 
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>> tightElePhotonPairList_; 
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>> tightIsoElePhotonPairList_; 

  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>> looseMuPhotonPairList_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>> tightMuPhotonPairList_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate>>> tightIsoMuPhotonPairList_; 

  double fGridRhoFastjetAll_;
  bool searchedEle_ {false}, 
       searchedMu_ {false}, 
       searchedPhoton_ {false};
  int nbJets_ {0};
};
#endif
