#ifndef __MultiLeptonAnalysis__hh
#define __MultiLeptonAnalysis__hh

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
#include "ZCandidate.h"
#include "AnaUtil.h"

class GenAnaBase;

class MultiLeptonAnalysis: public PhysicsObjSelector {
    
public:
  enum class EventType {
    unkwn = -1, mmem = 0, eeem
  };
  MultiLeptonAnalysis();
  virtual ~MultiLeptonAnalysis();
  
  virtual void eventLoop() final;  // the main analysis 
  virtual bool beginJob() override;
  virtual void endJob() override;
  virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  virtual void printJob(std::ostream& os=std::cout) const override;
  virtual void bookHistograms() override;
    
  void clearLists();

  // find an opposite signed lepton-lepton pair with different flavors
  template <typename T1, typename T2>
    std::bitset<6> leptonPairSelector(const std::vector<std::pair<T1, std::vector<vhtm::PackedPFCandidate>>>& lep1PhotonPairList, 
				      const std::vector<std::pair<T2, std::vector<vhtm::PackedPFCandidate>>>& lep2PhotonPairList,  
				      const ZCandidate& eventZ, 
				      std::vector<ZCandidate>& leptonPairList,
				      double modelMass=-1.0, bool verbose=false) 
  {
    int flags [] {0,0,0,0,0,0};
#if 0
    setiosflags(std::ios::fixed);
    std::cout << std::setprecision(6);
#endif
    if (verbose)
      std::cout << "==> #lepton (1): " << lep1PhotonPairList.size()
		<< ", #lepton (2): " << lep2PhotonPairList.size()
		<< std::endl;
    for (size_t i = 0; i < lep1PhotonPairList.size(); ++i) {
      const auto& ielem = lep1PhotonPairList[i];
      
      TLorentzVector lep1P4(HZZ4lUtil::getP4(ielem.first));
      if (AnaUtil::sameObject(lep1P4, eventZ.l1P4) || AnaUtil::sameObject(lep1P4, eventZ.l2P4)) {
        if (verbose)
	std::cout << "=> Lepton i:" << i << " is a part of the Z! (bit: 0)" 
	          << std::endl;
#if 0
        std::cout << "=> same Object" << std::endl;
	std::cout << "         Pt        Eta        Phi" << std::endl;
	std::cout << std::setw(11) << lep1P4.Pt() 
		  << std::setw(11) << lep1P4.Eta() 
		  << std::setw(11) << lep1P4.Phi()
		  << std::endl; 
	std::cout << std::setw(11) << eventZ.l1P4.Pt() 
		  << std::setw(11) << eventZ.l1P4.Eta() 
		  << std::setw(11) << eventZ.l1P4.Phi()
		  << std::endl; 
	std::cout << std::setw(11) << eventZ.l2P4.Pt() 
		  << std::setw(11) << eventZ.l2P4.Eta() 
		  << std::setw(11) << eventZ.l2P4.Phi()
		  << std::endl; 
#endif
	++flags[0];
        continue;
      }
      
      // Angular separation of the lepton with the leptons from Z
      if (lep1P4.DeltaR(eventZ.l1P4) < 0.1 || lep1P4.DeltaR(eventZ.l2P4) < 0.1) {
        if (verbose)
	std::cout << "=> Lepton " << i << " close to a lepton from Z!" 
		  << " dR: (1: " << lep1P4.DeltaR(eventZ.l1P4) << ", 2: " 
		  << lep1P4.DeltaR(eventZ.l2P4) << ") (bit: 1)"
	          << std::endl;
	++flags[1];
  	continue;
      }

      TLorentzVector lep1FsrP4;
      bool lep1HasFsr = HZZ4lUtil::fsrPhotonP4(ielem.second, lep1FsrP4); 
      
      for (size_t j = 0; j < lep2PhotonPairList.size(); ++j) {
	const auto& jelem = lep2PhotonPairList[j];
	
        // Distinct object
	TLorentzVector lep2P4(HZZ4lUtil::getP4(jelem.first));
	if (AnaUtil::sameObject(lep2P4, eventZ.l1P4) || AnaUtil::sameObject(lep2P4, eventZ.l2P4)) {
	  if (verbose)
	    std::cout << "=> Lepton j:" << j << " is a part of the Z! (bit: 2)" 
		      << std::endl;
	  ++flags[2];
	  continue;
        }

	// Angular separation 
	if (lep2P4.DeltaR(eventZ.l1P4) < 0.1 || lep2P4.DeltaR(eventZ.l2P4) < 0.1) {
          if (verbose)
	  std::cout << "=> Lepton " << j << " close to a lepton from Z!" 
		    << " dR: (1: " << lep2P4.DeltaR(eventZ.l1P4) << ", 2: " 
		    << lep2P4.DeltaR(eventZ.l2P4) << ") (bit: 3)"
		    << std::endl;
	  ++flags[3];
	  continue;
	}

        // Require oppositely charged leptons
	if (ielem.first.charge * jelem.first.charge > 0) {
          if (verbose)
	  std::cout << "=> Same charge leptons (qi*qj): " << (ielem.first.charge * jelem.first.charge)
		    << " (bit: 4)"
		    << std::endl;
	  ++flags[4];
	  continue;
        }
		
	// Angular separation between the two leptons being considered (move to job card)
	if (lep1P4.DeltaR(lep2P4) < AnaUtil::cutValue(evselCutMap(), "minDRLP")) {
          if (verbose)
	  std::cout << "=> Leptons are too close! i: " << i << ", j: " << j 
		    << " dR: " << lep1P4.DeltaR(lep2P4)
		    << " (bit: 5)"
		    << std::endl;
	  ++flags[5];
	  continue;
	}
	
	TLorentzVector lep2FsrP4;
	bool lep2HasFsr = HZZ4lUtil::fsrPhotonP4(jelem.second, lep2FsrP4); 

	ZCandidate llp;
	if (typeid(ielem.first) == typeid(vhtm::Muon)) 
	  llp.flavour = HZZ4lUtil::llType::mue;
	else if (typeid(ielem.first) == typeid(vhtm::Electron))
	  llp.flavour = HZZ4lUtil::llType::emu;
	else 
	  llp.flavour = -1;
	
	llp.l1Index = i;
	llp.l1P4 = lep1P4;
	llp.l1Charge = ielem.first.charge;
	llp.l1FsrP4 = lep1FsrP4;
	llp.l1RP4 = lep1P4 + lep1FsrP4;
	
	llp.l2Index = j;
	llp.l2P4 = lep2P4;
	llp.l2Charge = jelem.first.charge;
	llp.l2FsrP4 = lep2FsrP4;
	llp.l2RP4 = lep2P4 + lep2FsrP4;
	
	int whichLep;
	if (lep1HasFsr && lep2HasFsr) whichLep = 3;
	else if (lep1HasFsr)          whichLep = 1;
	else if (lep2HasFsr)          whichLep = 2;
	else                 	      whichLep = 0;
	llp.fsrWithLep = whichLep;
	llp.fsrPhoP4 = lep1FsrP4 + lep2FsrP4;
	
        TLorentzVector p4 = llp.l1RP4 + llp.l2RP4;
        llp.p4 = p4;
	llp.mass = p4.M();
	llp.massDiff = (modelMass > 0) ? std::fabs(llp.mass - modelMass) : -999;
        llp.dEtall = llp.l1RP4.Eta() - llp.l2RP4.Eta();
        llp.dPhill = TVector2::Phi_mpi_pi(llp.l1RP4.Phi() - llp.l2RP4.Phi());
        llp.dRll = llp.l1RP4.DeltaR(llp.l2RP4);	

	leptonPairList.push_back(llp);
      }
    } 
#if 0
    resetiosflags(std::ios::fixed);
#endif
    std::bitset<6> result;
    int i = 0;
    for (auto v: flags) {
      if (verbose) std::cout << "i: " << i << ", v: " << v << std::endl;
      if (v > 0) result.set(i, 1); 
      i++;
    }
    return result;
  }
  int findExtraLeptons(const ZCandidate& Z1, const ZCandidate& Z2);
  int findEventCategory(int nleptons, const std::vector<vhtm::Jet>& jetList, int nbjets,
			const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, bool verbose=true);

  EventType setEventType(const ZCandidate& ZCand, const ZCandidate& leptonPairCand);
  static bool hasJetPair(const std::vector<vhtm::Jet>& jetList);
  static bool dmComparator(const ZCandidate& a, const ZCandidate& b) {
    return (a.massDiff < b.massDiff);
  }
  static bool massComparator(const ZCandidate& a, const ZCandidate& b) {
    return (a.mass > b.mass);
  }

private:
  std::vector<vhtm::Vertex> vtxList_;
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
