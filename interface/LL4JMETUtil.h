#ifndef __LL4JMETUtil__hh
#define __LL4JMETUtil__hh

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"

#include "PhysicsObjects.h"
#include "AnaBase.h"
//#include "ZCandidate.h"

namespace LL4JMETUtil {
  const double MZnominal = 91.1876;

  enum ZType {
    mumu = 0, ee, 
  };
  
  enum llType {
    mue = 0, emu
  };
  enum lTauType {
    mutau = 0, etau
  };
  // must be defined inside the header to be effective
  template <class T> 
  TLorentzVector getP4(const T& obj) {
    TLorentzVector lv;
    lv.SetPtEtaPhiE(obj.pt, obj.eta, obj.phi, obj.energy);
    return lv;
  }
  double getEleRhoEffectiveArea(double eta);
  double getEleRhoEffectiveArea03(double etax);

  double computeMuonReliso(const vhtm::Muon& mu, 
			   const std::vector<vhtm::PackedPFCandidate>& fsrCandList, 
			   double vetoCone=0.01, 
			   double isoCone=0.4, 
			   bool verbose=false);
  double computeElectronReliso(const vhtm::Electron& ele, 
			       const std::vector<vhtm::PackedPFCandidate>& fsrCandList, 
			       double eventRho, 
			       double vetoCone=0.08, 
			       double isoCone=0.4, 
			       bool verbose=false);

  double pfiso(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum=0.0);
  double pfiso(const vhtm::Muon& mu, double fsrPhotonEtSum=0.0);
  double pfiso(const vhtm::PackedPFCandidate& cand);
  double pfiso03(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum=0.0);
  double pfiso03(const vhtm::Muon& mu, double fsrPhotonEtSum=0.0);

  bool jetpuMVAid(const vhtm::Jet& jet);
  void printP4(const TLorentzVector& lv, const std::string& tag, std::ostream& os=std::cout);
  bool electronBDT(const vhtm::Electron& electron);

  void scaleHistogram(const std::string& hname, double fac);
  void showEfficiency(const std::string& hname, 
		      const std::vector<std::string>& slist, 
		      const std::string& header, 
		      const std::string& tag="Events", std::ostream& os=std::cout);
  void showCount(const std::string& hname, 
		 const std::vector<std::string>& slist, 
		 const std::string& tag, int prec=0, std::ostream& os=std::cout);
  bool isLooseJet(const vhtm::Jet& jet);
  bool isTightJet(const vhtm::Jet& jet);
  bool fsrPhotonP4(const std::vector<vhtm::PackedPFCandidate>& fsrList, TLorentzVector& p4);
}
#endif
