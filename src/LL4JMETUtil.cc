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

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AnaUtil.h"
#include "LL4JMETUtil.h"

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

namespace LL4JMETUtil {
  double getEleRhoEffectiveArea(double etax) {
    double area = 0.0;
    double eta = std::fabs(etax);
    if      (eta >= 0.0   && 1.0)   area = 0.1752;
    else if (eta >= 1.0   && 1.479) area = 0.1862;
    else if (eta >= 1.479 && 2.0)   area = 0.1411;
    else if (eta >= 2.0   && 2.2)   area = 0.1534;
    else if (eta >= 2.2   && 2.3)   area = 0.1903;
    else if (eta >= 2.3   && 2.4)   area = 0.2243;
    else if (eta >= 2.4   && 5.0)   area = 0.2687;
    return area;
  }
  double getEleRhoEffectiveArea03(double etax) {
    double area = 0.0;
    double eta = std::fabs(etax);
    if (eta >= 0.0 && 1.0) area = 0.1752;
    else if (eta >= 1.0 && 1.479) area = 0.1862;
    else if (eta >= 1.479 && 2.0) area = 0.1411;
    else if (eta >= 2.0 && 2.2) area = 0.1534;
    else if (eta >= 2.2 && 2.3) area = 0.1903;
    else if (eta >= 2.3 && 2.4) area = 0.2243;
    else if (eta >= 2.4 && 5.0) area = 0.2687;
    return area;
  }
  double pfiso(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum) {
    return (ele.chargedHadronIso + std::max(0., ele.neutralHadronIso + ele.photonIso - fsrPhotonEtSum
  					    - getEleRhoEffectiveArea(std::fabs(ele.eta)) * eventRho));
  }
  double pfiso(const vhtm::Muon& mu, double fsrPhotonEtSum) {
    return (mu.sumChargedHadronPt + std::max(0., mu.sumNeutralHadronEt + mu.sumPhotonEt - fsrPhotonEtSum - 0.5 * mu.sumPUPt));
  }
  double pfiso(const vhtm::PackedPFCandidate& cand) {
    return (cand.isolationMap.at("c30").at(0) + 
	    cand.isolationMap.at("c30").at(2) + 
	    cand.isolationMap.at("c30").at(3) + 
	    cand.isolationMap.at("c30").at(4));
  }
  // note EG pog now uses scETA
  double pfiso03(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum) {
    return (ele.sumChargedHadronPt + 
            std::max(0., ele.sumNeutralHadronEt + ele.sumPhotonEt - fsrPhotonEtSum
	 	       - getEleRhoEffectiveArea03(std::fabs(ele.scEta)) * eventRho));
  }

  double pfiso03(const vhtm::Muon& mu, double fsrPhotonEtSum) {
    return mu.pfChargedHadIsoR03 + std::max(0.0, mu.pfNeutralHadIsoR03 + mu.pfPhotonIso03 
					    - fsrPhotonEtSum - 0.5 * mu.sumPUPt03);
  }
  // Isolation in new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone of all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons 
  double computeElectronReliso(const vhtm::Electron& ele, 
			       const std::vector<vhtm::PackedPFCandidate>& fsrCandList, 
			       double eventRho, 
			       double vetoCone, 
			       double isoCone, 
			       bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4(getP4(ele));
    for (const auto& v: fsrCandList) {
      TLorentzVector fsrP4(getP4(v));
      double dR = lP4.DeltaR(fsrP4);
      if ((std::fabs(ele.scEta) < 1.479 || dR > vetoCone) && dR < isoCone)
  	phoEtSum += fsrP4.Et();
    }
    double iso = (std::fabs(isoCone - 0.3) < 1.0e-10) 
      ? pfiso03(ele, eventRho, phoEtSum) 
      : pfiso(ele, eventRho, phoEtSum);

    if (verbose) {
      cout << "electron isolation: " << endl;
      cout << "      iso    lepPt  chHadPt neuHadEt photonEt    fsrEt  effArea eventRho" << endl;
      cout << setprecision(3) 
           << setw(9) << iso 
           << setw(9) << lP4.Pt()
  	   << setw(9) << ele.chargedHadronIso
  	   << setw(9) << ele.neutralHadronIso 
  	   << setw(9) << ele.photonIso
  	   << setw(9) << phoEtSum
  	   << setw(9) << getEleRhoEffectiveArea(std::fabs(lP4.Eta()))
           << setw(9) << eventRho
  	   << endl;
    }
    return iso/lP4.Pt();
  }
  // Isolation with new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // dR > 0.01 for muons
  double computeMuonReliso(const vhtm::Muon& mu, 
			   const std::vector<vhtm::PackedPFCandidate>& fsrCandList, 
			   double vetoCone, 
			   double isoCone, 
			   bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4(getP4(mu));
    for (const auto& v: fsrCandList) {
      TLorentzVector fsrP4(getP4(v));
      double dR = lP4.DeltaR(fsrP4);
      if (dR > vetoCone && dR < isoCone)
  	phoEtSum += fsrP4.Et();
    }
    double iso = (std::fabs(isoCone - 0.3) < 1.0e-10) 
      ? pfiso03(mu, phoEtSum) 
      : pfiso(mu, phoEtSum);

    if (verbose) {
      cout << "muon isolation: " << endl;
      cout << "      iso    lepPt  chHadPt neuHadEt photonEt    fsrEt     PUPt" << endl;
      cout << setprecision(3) 
           << setw(9) << iso 
           << setw(9) << lP4.Pt()
  	   << setw(9) << mu.sumChargedHadronPt
  	   << setw(9) << mu.sumNeutralHadronEt
  	   << setw(9) << mu.sumPhotonEt
  	   << setw(9) << phoEtSum
  	   << setw(9) << mu.sumPUPt
  	   << endl;
    }
    
    return iso/lP4.Pt();
  }  
  bool jetpuMVAid(const vhtm::Jet& jet) {
    float jpumva = jet.jpumva;
    double pt = jet.pt;
    double eta = std::fabs(jet.eta);
    bool passPU = true;
    if (pt > 20) {
      if (eta > 3.) {
  	if (jpumva <= -0.45) passPU = false;
      }
      else if (eta > 2.75) {
  	if (jpumva <= -0.55) passPU = false;
      }
      else if (eta > 2.5) {
  	if (jpumva <= -0.6) passPU = false;
      }
      else if (jpumva <= -0.63) passPU = false;
    }
    else {
      if (eta > 3.) {
  	if (jpumva <= -0.95) passPU = false;
      }
      else if (eta > 2.75) {
  	if (jpumva <= -0.94) passPU = false;
      }
      else if (eta > 2.5) {
  	if (jpumva <= -0.96) passPU = false;
      }
      else if (jpumva <= -0.95) passPU = false;
    }
    return passPU;
  }
  bool isLooseJet(const vhtm::Jet& jet) {
    bool centralCut = (std::fabs(jet.eta) <= 2.4) 
      ? (jet.chargedHadronEnergyFraction > 0 && 
  	 jet.chargedMultiplicity > 0 && 
  	 jet.chargedEmEnergyFraction < 0.99)
      : true;
    
    return (jet.neutralHadronEnergyFraction < 0.99 && 
  	    jet.neutralEmEnergyFraction < 0.99 &&
  	    (jet.chargedMultiplicity + jet.neutralMultiplicity) > 1 &&
  	    jet.muonEnergyFraction < 0.8 &&
  	    centralCut);
  }
  bool isTightJet(const vhtm::Jet& jet) {
    bool centralCut = (std::fabs(jet.eta) <= 2.4)
      ? (jet.chargedHadronEnergyFraction > 0 && 
  	 jet.chargedMultiplicity > 0 && 
  	 jet.chargedEmEnergyFraction < 0.9)
      : true;
    
    return (jet.neutralHadronEnergyFraction < 0.9 && 
  	    jet.neutralEmEnergyFraction < 0.9 && 
  	    (jet.chargedMultiplicity + jet.neutralMultiplicity) > 1 &&
  	    jet.muonEnergyFraction < 0.8 && 
  	    centralCut);
  }
  bool electronBDT(const vhtm::Electron& electron) {
    double scEta = std::fabs(electron.scEta);
    double elePt = electron.pt;
    double BDT = electron.BDT;

    bool isBDT = (elePt <= 10 && ((scEta < 0.8                    && BDT > -0.211)   ||
  				 ((scEta >= 0.8 && scEta < 1.479) && BDT > -0.396)   ||
  				  (scEta >= 1.479                 && BDT > -0.215))) ||
                 (elePt >  10 && ((scEta < 0.8                    && BDT > -0.870)   ||
  				 ((scEta >= 0.8 && scEta < 1.479) && BDT > -0.838)   ||
  				  (scEta >= 1.479                 && BDT > -0.763)));
    return isBDT;
  }
  
  void showEfficiency(const string& hname, 
		      const std::vector<std::string>& slist, 
		      const string& header, 
		      const string& tag, 
		      std::ostream& os) 
  {
    os << ">>> " << header << " Efficiency" << endl;
    TH1 *h = AnaUtil::getHist1D(hname);
    if (h != nullptr) {
      os << setw(64) << "CutFlow"
  	   << setw(13) << tag
  	   << setw(10) << "AbsEff"
  	   << setw(10) << "RelEff"
  	   << endl;
      os.precision(3);
      int nbins = h->GetNbinsX();
      for (int i = 1; i <= nbins; ++i) {
        double cont  = h->GetBinContent(1);
        double conti = h->GetBinContent(i);
        double contj = h->GetBinContent(i-1);
  	os << setw(64) << slist[i-1]
  	   << setprecision(0) 
	   << setw(13) << conti
	   << setprecision(5) 
	   << setw(10) << ((conti > 0) ? conti/cont : 0.0)
	   << setw(10) << ( i == 1 ? 1.0 :(contj > 0) ? conti/contj : 0.0)
	   << endl;
      }
    }
  }
  void scaleHistogram(const string& hname, double fac) {
    TH1 *h = AnaUtil::getHist1D(hname);
    if (h != nullptr) {
      int nbins = h->GetNbinsX();
      for (int i = 1; i <= nbins; ++i) {
        double cont = h->GetBinContent(i) * fac;
        double err  = h->GetBinError(i) * fac;
        h->SetBinContent(i, cont);
        h->SetBinError(i, err);
      }
    }
  }
  void showCount(const string& hname, 
		 const std::vector<std::string>& slist, 
		 const std::string& tag, 
		 int prec, 
		 std::ostream& os) 
  {
    os << ">>> " << tag << endl;
    TH1 *h = AnaUtil::getHist1D(hname);
    if (h == nullptr) return;

    os.precision(prec);
    os << setw(16) << " Total"
       << setw(10) << h->Integral()
       << endl;
    int nbins = h->GetNbinsX();
    os.precision(prec);
    for (int i = 1; i <= nbins; ++i)
      os << setw(16) << slist[i-1]
	 << setw(10) << h->GetBinContent(i)
	 << endl;
  }
  void printP4(const TLorentzVector& lv, const string& tag, std::ostream& os) {
    os << setprecision(3);
    os << tag << " = (" 
       << setw(7) << lv.Pt()  << "," 
       << setw(7) << lv.Eta() << "," 
       << setw(7) << lv.Phi() << "," 
       << setw(7) << lv.Energy() << ")" 
       << endl;
  }
  bool fsrPhotonP4(const vector<vhtm::PackedPFCandidate>& fsrList, TLorentzVector& p4) {
    bool hasFsr = false;
    p4.SetPtEtaPhiE(0., 0., 0., 0);  
    if (!fsrList.empty()) {
      const auto& cand = fsrList[0]; // only the first one is important
      p4.SetPtEtaPhiE(cand.pt, cand.eta, cand.phi, cand.energy);
      hasFsr = true;
    }
    return hasFsr;
  }
}
