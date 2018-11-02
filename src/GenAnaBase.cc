#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <utility> 

#include "AnaUtil.h"
#include "LL4JMETUtil.h"
#include "GenAnaBase.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::setw;
using std::string;

GenAnaBase::~GenAnaBase() {
}
void GenAnaBase::dumpEvent(std::ostream& os) const {
  int ngenp = genPList_->size();
  if (!ngenp) return;

  os << setprecision(2);
  os << " -- # GenParticle: " << ngenp << endl;
  os << setw(3)  << "ind"
     << setw(4)  << "st"
     << setw(9)  << "pdgId"
     << setw(10) << "eta"
     << setw(7)  << "phi"
     << setw(8)  << "pt"
     << setw(9)  << "energy"
     << setw(6)  << "moInd"
     << setw(32) << "moID"
     << setw(64) << "daughterID"
     << std::endl;
  int indx {0};
  for (const auto& gp: *genPList_) {
    std::ostringstream mID;
    const auto& m = gp.motherIndices;
    mID << "(";
    for (const auto mi: m) {
      if (mi >= ngenp) continue;
      const auto& mgp = genPList_->at(mi);
      // skip low energy partons 
      if (std::abs(mgp.pdgId) == 21 && mgp.energy <= 5) continue;
      mID << " " << mgp.pdgId; 
    }
    mID << ")";
    string ms = mID.str();
    if (!ms.length()) ms = " -";
    
    std::ostringstream dID;
    const auto& d = gp.daughtIndices;
    dID << "(";
    for (const auto di: d) {
      if (di >= ngenp) continue;
      const auto& dgp = genPList_->at(di);
      // skip low energy partons 
      if (std::abs(dgp.pdgId) == 21 && dgp.energy <= 5) continue;
      dID << " " << dgp.pdgId; 
    }
    dID << ")";
    string ds = dID.str();
    if (!ds.length()) ds = " -";

    os << setw(3)  << indx++
       << setw(4)  << gp.status
       << setw(9)  << gp.pdgId
       << setw(10) << gp.eta
       << setw(7)  << gp.phi
       << setw(8)  << gp.pt
       << setw(9)  << gp.energy
       << setw(6)  << gp.motherIndex 
       << setw(32) << ms 
       << ds
       << endl;
  }
}
/*int GenAnaBase::getMotherId(const vhtm::GenParticle& gp, int& mmid) const {
  //std::cout<<"debug_1"<<std::endl;
  int pdgid = gp.pdgId;
  auto m = gp.motherIndices;
  if (m.size() < 1) return -1;
  //std::cout<<"debug_2"<<std::endl;
  int indx = m.at(0);
  std::cout<<"Mother's indx: "<<indx<<std::endl;;
  const auto& mgp = genPList_->at(indx);
  mmid = mgp.pdgId;
  std::cout<<"Mother's pid: "<<mmid<<std::endl;
  //std::cout<<"debug_3"<<std::endl;
  while (mmid == pdgid) {
    //std::cout<<"debug_4"<<std::endl;
    m = mgp.motherIndices;
    indx = m.at(0);
    std::cout<<"GMother's indx: "<<indx<<std::endl;;
    const auto& mgp = genPList_->at(indx);
    mmid = mgp.pdgId;
    std::cout<<"Gmoms pid: "<<mmid<<std::endl;
  }
  return indx;
  }*/
int GenAnaBase::getMotherId(const vhtm::GenParticle& gp, int& mmid) const {
  int pdgid = gp.pdgId;
  auto m = gp.motherIndices;
  if (m.size() < 1) return -1;
  int indx = m.at(0);
  vhtm::GenParticle mgp = genPList_->at(indx);
  mmid = mgp.pdgId;
  while (mmid == pdgid) {
    m = mgp.motherIndices;
    indx = m.at(0);
    mgp = genPList_->at(indx);
    mmid = mgp.pdgId;
  }
  return indx;
}
