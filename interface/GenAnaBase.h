#ifndef __GenAnaBase__hh
#define __GenAnaBase__hh

#include <memory>
#include <vector>
#include <sstream>

#include <TFile.h>

#include "PhysicsObjects.h"

class GenAnaBase {
    
public:
  virtual ~GenAnaBase();

  void dumpEvent(std::ostream& os=std::cout) const;
  int getMotherId(const vhtm::GenParticle& gp, int& mmid) const;
  void setEvent(const std::vector<vhtm::GenParticle>* genPList) {genPList_ = genPList;}

  virtual void bookHistograms(std::unique_ptr<TFile>& histf) = 0;
  virtual void analyze(const std::unique_ptr<TFile>& histf) = 0;
  virtual bool filter() const = 0;
    
public:
  const std::vector<vhtm::GenParticle>* genPList_;
};
#endif
