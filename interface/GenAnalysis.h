#ifndef __GenAnalysis__hh
#define __GenAnalysis__hh

//#include "AnaBase.h"
#include "GenAnaBase.h"

class GenAnalysis: public GenAnaBase {
    
public:
  enum EventType {
    mmmt = 0, mmet, eemt, eeet
  };
  GenAnalysis();
  ~GenAnalysis();

  virtual void bookHistograms(std::unique_ptr<TFile>& histf) override;
  void clearLists();
  virtual void analyze(const std::unique_ptr<TFile>& histf) override;
  virtual bool filter() const override;

  bool genOk() const {return genPass_;}

  std::vector <vhtm::GenParticle> leptonBox;
  std::vector <vhtm::GenParticle> neutrinoBox;
  std::vector <vhtm::GenParticle> quarkBox;
  std::vector <vhtm::GenParticle> hP_Decay;
  std::vector <vhtm::GenParticle> hM_Decay;
    
public:
  bool genPass_ {false};
};
#endif
