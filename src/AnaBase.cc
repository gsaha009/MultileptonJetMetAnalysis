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
#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "AnaUtil.h"
#include "AnaBase.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::setiosflags;
using std::resetiosflags;

using namespace vhtm;

// -----------
// Constructor
// -----------
AnaBase::AnaBase()
#ifdef TRUE_CPP14
  : chain_(std::make_unique<TChain>("treeCreator/vhtree")),
#else
  : chain_(std::unique_ptr<TChain>(new TChain("treeCreator/vhtree"))),
#endif
    eventList_(new vector<vhtm::Event>()),
    vertexList_(new vector<vhtm::Vertex>()),    
    genEventList_(new vector<vhtm::GenEvent>()),   
    tauList_(new vector<vhtm::Tau>()),    
    electronList_(new vector<vhtm::Electron>()),
    muonList_(new vector<vhtm::Muon>()),   
    photonList_(new vector<vhtm::Photon>()),
    packedPFCandidateList_(new vector<vhtm::PackedPFCandidate>()),
    jetList_(new vector<vhtm::Jet>()),    
    metList_(new vector<vhtm::MET>()),    
    corrmetList_(new vector<vhtm::MET>()),    
    puppimetList_(new vector<vhtm::MET>()),    
    genParticleList_(new vector<vhtm::GenParticle>()),  
    genJetList_(new vector<vhtm::GenJet>()),   
    genMetList_(new vector<vhtm::GenMET>()),
    triggerObjList_(new vector<vhtm::TriggerObject>()),
    l1physbits_(new vector<int>()),
    l1techbits_(new vector<int>()),
    hltpaths_(new vector<string>()),
    hltresults_(new vector<int>()),
    hltprescales_(new vector<int>())
{
  cout << setiosflags(ios::fixed); 
  fileList_.clear();
  brList_.clear();
  puWtList_.clear();
  trigPathList_.clear();
}
// ----------
// Destructor
// ----------
AnaBase::~AnaBase() 
{
  clearEvent();

  if (eventList_ != nullptr)             delete eventList_;
  if (vertexList_ != nullptr)            delete vertexList_;
  if (genEventList_ != nullptr)          delete genEventList_;
  if (tauList_ != nullptr)               delete tauList_;
  if (electronList_ != nullptr)          delete electronList_;
  if (muonList_ != nullptr)              delete muonList_;
  if (photonList_ != nullptr)            delete photonList_;
  if (packedPFCandidateList_ != nullptr) delete packedPFCandidateList_;
  if (jetList_ != nullptr)               delete jetList_;
  if (metList_ != nullptr)               delete metList_;
  if (corrmetList_ != nullptr)           delete corrmetList_;
  if (puppimetList_ != nullptr)          delete puppimetList_;
  if (genParticleList_ != nullptr)       delete genParticleList_;
  if (genJetList_ != nullptr)            delete genJetList_;
  if (genMetList_ != nullptr)            delete genMetList_;
  if (triggerObjList_ != nullptr)        delete triggerObjList_;
}
// ------------------------
// Clear the clones arrays
// ------------------------
void AnaBase::clearEvent() 
{
  // do not own the objects
  eventList_->clear();
  genEventList_->clear();
  vertexList_->clear();
  tauList_->clear();
  electronList_->clear();
  muonList_->clear();
  photonList_->clear();
  packedPFCandidateList_->clear();
  jetList_->clear();
  metList_->clear();
  corrmetList_->clear();
  puppimetList_->clear();
  genParticleList_->clear();
  genJetList_->clear();
  genMetList_->clear();
  triggerObjList_->clear();

  l1physbits_->clear();
  l1techbits_->clear();
  hltpaths_->clear();
  hltresults_->clear();
  hltprescales_->clear();
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaBase::beginJob() 
{ 
  if (isMC_ && usePUWt_ && !readPileUpHist()) return false;

  // Open the output ROOT file
  //histf_ = std::make_unique<TFile>(histFile_.c_str(), "RECREATE");

  TFile* f = TFile::Open(histFile_.c_str(), "RECREATE");
  histf_.reset(std::move(f));

  setAddresses();
  nEvents_ = static_cast<int>(chain_->GetEntries()); 
  if (nEvents_ <= 0) {
    cerr << "******* nEvents = " << nEvents_ << ", returning!" << endl;
    return false;
  }
  if (maxEvt_ > 0) nEvents_ = std::min(nEvents_, maxEvt_);
  cout << " >>> nEvents = " << nEvents_ << endl;

  openFiles();

  return true;
}
void AnaBase::endJob() 
{
  //  closeFiles();
}
void AnaBase::closeHistFile(){
  histf_->cd();
  histf_->Write();
  histf_->Close();
}
// ------------------------------------
// Get Run number for the present event
// ------------------------------------
int AnaBase::getRunNumber() const 
{
  return eventList_->at(0).run;
}    
double AnaBase::lumiWt(double evtWeightSum, bool verbose) const 
{
  double nevt = (evtWeightSum > -1) ? evtWeightSum : AnaUtil::cutValue(lumiWtMap(), "nevents");
  if (!verbose) 
    std::cout << "-- intLumi: " << AnaUtil::cutValue(lumiWtMap(), "intLumi") 
	      << " xsec: " << AnaUtil::cutValue(lumiWtMap(), "xsec") 
	      << " nevt: " << nevt << std::endl; 
  return (AnaUtil::cutValue(lumiWtMap(), "intLumi") * AnaUtil::cutValue(lumiWtMap(), "xsec") / nevt); 
}

// ---------------------------------
// Add input Root files to the chain
// ---------------------------------
int AnaBase::setInputFile(const string& fname) 
{
  auto found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())) {
    cerr << ">>> Warning: File <<" << fname << ">> was not found!!" << endl;
    return static_cast<int>(chain_->GetEntries()); 
  }
  chain_->AddFile(fname.c_str(), -1);
  return static_cast<int>(chain_->GetEntries()); 
}
// ---------------------------------------
// Get total number of events in the chain
// --------------------------------------
int AnaBase::getEntries() const 
{
  return static_cast<int>(chain_->GetEntries());
}

// ------------------------------------------------------
// Open the output file with a global filehandle, C++ way
// ------------------------------------------------------
bool AnaBase::openFiles() 
{
  fLog_.open(logFile_.c_str(), ios::out);
  if (!fLog_) {
    cerr << "File: " << logFile_ << " could not be opened!" << endl;
    return false;
  }
  fLog_ << setiosflags(ios::fixed);

  evLog_.open(evFile_.c_str(), ios::out);
  if (!evLog_) {
    cerr << "File: " << evFile_ << " could not be opened!" << endl;
    return false;
  }
  evLog_ << setiosflags(ios::fixed);

  selEvLog_.open(selEvFile_.c_str(), ios::out);
  if (!selEvLog_) {
    cerr << "File: " << selEvFile_ << " could not be opened!" << endl;
    return false;
  }
  selEvLog_ << setiosflags(ios::fixed);

  return true;
}
// ------------------------
// Close the output file
// ------------------------
void AnaBase::closeFiles() 
{
  if (fLog_) {
    fLog_ << resetiosflags(ios::fixed); 
    fLog_.close();
  }
  if (evLog_) {
    evLog_ << resetiosflags(ios::fixed); 
    evLog_.close();
  }
  if (selEvLog_) {
    selEvLog_ << resetiosflags(ios::fixed); 
    selEvLog_.close();
  }
  closeHistFile();
}
void AnaBase::setAddresses() 
{
  if (branchFound("Event"))    chain_->SetBranchAddress("Event", &eventList_);
  if (branchFound("Vertex"))   chain_->SetBranchAddress("Vertex", &vertexList_);
  if (branchFound("Tau"))      chain_->SetBranchAddress("Tau", &tauList_);
  if (branchFound("Electron")) chain_->SetBranchAddress("Electron", &electronList_);
  if (branchFound("Muon"))     chain_->SetBranchAddress("Muon", &muonList_);
  if (branchFound("Photon"))   chain_->SetBranchAddress("Photon", &photonList_);
  if (readPFObject_ && branchFound("PackedPFCandidate")) 
                               chain_->SetBranchAddress("PackedPFCandidate", &packedPFCandidateList_);
  if (branchFound("Jet"))      chain_->SetBranchAddress("Jet", &jetList_);
  if (branchFound("MET"))      chain_->SetBranchAddress("MET", &metList_);
  if (branchFound("corrMET"))  chain_->SetBranchAddress("corrMET", &corrmetList_);
  if (branchFound("puppiMET")) chain_->SetBranchAddress("puppiMET", &puppimetList_);
  if (readTrigObject_ && branchFound("TriggerObject")) 
                               chain_->SetBranchAddress("TriggerObject", &triggerObjList_);
  if (isMC_) {
    if (branchFound("GenEvent"))      chain_->SetBranchAddress("GenEvent", &genEventList_);
    if (branchFound("GenJet"))      chain_->SetBranchAddress("GenJet", &genJetList_);
    if (readGenInfo_) {
      if (branchFound("GenParticle")) chain_->SetBranchAddress("GenParticle", &genParticleList_);
      //      if (branchFound("GenJet"))      chain_->SetBranchAddress("GenJet", &genJetList_);
    }
  }
  // Now the trigger variables
  if (branchFound("l1physbits"))   chain_->SetBranchAddress("l1physbits", &l1physbits_);
  if (branchFound("l1techbits"))   chain_->SetBranchAddress("l1techbits", &l1techbits_);
  if (branchFound("hltresults"))   chain_->SetBranchAddress("hltresults", &hltresults_);
  if (branchFound("hltprescales")) chain_->SetBranchAddress("hltprescales", &hltprescales_);
  if (branchFound("hltpaths"))     chain_->SetBranchAddress("hltpaths", &hltpaths_);
}
bool AnaBase::branchFound(const string& b)
{
  TBranch* branch = chain_->GetBranch(b.c_str());  // Get branch pointer                                                                             
  if (branch == nullptr) {
    cout << ">>> SetBranchAddress: <" << b << "> not found!" << endl;
    return false;
  }
  cout << ">>> SetBranchAddress: <" << b << "> found!" << endl;
  brList_.push_back(b);
  return true;
}
int AnaBase::getEntry(int lflag) const
{
  int nbytes = 0;
  for (const auto& v: brList_) {
    TBranch* branch = chain_->GetBranch(v.c_str());
    if (branch == nullptr) {
      cout << ">>> Branch: " << v << " not found!" << endl;
      continue;
    }
    nbytes += branch->GetEntry(lflag);
  }
  return nbytes;
}
// not used yet
void AnaBase::enableBranches() 
{
  chain_->SetBranchStatus("*", kFALSE); // Disable all branches
  for (const auto& v: brList_)
    chain_->SetBranchStatus(v.c_str(), kTRUE);
}
bool AnaBase::readJob(const string& jobFile, int& nFiles)
{
  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "==> Input File: <<" << jobFile << ">> could not be opened!" << endl;
    return false;
  }

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    vector<string> tokens;
    AnaUtil::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);

    const string& key   = tokens.at(0);
    const string& value = tokens.at(1);
    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      AnaUtil::tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }
    else if (key == "readTrigObject") 
      readTrigObject_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "readGenInfo") 
      readGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "readPFObject") 
      readPFObject_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "useTrigger") 
      useTrigger_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "useLumiWt") 
      useLumiWt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "usePUWt") 
      usePUWt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "logFile")
      logFile_ = value;
    else if (key == "eventFile")
      evFile_  = value;
    else if (key == "selEventFile")
      selEvFile_  = value;
    else if (key == "logOption") 
      logOption_ = strtol(value.c_str(), NULL, 2);
    else if (key == "maxEvent") 
      maxEvt_ = std::stoi(value.c_str());
    else if (key == "startEvent") 
      firstEvt_ = std::stoi(value.c_str());
    else if (key == "endEvent") 
      lastEvt_ = std::stoi(value.c_str());
    else if (key == "bunchX") 
      bunchCrossing_ = std::stoi(value.c_str());
    else if (key == "histFile") 
      histFile_ = value;
    else if (key == "puHistFile") 
      puHistFile_ = value;
    else if (key == "puHistogram") 
      puHistogram_ = value;
    else if (key == "useTrueNInt") 
      useTrueNInt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "trigPathList") 
      AnaUtil::buildList(tokens, trigPathList_);
    else if (key == "inputFile") 
      AnaUtil::buildList(tokens, fileList_);
    else if (key == "eventId" && tokens.size() == 4) 
      AnaUtil::buildMap(tokens, eventIdMap_);
    else {
      if (0) cout << "==> " << line << endl;
      AnaUtil::storeCuts(tokens, hmap_);
    }
  }
  // Close the file
  fin.close();

  if (!isMC_) usePUWt_ = false;
  if (!isSignal_) readGenInfo_ = false;

  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    int nevt = setInputFile(fname);
    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }
  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }

  return true;
}
void AnaBase::printJob(ostream& os) const
{
  os << "       datatype: " << ((isMC_) ? "mc" : "data") << endl
     << "    readGenInfo: " << readGenInfo_ << endl
     << " readTrigObject: " << readTrigObject_ << endl
     << "   readPFObject: " << readPFObject_ << endl
     << "        logFile: " << logFile_ << endl 
     << "      eventFile: " << evFile_ << endl
     << "       histFile: " << histFile_ << endl
     << "      selEvFile: " << selEvFile_ << endl
     << "      useLumiWt: " << std::boolalpha << useLumiWt_ << endl
     << "        usePUWt: " << std::boolalpha << usePUWt_ << endl
     << "     puHistFile: " << puHistFile_ << endl
     << "    puHistogram: " << puHistogram_ << endl
     << "    useTrueNInt: " << std::boolalpha << useTrueNInt_ << endl
     << "     useTrigger: " << std::boolalpha << useTrigger_ << endl
     << "      logOption: " << logOption_ << endl
     << "       maxEvent: " << maxEvt_ 
     << endl;

  // Trigger Path List
  if (useTrigger_) 
    AnaUtil::showList(trigPathList_, ">>> INFO. Trigger Paths used:", os);

  // InputFiles
  if (chain_ != nullptr) {
    TObjArray *fileElements = chain_->GetListOfFiles();
    os << ">>> INFO. nFiles: " << fileElements->GetEntries() 
       << ", Files to analyse:" 
       << endl;
    TIter next(fileElements);
    TChainElement *chEl = 0;
    while (( chEl = dynamic_cast<TChainElement*>(next()) ))
      os << chEl->GetTitle() 
         << endl;
  }
  else
    AnaUtil::showList(fileList_, ">>> INFO. inputFiles:", os);

  // EventID 
  AnaUtil::showMap<string, int>(eventIdMap_, "Event List:", os);

  // Cuts
  AnaUtil::showCuts(hmap_, os);
}
// Collect Object information
// event vertex  
void AnaBase::findVtxInfo(vector<Vertex>& list, Options& op, ostream& os) {
  if (nvertex() < 1) return;

  if (op.verbose)
    os << "=>> Vertices: " << nvertex() << endl
       << "indx     ndf     dxy       z  chi2    sbit" 
       << endl; 
  // already in correct order, no need to sort
  int indx {0};
  for (const auto& vtx: *vertexList_) {
    double dxy = std::sqrt(pow(vtx.x, 2) + pow(vtx.y, 2));
    int sbit {0};
    if (vtx.ndf <= AnaUtil::cutValue(vtxCutMap(), "ndf"))      sbit |= (1 << 0);
    if (vtx.rho > AnaUtil::cutValue(vtxCutMap(), "Rho") )      sbit |= (1 << 1);
    if (std::abs(vtx.z) > AnaUtil::cutValue(vtxCutMap(), "z")) sbit |= (1 << 2);
    if (vtx.isfake)                                            sbit |= (1 << 3);
    if (op.verbose) {
      bool pp = (op.printselected && sbit) ? false : true;
      if (pp) {
        os << setprecision(2)
           << setw(4) << indx
           << setw(8) << vtx.ndf
           << setw(8) << dxy
           << setw(8) << vtx.z 
           << setw(8) << vtx.chi2;
        AnaUtil::bit_print(sbit, 8, os);
      }
    }
    if (op.usesbit && sbit) continue;
    list.push_back(vtx);
  }
}
void AnaBase::findTriggerObjectInfo(vector<TriggerObject>& list) {
  for (const auto& obj: *triggerObjList_)
    list.push_back(obj);
}
TVector3 AnaBase::findLeptonVtx(int index, bool& isGoodVtx) {
  const Vertex& vtx = vertexList_->at(index);

  isGoodVtx = false;
  double dxy = std::sqrt(pow(vtx.x, 2) + pow(vtx.y, 2));
  if (vtx.ndf > AnaUtil::cutValue(vtxCutMap(), "ndf") && 
      dxy < AnaUtil::cutValue(vtxCutMap(), "dxy") && 
      std::fabs(vtx.z) < AnaUtil::cutValue(vtxCutMap(), "z")) isGoodVtx = true;

  TVector3 v(vtx.x, vtx.y, vtx.z); 
  return v;
}
int AnaBase::getMotherId(const GenParticle& gp, int& mmid) const {
  int pdgid = gp.pdgId;
  auto m = gp.motherIndices;
  if (m.size() < 1) return -1;
  int indx = m.at(0);
  auto mgp = genParticleList_->at(indx);
  mmid = mgp.pdgId;
  while (mmid == pdgid) {
    m = mgp.motherIndices;
    indx = m.at(0);
    mgp = genParticleList_->at(indx);
    mmid = mgp.pdgId;
  } 
  return indx;
}
int AnaBase::getMotherIdForQ(const GenParticle& gp, int& mmid) const {
  int pdgid = gp.pdgId;
  auto m = gp.motherIndices;
  if (m.size() < 1) return -1;
  int indx = m.at(0);
  if (indx < 0) return -1;
  //cout<<indx<<endl;
  auto mgp = genParticleList_->at(indx);
  mmid = mgp.pdgId;
  if (std::abs(mmid) == 2212) return -1;
  while (!(std::abs(mmid) == 24 || std::abs(mmid) == 23 || std::abs(mmid) == 37 || std::abs(mmid) == 9900041)) {
    m = mgp.motherIndices;
    indx = m.at(0);
    //cout<<"indx: "<<indx<<endl;
    if (indx < 0) return -1;
    mgp = genParticleList_->at(indx);
    mmid = mgp.pdgId;
    m.clear();
    if (std::abs(mmid) == 2212) return -1;
  }
  return indx;
}
void AnaBase::dumpGenInfo(ostream& os) const {
  if (!ngenparticle()) return;

  os << setprecision(2);
  os << " -- # GenParticle: " << ngenparticle() << endl;
  os << "indx    status    pdgId     eta      phi      pt     energy  moIndx"
     << "      moID                   daughterID"
     << endl;
  int indx {0};
  for (const auto& gp: *genParticleList_) {
    std::ostringstream mID;
    const auto& m = gp.motherIndices;
    for (int mi: m) {
      if (mi >= ngenparticle()) continue;
      const auto& mgp = genParticleList_->at(mi);
      mID << " " << mgp.pdgId; 
    }
    string ms = mID.str();
    if (!ms.length()) ms = " -";
    
    std::ostringstream dID;
    const auto& d = gp.daughtIndices;
    for (int di: d) {
      if (di >= ngenparticle()) continue;
      const auto& dgp = genParticleList_->at(di);
      auto energy = dgp.energy;
      auto pdgid = dgp.pdgId;
      if (std::abs(pdgid) == 21 && energy <= 10) continue;
      dID << " " << dgp.pdgId; 
    }
    auto ds = dID.str();
    if (!ds.length()) ds = " -";
    os << setw(4)  << indx++
       << setw(8)  << gp.status
       << setw(10) << gp.pdgId
       << setw(10) << gp.eta
       << setw(9)  << gp.phi
       << setw(9)  << gp.pt
       << setw(9)  << gp.energy
       << setw(8)  << gp.motherIndex 
       << setw(10) << ms 
       << ds
       << endl;
  }
}
bool AnaBase::readPileUpHist(bool verbose) {
  size_t found = puHistFile_.find(".root");
  if (found == string::npos) {
    cerr << ">>> Warning: <<" << puHistFile_ << ">> does not have .root extension!!" << endl;
    return false;
  }

  const char* fname = gSystem->ExpandPathName(puHistFile_.c_str());
  if (gSystem->AccessPathName(fname)) {
    cerr << ">>> Warning: File <<" << puHistFile_ << ">> not found!!" << endl;
    return false;
  }

  TFile file(fname);
  TH1F *h = dynamic_cast<TH1F*>(file.Get(puHistogram_.c_str()));
  if (!h) {
    cerr << ">>> Warning: Histogram <<" << puHistogram_ << ">> not found!!" << endl;
    return false;
  }

  int nx = h->GetXaxis()->GetNbins();
  std::cout.precision(2);
  for (int i = 0; i < nx; ++i) {
    auto wt = h->GetBinContent(i);
    if (verbose) 
    cout << setw(3) << i 
	 << setw(6) << wt
	 << endl;
    puWtList_.push_back(wt);
  }
  return true;
}
double AnaBase::wtPileUp(int& nPU, bool verbose) const {
  nPU = 0;
  const Event& evt = eventList_->at(0);
  //cout << "PU bx: " << evt.bunchCrossing.size() << endl;
  for (size_t i = 0; i < evt.bunchCrossing.size(); ++i) {
    if (verbose)
    cout << setw(3) << i 
         << setw(4) << evt.bunchCrossing.at(i) 
         << setw(4) << evt.trueNInt.at(i) 
         << setw(4) << evt.nPU.at(i) 
         << endl;
    if (evt.bunchCrossing.at(i) == 0) { // in-time
      nPU = (useTrueNInt_) ? evt.trueNInt.at(i) : evt.nPU.at(i);
      break;
    }
  }

  int nbins = puWtList_.size();
  if (nPU < 0) nPU = 0;
  if (nPU >= nbins) nPU = nbins - 1;

  return puWtList_.at(nPU);
}

bool AnaBase::isTriggered(bool check_prescale, bool verbose) const {
  bool flag {false};
  for (size_t i = 0; i < hltpaths_->size(); ++i) {
    string str = (*hltpaths_).at(i);
    bool found = false;
    for (const auto& v: trigPathList_) {
      if (str.find(v) != string::npos) {
	found = true;
        break;
      }
    }
    if (!found) continue;

    int prescl = (check_prescale) ? (*hltprescales_).at(i) : 1;
    int result = (*hltresults_).at(i);
    if (verbose) cout << ">>> HLT Path = " << str 
                      << ", fired=" << result
                      << ", prescale=" << prescl
                      << endl;
    if (result == 1 && prescl == 1) {
      flag = true;
      break;
    }
  }
  return flag;
}
void AnaBase::dumpTriggerPaths(ostream& os, bool check_prescale) const 
{
  os << "=> Trigger paths" << endl;
  os << setw(96) << "Path"
     << setw(8) << "prescl"
     << setw(8) << "result"
     << endl;
  for (uint i = 0; i < hltpaths_->size(); ++i) {
    const string& path_name = (*hltpaths_).at(i);
    int prescale     = (*hltprescales_).at(i);  
    int result       = (*hltresults_).at(i);  
    if ((check_prescale && prescale != 1) || result != 1) continue;
    os << setw(96) << path_name 
       << setw(8) << prescale 
       << setw(8) << result
       << endl;
  }
}
void AnaBase::dumpTriggerObjectInfo(const vector<TriggerObject>& list, ostream& os) const
{
  os << setprecision(2);
  os << "=> TriggerObjects: " << list.size() << endl;
  os << "Indx     Eta     Phi      Pt  Energy            =Trigger path list=" << endl;
  int indx {0};
  for (const auto& tobj: list) {
    os << setw(4) << indx 
       << setw(8) << tobj.eta
       << setw(8) << tobj.phi
       << setw(8) << tobj.pt
       << setw(8) << tobj.energy
       << endl;
    map<string, uint> path_list = tobj.pathList;
    for (auto jt = path_list.begin(); jt != path_list.end(); ++jt) {
      os << "\t\t\t\t\t" << jt->first << " " << jt->second << endl;
    }
  }
}
bool AnaBase::matchTriggerPath(const vector<string>& pathList, const string& path) const {
  bool result = false;
  for (const auto& lname: pathList) {
    if (path.find(lname) != string::npos) {
      result = true;
      break;
    }
  }
  return result;
}
double AnaBase::matchTriggerObject(const vector<TriggerObject>& trigObjList, 
                                   const TLorentzVector& obj, 
                                   const string& trigPath, 
                                   int trig_skip, 
                                   double maxPtDiff,
                                   int& trig_indx) const
{
  double dRmin {999}; 
  trig_indx = -1;
  double obj_pt = obj.Pt();
  int indx {0};
  for (auto it = trigObjList.begin(); it != trigObjList.end(); ++it,++indx) {
    if (indx == trig_skip) continue;
    const TriggerObject& trigObj = (*it);
    const map<string, uint>& path_map = trigObj.pathList;
    bool matched = false;
    for (const auto& el: path_map) {
      string path = el.first;
      int flag = el.second;
      if (path.find(trigPath) != string::npos && flag == 1) {
        matched = true;
        break;
      }
    }
    if (!matched) continue;

    // check deltaPt
    if (std::abs(trigObj.pt - obj_pt) > maxPtDiff) continue;

    TLorentzVector trigTL;
    trigTL.SetPtEtaPhiE(trigObj.pt, trigObj.eta, trigObj.phi, trigObj.energy);
    double dR = AnaUtil::deltaR(obj, trigTL);
    if (dR < dRmin) {
      dRmin = dR;
      trig_indx = indx;
    }
  }
  return dRmin;
}
