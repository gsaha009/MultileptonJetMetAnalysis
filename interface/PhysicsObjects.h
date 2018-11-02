#ifndef __AnalysisSpace_TreeMaker_PhysicsObjects_h
#define __AnalysisSpace_TreeMaker_PhysicsObjects_h

#include <vector>
#include <map>
#include <string>
#include "TLorentzVector.h"
#include "TObject.h"

namespace vhtm {
  class Candidate: public TObject {
  public:
    Candidate();
    Candidate(float pt, float eta, float phi);
    virtual ~Candidate() {}

    float pt;
    float eta;
    float phi;
   
    ClassDef(Candidate, 1)
  };

  class PackedPFCandidate: public TObject {
  public:
    PackedPFCandidate();
    virtual ~PackedPFCandidate() {}

    float pt;
    float eta;
    float phi;
    float phiAtVtx;
    float energy;
    bool trackHighPurity;
    
    int pdgId;
    int charge;
    
    double vx;
    double vy;
    double vz;
   
    int fromPV;
    // w.r.t PV (essentially pv[0])
    float dxy;
    float dz;
    float dxyError;
    float dzError;

    // wrt PV with highest sumPtTracks
    float dxyEV;
    float dzEV;
    float dzAssociatedPV;

    int numberOfHits;
    int numberOfPixelHits;
    int pixelLayersWithMeasurement;
    int stripLayersWithMeasurement;
    int lostInnerHits;

    std::map<std::string, std::vector<double>> isolationMap;

    ClassDef(PackedPFCandidate, 1)
  };
 
  class Event: public TObject {
  public:
    Event();
    virtual ~Event() {}
  
    unsigned long run;
    unsigned long event;
    unsigned long lumis;
    unsigned long bunch;
    unsigned long orbit;
    double time;
    bool isdata;
  
    bool hasPrimaryVertex;
    bool passHBHENoiseFilter;
  
    double rho;
    double rhoNeutral;
    int nvtx; 

    double fGridRhoAll;
    double fGridRhoFastjetAll;
    double fGridRhoFastjetAllCalo;
    double fGridRhoFastjetCentralCalo;
    double fGridRhoFastjetCentralChargedPileUp;
    double fGridRhoFastjetCentralNeutral;
    std::vector<int> nPU;
    std::vector<int> bunchCrossing;
    std::vector<int> trueNInt;

    ClassDef(Event, 1)
  };
  class GenEvent: public TObject {
  public:
    GenEvent();
    virtual ~GenEvent() {}
  
    double evtWeight;
    unsigned int processID;
    double ptHat;
       int nMEPartons;
    double qScale;
    double alphaQCD;
    double alphaQED;

    double qScalePdf;
    double x1Pdf;
    double x2Pdf;
       int id1Pdf;
       int id2Pdf;

    double lheHt;
       int lheNOutPartons;

    std::vector<double> pdfWeightList;
    std::vector<double> alphasWeightList;
    std::vector<double> qcdScaleWeightList;

    std::vector<double> pdfWeights;
  
    ClassDef(GenEvent, 1)
  };
  
  class Electron: public TObject {
  public:
    Electron();
    ~Electron() {}

    double eta;
    double phi;
    double pt;
    bool ecalDriven;
    bool isGap;
    bool hasGsfTrack;
    double trackPt;
    double energy;
    double caloEnergy;
    int charge;
    int pixHits;
    int trkHits;
    int nValidHits;
    double trkD0;
    double trkDz;
    // ID variables
    double hoe;
    double eop;
    double sigmaEtaEta;
    double sigmaIEtaIEta;
    double deltaPhiTrkSC;
    double deltaEtaTrkSC;
    int classif;
  
    // Vertex
    double vx;
    double vy;
    double vz;

    // SC associated with electron
    double scEn;
    double scEta;
    double scPhi;
    double scET;
    double scRawEnergy;
    float BDT;  
    float BDTpreComp;  
    // Vertex association variables
    double dxyPV;
    double dzPV;
    double vtxDist3D;
    int vtxIndex;
    double vtxDistZ;
    float relIso;
    float pfRelIso;
    // PFlow isolation variable
    float chargedHadronIso;
    float neutralHadronIso;
    float photonIso;

    float sumChargedHadronPt;
    float sumNeutralHadronEt;
    float sumPhotonEt;
    float sumPUPt;

    int missingHits;

    double dB;
    double edB;    

    double dB3D;
    double edB3D;    

    int nBrems;
    float fbrem;
  
    float hasMatchedConv;
    bool mvaPreselection;
    bool isTriggerElectron;
    int fidFlag;
    std::map<std::string, float> idmap;
    bool passMediumId;
    bool passTightId;
    int mvaCategory;
    float ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values;
    std::map<std::string,std::vector<double>> isolationMap;
    std::vector<TLorentzVector> associatedPackedPFCandidatesP4;

    ClassDef(Electron, 1)
  };
  class GenParticle: public TObject {
  public:
    GenParticle();
    ~GenParticle() {}
  
    double eta;
    double phi;
    double p;
    double px;
    double py;
    double pz;
    double pt;
    double energy;
    int pdgId;
    double vx;
    double vy;
    double vz;
    int status;
    double charge;
    int numDaught;
    int numMother;
    int motherIndex;
    std::vector<int> motherIndices;
    std::vector<int> daughtIndices;
  
    ClassDef(GenParticle, 1)
  };
  class GenJet: public TObject {
  public:
    GenJet();
    ~GenJet() {}
  
    double eta;
    double phi;
    double p;
    double pt;
    double energy;
    double emf;
    double hadf;
  
    ClassDef(GenJet, 1)
  };
  class MET: public TObject {
  public:
    MET();
    ~MET() {}
  
    double met;
    double metphi;
    double sumet;
    
    double metJESUp;
    double metJESDn;

    double metuncorr;
    double metphiuncorr;
    double sumetuncorr;
  
    ClassDef(MET, 1)
  };
  class Tau: public TObject {
  public:
    Tau();
    ~Tau() {}
  
    double eta;
    double phi;
    double pt;
    double energy;
    int charge;
    double mass;
  
    double dxyPV;
    double dxyPVError;
    double dzPV;
    double dzPVError;
    int vtxIndex;
    double vtxDxy;
    double vtxDz;

    // Leading particle pT
    double leadChargedParticlePt;
    double leadNeutralParticlePt;
    double leadParticlePt;

    std::vector<vhtm::Candidate> sigChHadList;
    std::vector<vhtm::Candidate> sigNeHadList;
    std::vector<vhtm::Candidate> sigGammaList;
    std::vector<vhtm::Candidate> isoChHadList;
    std::vector<vhtm::Candidate> isoNeHadList;
    std::vector<vhtm::Candidate> isoGammaList;

    float ptSumChargedHadronsIsoCone;
    float ptSumNeutralHadronsIsoCone;
    float ptSumPhotonsIsoCone;

    // tau id. discriminators
    float decayModeFinding;
    float decayModeFindingNewDMs;
    //float decayModeFindingOldDMs;
    int decayMode;
    
    // discriminators against electrons/muons
    float againstMuonLoose3;
    float againstMuonTight3;
    
    float againstElectronVLooseMVA;
    float againstElectronLooseMVA;
    float againstElectronMediumMVA;
    float againstElectronTightMVA;
    float againstElectronVTightMVA;
    
    float byLooseCombinedIsolationDeltaBetaCorr3Hits;
    float byMediumCombinedIsolationDeltaBetaCorr3Hits;
    float byTightCombinedIsolationDeltaBetaCorr3Hits;
    float byCombinedIsolationDeltaBetaCorrRaw3Hits;
    
    float byVLooseIsolationMVArun2v1DBoldDMwLT;
    float byLooseIsolationMVArun2v1DBoldDMwLT;
    float byMediumIsolationMVArun2v1DBoldDMwLT;
    float byTightIsolationMVArun2v1DBoldDMwLT;
    float byVTightIsolationMVArun2v1DBoldDMwLT;
    
    float chargedIsoPtSum;
    float neutralIsoPtSum;
    float puCorrPtSum;
    
     // kinematic variables for PFJet associated to PFTau
    double jetPt;
    double jetEta;
    double jetPhi;
     float emFraction;
    double vx;
    double vy;
    double vz;
  
    double zvertex;
    double dxySig;
  
    ClassDef(Tau, 1)
  };
  class Muon: public TObject {
  public:
    Muon();
    ~Muon() {}
    bool isGlobalMuon;
    bool isTrackerMuon;
    bool isPFMuon;
    bool isghostCleaned;
    bool passTrackerhighPtid;
    bool isLooseMuon;
    bool isMediumMuon;
    bool isGoodMedMuon;
    bool isTightMuon;
    float chi2LocalPosition;
    float trkKink;

    double eta;
    double phi;
    double pt;

    double p;
    double energy;
    int charge;
    double trkD0;

    double trkDz;
    int muonBestTrackType;  
   
    double globalChi2;
    double tkNChi2;

    float pfChargedIsoR03;
    float pfChargedHadIsoR03;
    float pfNeutralHadIsoR03;
    float pfPhotonIso03;    
    float sumPUPt03;
    float pfRelIso03;

    float sumChargedParticlePt;
    float sumChargedHadronPt;
    float sumNeutralHadronEt;
    float sumPhotonEt;
    float sumPUPt;
    float pfRelIso04;

    int passID;
    double dxyPV;
    double dzPV;
    double vtxDist3D;
    int vtxIndex;
    double vtxDistZ;
    int pixHits;
    int trkHits;
    int muoHits;
    int matches;
    int trackerLayersWithMeasurement;

    double vx;
    double vy;
    double vz;

    double dB;  // 2D
    double edB;    
  
    double dB3D;
    double edB3D;    

    // UW Recommendation
    bool isGlobalMuonPromptTight;
    bool isAllArbitrated;
    int nChambers;
    int nMatches;
    int nMatchedStations;

    unsigned int stationMask;
    unsigned int stationGapMaskDistance;
    unsigned int stationGapMaskPull;

    bool muonID;
   
    int nSegments;

    std::map<std::string,std::vector<double>> isolationMap;

    ClassDef(Muon, 1)
  };
  class Jet: public TObject {
  public:
    Jet();
    ~Jet() {}
    double eta;
    double phi;
    double pt;
    double pt_raw;
    double energy;
    double energy_raw;
    double jecUnc;
    double resJEC;
    int partonFlavour;

    float chargedEmEnergyFraction;
    float chargedHadronEnergyFraction;
    float chargedMuEnergyFraction;
    float electronEnergyFraction;
    float muonEnergyFraction;
    float neutralEmEnergyFraction;
    float neutralHadronEnergyFraction;
    float photonEnergyFraction;
    int chargedHadronMultiplicity;
    int chargedMultiplicity;
    int electronMultiplicity;
    int muonMultiplicity;
    int neutralHadronMultiplicity;
    int neutralMultiplicity;
    int photonMultiplicity;
    int nConstituents;

    float combinedSecondaryVertexBTag;
    float combinedInclusiveSecondaryVertexBTag;
    float combinedInclusiveSecondaryVertexV2BJetTags;
    float pfCombinedInclusiveSecondaryVertexV2BJetTags;

    std::map<std::string, float> discrimap;
    
    float jpumva;

    int passLooseID;
    int passTightID;
  
    ClassDef(Jet, 1)
  };
  class Vertex: public TObject {
  public:
    Vertex();
    virtual ~Vertex() {}
  
    double x;
    double y;
    double z;
    double xErr;
    double yErr;
    double zErr;
    double rho;
    double chi2;
    double ndf;
    bool isfake;
    bool isvalid;

    ClassDef(Vertex, 1)
  };
  class GenMET: public TObject {
  public:
    GenMET();
    virtual ~GenMET() {}
  
    double met;
    double metphi;
    double sumet;
  
    ClassDef(GenMET, 1)
  };
  class TriggerObject: public TObject {
  public:
    TriggerObject();
    virtual ~TriggerObject() {}
  
    double energy;
    double pt;
    double eta;
    double phi;
  
    std::map<std::string, unsigned int> pathList;
  
    ClassDef(TriggerObject, 1)
  };
  class Photon : public TObject {
  public:
    Photon();
    virtual ~Photon() {}
    
    double et;
    double eta;
    double clusterEta;
    double phi;
    double clusterPhi;
    double energy;
    double theta; 
    double vx;
    double vy;
    double vz;
    
    double scEnergy;
    double scEta;
    double scPhi;
    double scSize;
    double scEtaWidth;
    double scPhiWidth;
    double scEt;
    double scRawEnergy;
    double scx;
    double scy;
    double scz; 
    double isoEcalRecHit03;
    double isoHcalRecHit03;
    double isoSolidTrkCone03;
    double isoHollowTrkCone03;
    int nTrkSolidCone03;
    int nTrkHollowCone03;
    
    double isoEcalRecHit04;
    double isoHcalRecHit04;
    double isoSolidTrkCone04;
    double isoHollowTrkCone04;
    int nTrkSolidCone04;
    int nTrkHollowCone04;
    
    bool isEB;
    bool isEE; 
    bool isEBGap;
    bool isEEGap;
    bool isEBEEGap;
    int fidFlag;
    
    bool passElectronVeto;
    bool hasPixelSeed;
    double ecalIso;
    double hcalIso;
    double trackIso;
    double chargedHadIso;
    double neutralHadIso;
    double photonIso;
    
    double r9;
    double hoe;
    double sigmaEtaEta;
    double sigmaIEtaIEta;
    double e1x5;
    double e2x5; 
    double e3x3;
    double e5x5;
    double r1x5;
    double r2x5;
    double maxEnergyXtal;
    
    bool hasConversionTracks;
    int nTracks;
    bool isConverted;
    double pairInvMass;
    double pairCotThetaSeparation;
    double pairPx;
    double pairPy;
    double pairPz;
    double conv_vx;
    double conv_vy;
    double conv_vz;
    double eovp;
    double zpv;
    double distOfMinApproach;
    double dPhiTracksAtVtx;
    double dPhiTracksAtEcal;
    double dEtaTracksAtEcal;  
    
    std::map<std::string,std::vector<double>> isolationMap;

    ClassDef(Photon, 1)
  };
}
#endif
