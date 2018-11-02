#include <climits>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "TDirectory.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "AnaUtil.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::setw;

using std::string;
using std::vector;
using std::pair;
using std::map; 
using std::unordered_map; 

namespace AnaUtil {
  void tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    
    while (string::npos != pos || string::npos != lastPos)  {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  }
  void bit_print(int value, int pos, ostream& os) {
    static const int INT_BIT = 4*CHAR_BIT; 
    int i, mask = 1 << 31; 
    
    if (pos > INT_BIT) i = INT_BIT; 
    for (i = 1; i <= (INT_BIT - pos); ++i) { 
      value <<= 1; 
    } 
    os.put(' ');
    for (i = 1; i <= pos; ++i) { 
      os.put(((value & mask) == 0) ? '0' : '1'); 
      value <<= 1; 
      if ((INT_BIT - pos + i) % CHAR_BIT == 0 && i != INT_BIT) os.put(' '); 
    } 
    os << endl; 
  }
  double deltaPhi(double phia, double phib) {
    double dphi = phia - phib;
    while (dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2 * TMath::Pi();
    return dphi;
  }
  double deltaPhi(const TLorentzVector& a, const TLorentzVector& b) {
    return deltaPhi(a.Phi(), b.Phi());
  }
  double deltaR(const TLorentzVector& a, const TLorentzVector& b) {
    double dphi = deltaPhi(a,b);
    double deta = a.Eta() - b.Eta();
    return std::sqrt(dphi * dphi + deta * deta);
  }
  bool sameObject(const TLorentzVector& lv1, const TLorentzVector& lv2) {
    //return (std::fabs(lv1.Pt() - lv2.Pt()) < 1e-10 && AnaUtil::deltaR(lv1, lv2) < 1e-10);
    return (std::fabs(lv1.Pt() - lv2.Pt()) < 1.0e-08 && lv1.DeltaR(lv2) < 1.0e-06);
  }
  double cutValue(const map<string, double>& m, const string& mkey) {
    if (m.find(mkey) == m.end()) {
      cerr << ">>> key: " << mkey << " not found in the map!" << endl;
      for (auto const& el: m)  
        cerr << el.first << ": " << setw(7) << el.second << endl;
      return -999;
    }
    //assert(m.find(cname) != m.end());
    return m.find(mkey)->second;
  }
  const map<string, double>& cutMap(const map<string, map<string, double>>& hmap, const string& mkey) {
    assert(hmap.find(mkey) != hmap.end());
    return hmap.find(mkey)->second;
  }
  void buildList(const vector<string>& tokens, vector<string>& list) {
    for (vector<string>::const_iterator it  = tokens.begin()+1; 
                                        it != tokens.end(); ++it) {
      list.push_back(*it);       
    }
  }
  void buildMap(const vector<string>& tokens, map<string, int>& hmap) {
    string key = tokens.at(1) + "-" + tokens.at(2) + "-" + tokens.at(3);
    hmap.insert({key, 1});
  }
  void buildMap(const vector<string>& tokens, unordered_map<string, int>& hmap) {
    string key = tokens.at(1) + "-" + tokens.at(2) + "-" + tokens.at(3);
    hmap.insert({key, 1});
  }
  void storeCuts(const vector<string>& tokens, map<string, map<string, double>>& hmap) {
    const string& key = tokens.at(0);
    auto pos = hmap.find(key);
    if (pos == hmap.end()) {
      map<string, double> m;
      for (auto it = tokens.begin()+1; it != tokens.end(); ++it) {
        // Split the line into words
        vector<string> cutstr;
        tokenize(*it, cutstr, "=");
        if (cutstr.size() < 2) continue;
        m.insert({cutstr.at(0), atof(cutstr.at(1).c_str())});
      }
      if (!m.empty()) hmap.insert({key, m});
    }
  }
  void showCuts(const map<string, map<string, double>>& hmap, ostream& os) {
    for (auto const& el: hmap) {
      os << ">>> " << el.first << endl; 
      auto const& cutm = el.second;
      os << std::setprecision(2);
      for (auto const& il: cutm)  
        os << setw(16) << il.first << ": " 
           << setw(7)  << il.second << endl;
    }
  }
  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  TH1* getHist1D(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (obj == nullptr) {
      /*      cerr << "**** getHist1D: Histogram for <" << hname 
	   << "> not found! (" 
	   << __FILE__ << ":" << __LINE__ << ")" 
	   << endl;*/
      return nullptr;
    }
    TH1* h = nullptr;
    if (obj->InheritsFrom("TH1D"))
      h = dynamic_cast<TH1D*>(obj);
    else if (obj->InheritsFrom("TH1C"))
      h = dynamic_cast<TH1C*>(obj);
    else if (obj->InheritsFrom("TH1K"))
      h = dynamic_cast<TH1K*>(obj);
    else if (obj->InheritsFrom("TH1S"))
      h = dynamic_cast<TH1S*>(obj);
    else if (obj->InheritsFrom("TH1I"))
      h = dynamic_cast<TH1I*>(obj);
    else
      h = dynamic_cast<TH1F*>(obj);
    
    if (h == nullptr) {
      cerr << "**** getHist1D: <" << hname 
  	 << "> may not be a 1D Histogram! (" 
  	 << __FILE__ << ":" << __LINE__ << ")" 
  	 << endl;
    }
    return h;
  }
  TH1* getHist1D(const string& hname) {
    return getHist1D(hname.c_str());
  }
  
  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  TH2* getHist2D(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (obj == nullptr) {
      cerr << "**** getHist2D: Histogram for <" << hname 
  	 << "> not found! (" 
  	 << __FILE__ << ":" << __LINE__ << ")" 
  	 << endl;
      return nullptr;
    }
    
    TH2* h = nullptr;
    if (obj->InheritsFrom("TH2D"))
      h = dynamic_cast<TH2D*>(obj);
    else if (obj->InheritsFrom("TH2C"))
      h = dynamic_cast<TH2C*>(obj);
    else if (obj->InheritsFrom("TH2S"))
      h = dynamic_cast<TH2S*>(obj);
    else if (obj->InheritsFrom("TH2I"))
      h = dynamic_cast<TH2I*>(obj);
    else
      h = dynamic_cast<TH2F*>(obj);
    
    if (h == nullptr) {
      cerr << "**** getHist2D: <" << hname 
  	 << "> may not be a 2D Histogram! (" 
  	 << __FILE__ << ":" << __LINE__ << ")" 
  	 << endl;
    }
    return h;
  }
  TH2* getHist2D(const string& hname) {
    return getHist2D(hname.c_str());
  }
  // ---------------------------------------------
  // Convenience routine for filling 3D histograms
  // ---------------------------------------------
  TH3* getHist3D(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (obj == nullptr) {
      cerr << "**** getHist3D: Histogram for <" << hname 
  	 << "> not found! (" 
  	 << __FILE__ << ":" << __LINE__ << ")" 
  	 << endl;
      return nullptr;
    }
    
    TH3* h = nullptr;
    if (obj->InheritsFrom("TH3D"))
      h = dynamic_cast<TH3D*>(obj);
    else if (obj->InheritsFrom("TH3C"))
      h = dynamic_cast<TH3C*>(obj);
    else if (obj->InheritsFrom("TH3S"))
      h = dynamic_cast<TH3S*>(obj);
    else if (obj->InheritsFrom("TH3I"))
      h = dynamic_cast<TH3I*>(obj);
    else
      h = dynamic_cast<TH3F*>(obj);
    
    if (h == nullptr) {
      cerr << "**** getHist3D: <" << hname 
  	 << "> may not be a 3D Histogram! (" 
  	 << __FILE__ << ":" << __LINE__ << ")" 
  	 << endl;
    }
    return h;
  }
  TH3* getHist3D(const string& hname) {
    return getHist3D(hname.c_str());
  }
  
  // --------------------------------------------------
  // Convenience routine for filling profile histograms
  // --------------------------------------------------
  TProfile* getProfile(const char* hname) {
    TProfile *h = dynamic_cast<TProfile*>(gDirectory->GetList()->FindObject(hname));
    if (h == nullptr) {
      cerr << "**** getProfile: Profile Histogram for <" << hname 
  	 << "> not found! (" 
  	 << __FILE__ << ":" << __LINE__ << ")" 
  	 << endl;
    }
    return h;
  }
  TProfile* getProfile(const string& hname) {
    return getProfile(hname.c_str());
  }
  
  bool fillProfile(const char *hname, float xvalue, float yvalue, double w) {
    TProfile *h = getProfile(hname);
    if (h == nullptr) return false;
    
    h->Fill(xvalue, yvalue, w);
    return true;
  }
  bool fillProfile(const string& hname, float xvalue, float yvalue, double w) {
    return fillProfile(hname.c_str(), xvalue, yvalue, w);
  }
}
