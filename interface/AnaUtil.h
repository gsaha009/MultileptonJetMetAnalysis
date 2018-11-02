#ifndef __ANAUTIL__HH
#define __ANAUTIL__HH

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TH1K.h"

using std::ostream;

namespace AnaUtil {
  // Templated functioned must be defined in the header itself
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");
  void bit_print(int value, int pos=32, ostream& os=std::cout);
  template <typename T>
  T deltaPhiT(T phi1, T phi2) {
    T result = phi1 - phi2;
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return result;
  }  
  double deltaPhi(double phia, double phib);
  double deltaPhi(const TLorentzVector& a, const TLorentzVector& b);
  double deltaR(const TLorentzVector& a, const TLorentzVector& b);
  bool sameObject(const TLorentzVector& lv1, const TLorentzVector& lv2);
  double cutValue(const std::map<std::string, double>& m, const std::string& cname);
  const std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);
  void buildList(const std::vector<std::string>& tokens, std::vector<std::string>& list);
  void buildMap(const std::vector<std::string>& tokens, std::map<std::string, int>& hmap);
  void buildMap(const std::vector<std::string>& tokens, std::unordered_map<std::string, int>& hmap);
  void storeCuts(const std::vector<std::string>& tokens, std::map<std::string, std::map<std::string, double>>& hmap);
  void showCuts(const std::map<std::string, std::map<std::string, double> >& hmap, ostream& os=std::cout);
  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  TH1* getHist1D(const char* hname);
  TH1* getHist1D(const std::string& hname);
  template <class T>
  bool fillHist1D(const char* hname, T value, double w=1.0) {
    TH1* h = getHist1D(hname);
    if (h == nullptr) return false;
    h->Fill(value, w);
    return true;
  }
  template <class T>
    bool fillHist1D(const std::string& hname, T value, double w=1.0) {
    return fillHist1D(hname.c_str(), value, w);
  }

  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  TH2* getHist2D(const char* hname);
  TH2* getHist2D(const std::string& hname);
  template <class T1, class T2>
  bool fillHist2D(const char* hname, T1 xvalue, T2 yvalue, double w=1.0) {
    TH2* h = getHist2D(hname);
    if (h == nullptr) return false;
    h->Fill(xvalue, yvalue, w);
    return true;
  }
  template <class T1, class T2>
  bool fillHist2D(const std::string& hname, T1 xvalue, T2 yvalue, double w=1.0) {
    return fillHist2D(hname.c_str(), xvalue, yvalue, w);
  }
  // ---------------------------------------------
  // Convenience routine for filling 3D histograms
  // ---------------------------------------------
  TH3* getHist3D(const char* hname);
  TH3* getHist3D(const std::string& hname);
  template <class T1, class T2, class T3>
  bool fillHist3D(const char* hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) {
    TH3* h = getHist3D(hname);
    if (h == nullptr) return false;
    h->Fill(xvalue, yvalue, zvalue, w);
    return true;
  }
  template <class T1, class T2, class T3>
  bool fillHist3D(const std::string& hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) {
    return fillHist3D(hname.c_str(), xvalue, yvalue, zvalue, w);
  }

  // --------------------------------------------------
  // Convenience routine for filling profile histograms
  // --------------------------------------------------
  TProfile* getProfile(const char* hname);
  TProfile* getProfile(const std::string& hname);
  bool fillProfile(const char *hname, float xvalue, float yvalue, double w=1.0);
  bool fillProfile(const std::string& hname, float xvalue, float yvalue, double w=1.0);

  /* print_list_elements()
   * - prints optional C-string optcstr followed by
   * - all elements of the collection coll
   */
  template <class T>
  void showList(const T& coll, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << ", Total # = " << coll.size() << ":" << std::endl;
    for (auto const& v: coll)
      os << v << std::endl;
  }  
  template <class T1, class T2>
  void showMap(const std::map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << std::endl;
    for (auto const& k: m)
      os << k.first << std::endl;
  }
  template <class T1, class T2>
  void showMap(const std::unordered_map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << std::endl;
    for (auto const& k: m)
      os << k.first << std::endl;
  }
  /* copyList */
  template <class T>
  void copyList (const T& sourceColl, T& destColl) {
    destColl.clear();
    for (auto const& v: sourceColl)   
      destColl.push_back(v); 
  }
}
#endif
