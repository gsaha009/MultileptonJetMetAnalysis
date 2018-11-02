// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "interface/PhysicsObjects.h"
#include "interface/PhysicsObjects.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *pairlEstringcOunsignedsPintgR_Dictionary();
   static void pairlEstringcOunsignedsPintgR_TClassManip(TClass*);
   static void *new_pairlEstringcOunsignedsPintgR(void *p = 0);
   static void *newArray_pairlEstringcOunsignedsPintgR(Long_t size, void *p);
   static void delete_pairlEstringcOunsignedsPintgR(void *p);
   static void deleteArray_pairlEstringcOunsignedsPintgR(void *p);
   static void destruct_pairlEstringcOunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,unsigned int>*)
   {
      pair<string,unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,unsigned int>", "string", 190,
                  typeid(pair<string,unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOunsignedsPintgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,unsigned int>) );
      instance.SetNew(&new_pairlEstringcOunsignedsPintgR);
      instance.SetNewArray(&newArray_pairlEstringcOunsignedsPintgR);
      instance.SetDelete(&delete_pairlEstringcOunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOunsignedsPintgR);
      instance.SetDestructor(&destruct_pairlEstringcOunsignedsPintgR);

      ::ROOT::AddClassAlternate("pair<string,unsigned int>","pair<std::string,unsigned int>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const pair<string,unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const pair<string,unsigned int>*)0x0)->GetClass();
      pairlEstringcOunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEvhtmcLcLMuoncOvhtmcLcLMuongR_Dictionary();
   static void pairlEvhtmcLcLMuoncOvhtmcLcLMuongR_TClassManip(TClass*);
   static void *new_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p = 0);
   static void *newArray_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(Long_t size, void *p);
   static void delete_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p);
   static void deleteArray_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p);
   static void destruct_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<vhtm::Muon,vhtm::Muon>*)
   {
      pair<vhtm::Muon,vhtm::Muon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<vhtm::Muon,vhtm::Muon>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<vhtm::Muon,vhtm::Muon>", "string", 190,
                  typeid(pair<vhtm::Muon,vhtm::Muon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEvhtmcLcLMuoncOvhtmcLcLMuongR_Dictionary, isa_proxy, 4,
                  sizeof(pair<vhtm::Muon,vhtm::Muon>) );
      instance.SetNew(&new_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR);
      instance.SetNewArray(&newArray_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR);
      instance.SetDelete(&delete_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR);
      instance.SetDeleteArray(&deleteArray_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR);
      instance.SetDestructor(&destruct_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR);
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const pair<vhtm::Muon,vhtm::Muon>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEvhtmcLcLMuoncOvhtmcLcLMuongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const pair<vhtm::Muon,vhtm::Muon>*)0x0)->GetClass();
      pairlEvhtmcLcLMuoncOvhtmcLcLMuongR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEvhtmcLcLMuoncOvhtmcLcLMuongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLCandidate(void *p = 0);
   static void *newArray_vhtmcLcLCandidate(Long_t size, void *p);
   static void delete_vhtmcLcLCandidate(void *p);
   static void deleteArray_vhtmcLcLCandidate(void *p);
   static void destruct_vhtmcLcLCandidate(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Candidate*)
   {
      ::vhtm::Candidate *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Candidate >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Candidate", ::vhtm::Candidate::Class_Version(), "PhysicsObjects.h", 11,
                  typeid(::vhtm::Candidate), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Candidate::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Candidate) );
      instance.SetNew(&new_vhtmcLcLCandidate);
      instance.SetNewArray(&newArray_vhtmcLcLCandidate);
      instance.SetDelete(&delete_vhtmcLcLCandidate);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLCandidate);
      instance.SetDestructor(&destruct_vhtmcLcLCandidate);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Candidate*)
   {
      return GenerateInitInstanceLocal((::vhtm::Candidate*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Candidate*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLPackedPFCandidate(void *p = 0);
   static void *newArray_vhtmcLcLPackedPFCandidate(Long_t size, void *p);
   static void delete_vhtmcLcLPackedPFCandidate(void *p);
   static void deleteArray_vhtmcLcLPackedPFCandidate(void *p);
   static void destruct_vhtmcLcLPackedPFCandidate(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::PackedPFCandidate*)
   {
      ::vhtm::PackedPFCandidate *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::PackedPFCandidate >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::PackedPFCandidate", ::vhtm::PackedPFCandidate::Class_Version(), "PhysicsObjects.h", 24,
                  typeid(::vhtm::PackedPFCandidate), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::PackedPFCandidate::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::PackedPFCandidate) );
      instance.SetNew(&new_vhtmcLcLPackedPFCandidate);
      instance.SetNewArray(&newArray_vhtmcLcLPackedPFCandidate);
      instance.SetDelete(&delete_vhtmcLcLPackedPFCandidate);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLPackedPFCandidate);
      instance.SetDestructor(&destruct_vhtmcLcLPackedPFCandidate);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::PackedPFCandidate*)
   {
      return GenerateInitInstanceLocal((::vhtm::PackedPFCandidate*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::PackedPFCandidate*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLEvent(void *p = 0);
   static void *newArray_vhtmcLcLEvent(Long_t size, void *p);
   static void delete_vhtmcLcLEvent(void *p);
   static void deleteArray_vhtmcLcLEvent(void *p);
   static void destruct_vhtmcLcLEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Event*)
   {
      ::vhtm::Event *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Event >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Event", ::vhtm::Event::Class_Version(), "PhysicsObjects.h", 66,
                  typeid(::vhtm::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Event::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Event) );
      instance.SetNew(&new_vhtmcLcLEvent);
      instance.SetNewArray(&newArray_vhtmcLcLEvent);
      instance.SetDelete(&delete_vhtmcLcLEvent);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLEvent);
      instance.SetDestructor(&destruct_vhtmcLcLEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Event*)
   {
      return GenerateInitInstanceLocal((::vhtm::Event*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Event*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLGenEvent(void *p = 0);
   static void *newArray_vhtmcLcLGenEvent(Long_t size, void *p);
   static void delete_vhtmcLcLGenEvent(void *p);
   static void deleteArray_vhtmcLcLGenEvent(void *p);
   static void destruct_vhtmcLcLGenEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::GenEvent*)
   {
      ::vhtm::GenEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::GenEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::GenEvent", ::vhtm::GenEvent::Class_Version(), "PhysicsObjects.h", 98,
                  typeid(::vhtm::GenEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::GenEvent::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::GenEvent) );
      instance.SetNew(&new_vhtmcLcLGenEvent);
      instance.SetNewArray(&newArray_vhtmcLcLGenEvent);
      instance.SetDelete(&delete_vhtmcLcLGenEvent);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLGenEvent);
      instance.SetDestructor(&destruct_vhtmcLcLGenEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::GenEvent*)
   {
      return GenerateInitInstanceLocal((::vhtm::GenEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::GenEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLElectron(void *p = 0);
   static void *newArray_vhtmcLcLElectron(Long_t size, void *p);
   static void delete_vhtmcLcLElectron(void *p);
   static void deleteArray_vhtmcLcLElectron(void *p);
   static void destruct_vhtmcLcLElectron(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Electron*)
   {
      ::vhtm::Electron *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Electron >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Electron", ::vhtm::Electron::Class_Version(), "PhysicsObjects.h", 129,
                  typeid(::vhtm::Electron), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Electron::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Electron) );
      instance.SetNew(&new_vhtmcLcLElectron);
      instance.SetNewArray(&newArray_vhtmcLcLElectron);
      instance.SetDelete(&delete_vhtmcLcLElectron);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLElectron);
      instance.SetDestructor(&destruct_vhtmcLcLElectron);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Electron*)
   {
      return GenerateInitInstanceLocal((::vhtm::Electron*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Electron*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLGenParticle(void *p = 0);
   static void *newArray_vhtmcLcLGenParticle(Long_t size, void *p);
   static void delete_vhtmcLcLGenParticle(void *p);
   static void deleteArray_vhtmcLcLGenParticle(void *p);
   static void destruct_vhtmcLcLGenParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::GenParticle*)
   {
      ::vhtm::GenParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::GenParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::GenParticle", ::vhtm::GenParticle::Class_Version(), "PhysicsObjects.h", 214,
                  typeid(::vhtm::GenParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::GenParticle::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::GenParticle) );
      instance.SetNew(&new_vhtmcLcLGenParticle);
      instance.SetNewArray(&newArray_vhtmcLcLGenParticle);
      instance.SetDelete(&delete_vhtmcLcLGenParticle);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLGenParticle);
      instance.SetDestructor(&destruct_vhtmcLcLGenParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::GenParticle*)
   {
      return GenerateInitInstanceLocal((::vhtm::GenParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::GenParticle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLGenJet(void *p = 0);
   static void *newArray_vhtmcLcLGenJet(Long_t size, void *p);
   static void delete_vhtmcLcLGenJet(void *p);
   static void deleteArray_vhtmcLcLGenJet(void *p);
   static void destruct_vhtmcLcLGenJet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::GenJet*)
   {
      ::vhtm::GenJet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::GenJet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::GenJet", ::vhtm::GenJet::Class_Version(), "PhysicsObjects.h", 241,
                  typeid(::vhtm::GenJet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::GenJet::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::GenJet) );
      instance.SetNew(&new_vhtmcLcLGenJet);
      instance.SetNewArray(&newArray_vhtmcLcLGenJet);
      instance.SetDelete(&delete_vhtmcLcLGenJet);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLGenJet);
      instance.SetDestructor(&destruct_vhtmcLcLGenJet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::GenJet*)
   {
      return GenerateInitInstanceLocal((::vhtm::GenJet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::GenJet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLMET(void *p = 0);
   static void *newArray_vhtmcLcLMET(Long_t size, void *p);
   static void delete_vhtmcLcLMET(void *p);
   static void deleteArray_vhtmcLcLMET(void *p);
   static void destruct_vhtmcLcLMET(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::MET*)
   {
      ::vhtm::MET *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::MET >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::MET", ::vhtm::MET::Class_Version(), "PhysicsObjects.h", 256,
                  typeid(::vhtm::MET), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::MET::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::MET) );
      instance.SetNew(&new_vhtmcLcLMET);
      instance.SetNewArray(&newArray_vhtmcLcLMET);
      instance.SetDelete(&delete_vhtmcLcLMET);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLMET);
      instance.SetDestructor(&destruct_vhtmcLcLMET);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::MET*)
   {
      return GenerateInitInstanceLocal((::vhtm::MET*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::MET*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLTau(void *p = 0);
   static void *newArray_vhtmcLcLTau(Long_t size, void *p);
   static void delete_vhtmcLcLTau(void *p);
   static void deleteArray_vhtmcLcLTau(void *p);
   static void destruct_vhtmcLcLTau(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Tau*)
   {
      ::vhtm::Tau *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Tau >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Tau", ::vhtm::Tau::Class_Version(), "PhysicsObjects.h", 274,
                  typeid(::vhtm::Tau), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Tau::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Tau) );
      instance.SetNew(&new_vhtmcLcLTau);
      instance.SetNewArray(&newArray_vhtmcLcLTau);
      instance.SetDelete(&delete_vhtmcLcLTau);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLTau);
      instance.SetDestructor(&destruct_vhtmcLcLTau);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Tau*)
   {
      return GenerateInitInstanceLocal((::vhtm::Tau*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Tau*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLMuon(void *p = 0);
   static void *newArray_vhtmcLcLMuon(Long_t size, void *p);
   static void delete_vhtmcLcLMuon(void *p);
   static void deleteArray_vhtmcLcLMuon(void *p);
   static void destruct_vhtmcLcLMuon(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Muon*)
   {
      ::vhtm::Muon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Muon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Muon", ::vhtm::Muon::Class_Version(), "PhysicsObjects.h", 355,
                  typeid(::vhtm::Muon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Muon::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Muon) );
      instance.SetNew(&new_vhtmcLcLMuon);
      instance.SetNewArray(&newArray_vhtmcLcLMuon);
      instance.SetDelete(&delete_vhtmcLcLMuon);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLMuon);
      instance.SetDestructor(&destruct_vhtmcLcLMuon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Muon*)
   {
      return GenerateInitInstanceLocal((::vhtm::Muon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Muon*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLJet(void *p = 0);
   static void *newArray_vhtmcLcLJet(Long_t size, void *p);
   static void delete_vhtmcLcLJet(void *p);
   static void deleteArray_vhtmcLcLJet(void *p);
   static void destruct_vhtmcLcLJet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Jet*)
   {
      ::vhtm::Jet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Jet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Jet", ::vhtm::Jet::Class_Version(), "PhysicsObjects.h", 441,
                  typeid(::vhtm::Jet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Jet::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Jet) );
      instance.SetNew(&new_vhtmcLcLJet);
      instance.SetNewArray(&newArray_vhtmcLcLJet);
      instance.SetDelete(&delete_vhtmcLcLJet);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLJet);
      instance.SetDestructor(&destruct_vhtmcLcLJet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Jet*)
   {
      return GenerateInitInstanceLocal((::vhtm::Jet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Jet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLVertex(void *p = 0);
   static void *newArray_vhtmcLcLVertex(Long_t size, void *p);
   static void delete_vhtmcLcLVertex(void *p);
   static void deleteArray_vhtmcLcLVertex(void *p);
   static void destruct_vhtmcLcLVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Vertex*)
   {
      ::vhtm::Vertex *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Vertex >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Vertex", ::vhtm::Vertex::Class_Version(), "PhysicsObjects.h", 486,
                  typeid(::vhtm::Vertex), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Vertex::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Vertex) );
      instance.SetNew(&new_vhtmcLcLVertex);
      instance.SetNewArray(&newArray_vhtmcLcLVertex);
      instance.SetDelete(&delete_vhtmcLcLVertex);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLVertex);
      instance.SetDestructor(&destruct_vhtmcLcLVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Vertex*)
   {
      return GenerateInitInstanceLocal((::vhtm::Vertex*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Vertex*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLGenMET(void *p = 0);
   static void *newArray_vhtmcLcLGenMET(Long_t size, void *p);
   static void delete_vhtmcLcLGenMET(void *p);
   static void deleteArray_vhtmcLcLGenMET(void *p);
   static void destruct_vhtmcLcLGenMET(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::GenMET*)
   {
      ::vhtm::GenMET *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::GenMET >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::GenMET", ::vhtm::GenMET::Class_Version(), "PhysicsObjects.h", 505,
                  typeid(::vhtm::GenMET), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::GenMET::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::GenMET) );
      instance.SetNew(&new_vhtmcLcLGenMET);
      instance.SetNewArray(&newArray_vhtmcLcLGenMET);
      instance.SetDelete(&delete_vhtmcLcLGenMET);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLGenMET);
      instance.SetDestructor(&destruct_vhtmcLcLGenMET);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::GenMET*)
   {
      return GenerateInitInstanceLocal((::vhtm::GenMET*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::GenMET*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLTriggerObject(void *p = 0);
   static void *newArray_vhtmcLcLTriggerObject(Long_t size, void *p);
   static void delete_vhtmcLcLTriggerObject(void *p);
   static void deleteArray_vhtmcLcLTriggerObject(void *p);
   static void destruct_vhtmcLcLTriggerObject(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::TriggerObject*)
   {
      ::vhtm::TriggerObject *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::TriggerObject >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::TriggerObject", ::vhtm::TriggerObject::Class_Version(), "PhysicsObjects.h", 516,
                  typeid(::vhtm::TriggerObject), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::TriggerObject::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::TriggerObject) );
      instance.SetNew(&new_vhtmcLcLTriggerObject);
      instance.SetNewArray(&newArray_vhtmcLcLTriggerObject);
      instance.SetDelete(&delete_vhtmcLcLTriggerObject);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLTriggerObject);
      instance.SetDestructor(&destruct_vhtmcLcLTriggerObject);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::TriggerObject*)
   {
      return GenerateInitInstanceLocal((::vhtm::TriggerObject*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::TriggerObject*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_vhtmcLcLPhoton(void *p = 0);
   static void *newArray_vhtmcLcLPhoton(Long_t size, void *p);
   static void delete_vhtmcLcLPhoton(void *p);
   static void deleteArray_vhtmcLcLPhoton(void *p);
   static void destruct_vhtmcLcLPhoton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::vhtm::Photon*)
   {
      ::vhtm::Photon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::vhtm::Photon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("vhtm::Photon", ::vhtm::Photon::Class_Version(), "PhysicsObjects.h", 530,
                  typeid(::vhtm::Photon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::vhtm::Photon::Dictionary, isa_proxy, 4,
                  sizeof(::vhtm::Photon) );
      instance.SetNew(&new_vhtmcLcLPhoton);
      instance.SetNewArray(&newArray_vhtmcLcLPhoton);
      instance.SetDelete(&delete_vhtmcLcLPhoton);
      instance.SetDeleteArray(&deleteArray_vhtmcLcLPhoton);
      instance.SetDestructor(&destruct_vhtmcLcLPhoton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::vhtm::Photon*)
   {
      return GenerateInitInstanceLocal((::vhtm::Photon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::vhtm::Photon*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Candidate::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Candidate::Class_Name()
{
   return "vhtm::Candidate";
}

//______________________________________________________________________________
const char *Candidate::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Candidate*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Candidate::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Candidate*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Candidate::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Candidate*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Candidate::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Candidate*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr PackedPFCandidate::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PackedPFCandidate::Class_Name()
{
   return "vhtm::PackedPFCandidate";
}

//______________________________________________________________________________
const char *PackedPFCandidate::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::PackedPFCandidate*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PackedPFCandidate::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::PackedPFCandidate*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PackedPFCandidate::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::PackedPFCandidate*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PackedPFCandidate::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::PackedPFCandidate*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "vhtm::Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Event*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Event*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Event*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Event*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr GenEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenEvent::Class_Name()
{
   return "vhtm::GenEvent";
}

//______________________________________________________________________________
const char *GenEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenEvent*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Electron::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Electron::Class_Name()
{
   return "vhtm::Electron";
}

//______________________________________________________________________________
const char *Electron::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Electron*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Electron::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Electron*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Electron::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Electron*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Electron::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Electron*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr GenParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenParticle::Class_Name()
{
   return "vhtm::GenParticle";
}

//______________________________________________________________________________
const char *GenParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenParticle*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr GenJet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenJet::Class_Name()
{
   return "vhtm::GenJet";
}

//______________________________________________________________________________
const char *GenJet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenJet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenJet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenJet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenJet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenJet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenJet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenJet*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr MET::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MET::Class_Name()
{
   return "vhtm::MET";
}

//______________________________________________________________________________
const char *MET::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::MET*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MET::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::MET*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MET::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::MET*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MET::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::MET*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Tau::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Tau::Class_Name()
{
   return "vhtm::Tau";
}

//______________________________________________________________________________
const char *Tau::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Tau*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Tau::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Tau*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Tau::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Tau*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Tau::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Tau*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Muon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Muon::Class_Name()
{
   return "vhtm::Muon";
}

//______________________________________________________________________________
const char *Muon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Muon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Muon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Muon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Muon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Muon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Muon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Muon*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Jet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Jet::Class_Name()
{
   return "vhtm::Jet";
}

//______________________________________________________________________________
const char *Jet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Jet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Jet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Jet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Jet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Jet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Jet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Jet*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Vertex::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Vertex::Class_Name()
{
   return "vhtm::Vertex";
}

//______________________________________________________________________________
const char *Vertex::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Vertex*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Vertex::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Vertex*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Vertex::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Vertex*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Vertex::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Vertex*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr GenMET::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenMET::Class_Name()
{
   return "vhtm::GenMET";
}

//______________________________________________________________________________
const char *GenMET::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenMET*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenMET::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenMET*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenMET::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenMET*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenMET::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::GenMET*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr TriggerObject::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TriggerObject::Class_Name()
{
   return "vhtm::TriggerObject";
}

//______________________________________________________________________________
const char *TriggerObject::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::TriggerObject*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TriggerObject::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::TriggerObject*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TriggerObject::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::TriggerObject*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TriggerObject::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::TriggerObject*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace vhtm {
//______________________________________________________________________________
atomic_TClass_ptr Photon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Photon::Class_Name()
{
   return "vhtm::Photon";
}

//______________________________________________________________________________
const char *Photon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Photon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Photon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Photon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Photon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Photon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Photon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::vhtm::Photon*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,unsigned int> : new pair<string,unsigned int>;
   }
   static void *newArray_pairlEstringcOunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,unsigned int>[nElements] : new pair<string,unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOunsignedsPintgR(void *p) {
      delete ((pair<string,unsigned int>*)p);
   }
   static void deleteArray_pairlEstringcOunsignedsPintgR(void *p) {
      delete [] ((pair<string,unsigned int>*)p);
   }
   static void destruct_pairlEstringcOunsignedsPintgR(void *p) {
      typedef pair<string,unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class pair<string,unsigned int>

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<vhtm::Muon,vhtm::Muon> : new pair<vhtm::Muon,vhtm::Muon>;
   }
   static void *newArray_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<vhtm::Muon,vhtm::Muon>[nElements] : new pair<vhtm::Muon,vhtm::Muon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p) {
      delete ((pair<vhtm::Muon,vhtm::Muon>*)p);
   }
   static void deleteArray_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p) {
      delete [] ((pair<vhtm::Muon,vhtm::Muon>*)p);
   }
   static void destruct_pairlEvhtmcLcLMuoncOvhtmcLcLMuongR(void *p) {
      typedef pair<vhtm::Muon,vhtm::Muon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class pair<vhtm::Muon,vhtm::Muon>

namespace vhtm {
//______________________________________________________________________________
void Candidate::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Candidate.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Candidate::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Candidate::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLCandidate(void *p) {
      return  p ? new(p) ::vhtm::Candidate : new ::vhtm::Candidate;
   }
   static void *newArray_vhtmcLcLCandidate(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Candidate[nElements] : new ::vhtm::Candidate[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLCandidate(void *p) {
      delete ((::vhtm::Candidate*)p);
   }
   static void deleteArray_vhtmcLcLCandidate(void *p) {
      delete [] ((::vhtm::Candidate*)p);
   }
   static void destruct_vhtmcLcLCandidate(void *p) {
      typedef ::vhtm::Candidate current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Candidate

namespace vhtm {
//______________________________________________________________________________
void PackedPFCandidate::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::PackedPFCandidate.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::PackedPFCandidate::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::PackedPFCandidate::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLPackedPFCandidate(void *p) {
      return  p ? new(p) ::vhtm::PackedPFCandidate : new ::vhtm::PackedPFCandidate;
   }
   static void *newArray_vhtmcLcLPackedPFCandidate(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::PackedPFCandidate[nElements] : new ::vhtm::PackedPFCandidate[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLPackedPFCandidate(void *p) {
      delete ((::vhtm::PackedPFCandidate*)p);
   }
   static void deleteArray_vhtmcLcLPackedPFCandidate(void *p) {
      delete [] ((::vhtm::PackedPFCandidate*)p);
   }
   static void destruct_vhtmcLcLPackedPFCandidate(void *p) {
      typedef ::vhtm::PackedPFCandidate current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::PackedPFCandidate

namespace vhtm {
//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Event::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLEvent(void *p) {
      return  p ? new(p) ::vhtm::Event : new ::vhtm::Event;
   }
   static void *newArray_vhtmcLcLEvent(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Event[nElements] : new ::vhtm::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLEvent(void *p) {
      delete ((::vhtm::Event*)p);
   }
   static void deleteArray_vhtmcLcLEvent(void *p) {
      delete [] ((::vhtm::Event*)p);
   }
   static void destruct_vhtmcLcLEvent(void *p) {
      typedef ::vhtm::Event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Event

namespace vhtm {
//______________________________________________________________________________
void GenEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::GenEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::GenEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::GenEvent::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLGenEvent(void *p) {
      return  p ? new(p) ::vhtm::GenEvent : new ::vhtm::GenEvent;
   }
   static void *newArray_vhtmcLcLGenEvent(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::GenEvent[nElements] : new ::vhtm::GenEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLGenEvent(void *p) {
      delete ((::vhtm::GenEvent*)p);
   }
   static void deleteArray_vhtmcLcLGenEvent(void *p) {
      delete [] ((::vhtm::GenEvent*)p);
   }
   static void destruct_vhtmcLcLGenEvent(void *p) {
      typedef ::vhtm::GenEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::GenEvent

namespace vhtm {
//______________________________________________________________________________
void Electron::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Electron.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Electron::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Electron::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLElectron(void *p) {
      return  p ? new(p) ::vhtm::Electron : new ::vhtm::Electron;
   }
   static void *newArray_vhtmcLcLElectron(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Electron[nElements] : new ::vhtm::Electron[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLElectron(void *p) {
      delete ((::vhtm::Electron*)p);
   }
   static void deleteArray_vhtmcLcLElectron(void *p) {
      delete [] ((::vhtm::Electron*)p);
   }
   static void destruct_vhtmcLcLElectron(void *p) {
      typedef ::vhtm::Electron current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Electron

namespace vhtm {
//______________________________________________________________________________
void GenParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::GenParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::GenParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::GenParticle::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLGenParticle(void *p) {
      return  p ? new(p) ::vhtm::GenParticle : new ::vhtm::GenParticle;
   }
   static void *newArray_vhtmcLcLGenParticle(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::GenParticle[nElements] : new ::vhtm::GenParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLGenParticle(void *p) {
      delete ((::vhtm::GenParticle*)p);
   }
   static void deleteArray_vhtmcLcLGenParticle(void *p) {
      delete [] ((::vhtm::GenParticle*)p);
   }
   static void destruct_vhtmcLcLGenParticle(void *p) {
      typedef ::vhtm::GenParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::GenParticle

namespace vhtm {
//______________________________________________________________________________
void GenJet::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::GenJet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::GenJet::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::GenJet::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLGenJet(void *p) {
      return  p ? new(p) ::vhtm::GenJet : new ::vhtm::GenJet;
   }
   static void *newArray_vhtmcLcLGenJet(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::GenJet[nElements] : new ::vhtm::GenJet[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLGenJet(void *p) {
      delete ((::vhtm::GenJet*)p);
   }
   static void deleteArray_vhtmcLcLGenJet(void *p) {
      delete [] ((::vhtm::GenJet*)p);
   }
   static void destruct_vhtmcLcLGenJet(void *p) {
      typedef ::vhtm::GenJet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::GenJet

namespace vhtm {
//______________________________________________________________________________
void MET::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::MET.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::MET::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::MET::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLMET(void *p) {
      return  p ? new(p) ::vhtm::MET : new ::vhtm::MET;
   }
   static void *newArray_vhtmcLcLMET(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::MET[nElements] : new ::vhtm::MET[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLMET(void *p) {
      delete ((::vhtm::MET*)p);
   }
   static void deleteArray_vhtmcLcLMET(void *p) {
      delete [] ((::vhtm::MET*)p);
   }
   static void destruct_vhtmcLcLMET(void *p) {
      typedef ::vhtm::MET current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::MET

namespace vhtm {
//______________________________________________________________________________
void Tau::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Tau.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Tau::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Tau::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLTau(void *p) {
      return  p ? new(p) ::vhtm::Tau : new ::vhtm::Tau;
   }
   static void *newArray_vhtmcLcLTau(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Tau[nElements] : new ::vhtm::Tau[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLTau(void *p) {
      delete ((::vhtm::Tau*)p);
   }
   static void deleteArray_vhtmcLcLTau(void *p) {
      delete [] ((::vhtm::Tau*)p);
   }
   static void destruct_vhtmcLcLTau(void *p) {
      typedef ::vhtm::Tau current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Tau

namespace vhtm {
//______________________________________________________________________________
void Muon::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Muon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Muon::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Muon::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLMuon(void *p) {
      return  p ? new(p) ::vhtm::Muon : new ::vhtm::Muon;
   }
   static void *newArray_vhtmcLcLMuon(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Muon[nElements] : new ::vhtm::Muon[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLMuon(void *p) {
      delete ((::vhtm::Muon*)p);
   }
   static void deleteArray_vhtmcLcLMuon(void *p) {
      delete [] ((::vhtm::Muon*)p);
   }
   static void destruct_vhtmcLcLMuon(void *p) {
      typedef ::vhtm::Muon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Muon

namespace vhtm {
//______________________________________________________________________________
void Jet::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Jet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Jet::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Jet::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLJet(void *p) {
      return  p ? new(p) ::vhtm::Jet : new ::vhtm::Jet;
   }
   static void *newArray_vhtmcLcLJet(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Jet[nElements] : new ::vhtm::Jet[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLJet(void *p) {
      delete ((::vhtm::Jet*)p);
   }
   static void deleteArray_vhtmcLcLJet(void *p) {
      delete [] ((::vhtm::Jet*)p);
   }
   static void destruct_vhtmcLcLJet(void *p) {
      typedef ::vhtm::Jet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Jet

namespace vhtm {
//______________________________________________________________________________
void Vertex::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Vertex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Vertex::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Vertex::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLVertex(void *p) {
      return  p ? new(p) ::vhtm::Vertex : new ::vhtm::Vertex;
   }
   static void *newArray_vhtmcLcLVertex(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Vertex[nElements] : new ::vhtm::Vertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLVertex(void *p) {
      delete ((::vhtm::Vertex*)p);
   }
   static void deleteArray_vhtmcLcLVertex(void *p) {
      delete [] ((::vhtm::Vertex*)p);
   }
   static void destruct_vhtmcLcLVertex(void *p) {
      typedef ::vhtm::Vertex current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Vertex

namespace vhtm {
//______________________________________________________________________________
void GenMET::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::GenMET.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::GenMET::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::GenMET::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLGenMET(void *p) {
      return  p ? new(p) ::vhtm::GenMET : new ::vhtm::GenMET;
   }
   static void *newArray_vhtmcLcLGenMET(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::GenMET[nElements] : new ::vhtm::GenMET[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLGenMET(void *p) {
      delete ((::vhtm::GenMET*)p);
   }
   static void deleteArray_vhtmcLcLGenMET(void *p) {
      delete [] ((::vhtm::GenMET*)p);
   }
   static void destruct_vhtmcLcLGenMET(void *p) {
      typedef ::vhtm::GenMET current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::GenMET

namespace vhtm {
//______________________________________________________________________________
void TriggerObject::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::TriggerObject.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::TriggerObject::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::TriggerObject::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLTriggerObject(void *p) {
      return  p ? new(p) ::vhtm::TriggerObject : new ::vhtm::TriggerObject;
   }
   static void *newArray_vhtmcLcLTriggerObject(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::TriggerObject[nElements] : new ::vhtm::TriggerObject[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLTriggerObject(void *p) {
      delete ((::vhtm::TriggerObject*)p);
   }
   static void deleteArray_vhtmcLcLTriggerObject(void *p) {
      delete [] ((::vhtm::TriggerObject*)p);
   }
   static void destruct_vhtmcLcLTriggerObject(void *p) {
      typedef ::vhtm::TriggerObject current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::TriggerObject

namespace vhtm {
//______________________________________________________________________________
void Photon::Streamer(TBuffer &R__b)
{
   // Stream an object of class vhtm::Photon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(vhtm::Photon::Class(),this);
   } else {
      R__b.WriteClassBuffer(vhtm::Photon::Class(),this);
   }
}

} // namespace vhtm
namespace ROOT {
   // Wrappers around operator new
   static void *new_vhtmcLcLPhoton(void *p) {
      return  p ? new(p) ::vhtm::Photon : new ::vhtm::Photon;
   }
   static void *newArray_vhtmcLcLPhoton(Long_t nElements, void *p) {
      return p ? new(p) ::vhtm::Photon[nElements] : new ::vhtm::Photon[nElements];
   }
   // Wrapper around operator delete
   static void delete_vhtmcLcLPhoton(void *p) {
      delete ((::vhtm::Photon*)p);
   }
   static void deleteArray_vhtmcLcLPhoton(void *p) {
      delete [] ((::vhtm::Photon*)p);
   }
   static void destruct_vhtmcLcLPhoton(void *p) {
      typedef ::vhtm::Photon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::vhtm::Photon

namespace ROOT {
   static TClass *vectorlEvhtmcLcLVertexgR_Dictionary();
   static void vectorlEvhtmcLcLVertexgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLVertexgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLVertexgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLVertexgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLVertexgR(void *p);
   static void destruct_vectorlEvhtmcLcLVertexgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Vertex>*)
   {
      vector<vhtm::Vertex> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Vertex>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Vertex>", -2, "vector", 214,
                  typeid(vector<vhtm::Vertex>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLVertexgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Vertex>) );
      instance.SetNew(&new_vectorlEvhtmcLcLVertexgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLVertexgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLVertexgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLVertexgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLVertexgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Vertex> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Vertex>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLVertexgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Vertex>*)0x0)->GetClass();
      vectorlEvhtmcLcLVertexgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLVertexgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLVertexgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Vertex> : new vector<vhtm::Vertex>;
   }
   static void *newArray_vectorlEvhtmcLcLVertexgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Vertex>[nElements] : new vector<vhtm::Vertex>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLVertexgR(void *p) {
      delete ((vector<vhtm::Vertex>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLVertexgR(void *p) {
      delete [] ((vector<vhtm::Vertex>*)p);
   }
   static void destruct_vectorlEvhtmcLcLVertexgR(void *p) {
      typedef vector<vhtm::Vertex> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Vertex>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLTriggerObjectgR_Dictionary();
   static void vectorlEvhtmcLcLTriggerObjectgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLTriggerObjectgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLTriggerObjectgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLTriggerObjectgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLTriggerObjectgR(void *p);
   static void destruct_vectorlEvhtmcLcLTriggerObjectgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::TriggerObject>*)
   {
      vector<vhtm::TriggerObject> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::TriggerObject>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::TriggerObject>", -2, "vector", 214,
                  typeid(vector<vhtm::TriggerObject>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLTriggerObjectgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::TriggerObject>) );
      instance.SetNew(&new_vectorlEvhtmcLcLTriggerObjectgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLTriggerObjectgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLTriggerObjectgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLTriggerObjectgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLTriggerObjectgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::TriggerObject> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::TriggerObject>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLTriggerObjectgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::TriggerObject>*)0x0)->GetClass();
      vectorlEvhtmcLcLTriggerObjectgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLTriggerObjectgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLTriggerObjectgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::TriggerObject> : new vector<vhtm::TriggerObject>;
   }
   static void *newArray_vectorlEvhtmcLcLTriggerObjectgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::TriggerObject>[nElements] : new vector<vhtm::TriggerObject>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLTriggerObjectgR(void *p) {
      delete ((vector<vhtm::TriggerObject>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLTriggerObjectgR(void *p) {
      delete [] ((vector<vhtm::TriggerObject>*)p);
   }
   static void destruct_vectorlEvhtmcLcLTriggerObjectgR(void *p) {
      typedef vector<vhtm::TriggerObject> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::TriggerObject>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLTaugR_Dictionary();
   static void vectorlEvhtmcLcLTaugR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLTaugR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLTaugR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLTaugR(void *p);
   static void deleteArray_vectorlEvhtmcLcLTaugR(void *p);
   static void destruct_vectorlEvhtmcLcLTaugR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Tau>*)
   {
      vector<vhtm::Tau> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Tau>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Tau>", -2, "vector", 214,
                  typeid(vector<vhtm::Tau>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLTaugR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Tau>) );
      instance.SetNew(&new_vectorlEvhtmcLcLTaugR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLTaugR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLTaugR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLTaugR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLTaugR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Tau> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Tau>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLTaugR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Tau>*)0x0)->GetClass();
      vectorlEvhtmcLcLTaugR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLTaugR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLTaugR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Tau> : new vector<vhtm::Tau>;
   }
   static void *newArray_vectorlEvhtmcLcLTaugR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Tau>[nElements] : new vector<vhtm::Tau>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLTaugR(void *p) {
      delete ((vector<vhtm::Tau>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLTaugR(void *p) {
      delete [] ((vector<vhtm::Tau>*)p);
   }
   static void destruct_vectorlEvhtmcLcLTaugR(void *p) {
      typedef vector<vhtm::Tau> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Tau>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLPhotongR_Dictionary();
   static void vectorlEvhtmcLcLPhotongR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLPhotongR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLPhotongR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLPhotongR(void *p);
   static void deleteArray_vectorlEvhtmcLcLPhotongR(void *p);
   static void destruct_vectorlEvhtmcLcLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Photon>*)
   {
      vector<vhtm::Photon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Photon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Photon>", -2, "vector", 214,
                  typeid(vector<vhtm::Photon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Photon>) );
      instance.SetNew(&new_vectorlEvhtmcLcLPhotongR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLPhotongR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLPhotongR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLPhotongR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLPhotongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Photon> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Photon>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Photon>*)0x0)->GetClass();
      vectorlEvhtmcLcLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLPhotongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Photon> : new vector<vhtm::Photon>;
   }
   static void *newArray_vectorlEvhtmcLcLPhotongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Photon>[nElements] : new vector<vhtm::Photon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLPhotongR(void *p) {
      delete ((vector<vhtm::Photon>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLPhotongR(void *p) {
      delete [] ((vector<vhtm::Photon>*)p);
   }
   static void destruct_vectorlEvhtmcLcLPhotongR(void *p) {
      typedef vector<vhtm::Photon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Photon>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLPackedPFCandidategR_Dictionary();
   static void vectorlEvhtmcLcLPackedPFCandidategR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLPackedPFCandidategR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLPackedPFCandidategR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLPackedPFCandidategR(void *p);
   static void deleteArray_vectorlEvhtmcLcLPackedPFCandidategR(void *p);
   static void destruct_vectorlEvhtmcLcLPackedPFCandidategR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::PackedPFCandidate>*)
   {
      vector<vhtm::PackedPFCandidate> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::PackedPFCandidate>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::PackedPFCandidate>", -2, "vector", 214,
                  typeid(vector<vhtm::PackedPFCandidate>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLPackedPFCandidategR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::PackedPFCandidate>) );
      instance.SetNew(&new_vectorlEvhtmcLcLPackedPFCandidategR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLPackedPFCandidategR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLPackedPFCandidategR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLPackedPFCandidategR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLPackedPFCandidategR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::PackedPFCandidate> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::PackedPFCandidate>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLPackedPFCandidategR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::PackedPFCandidate>*)0x0)->GetClass();
      vectorlEvhtmcLcLPackedPFCandidategR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLPackedPFCandidategR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLPackedPFCandidategR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::PackedPFCandidate> : new vector<vhtm::PackedPFCandidate>;
   }
   static void *newArray_vectorlEvhtmcLcLPackedPFCandidategR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::PackedPFCandidate>[nElements] : new vector<vhtm::PackedPFCandidate>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLPackedPFCandidategR(void *p) {
      delete ((vector<vhtm::PackedPFCandidate>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLPackedPFCandidategR(void *p) {
      delete [] ((vector<vhtm::PackedPFCandidate>*)p);
   }
   static void destruct_vectorlEvhtmcLcLPackedPFCandidategR(void *p) {
      typedef vector<vhtm::PackedPFCandidate> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::PackedPFCandidate>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLMuongR_Dictionary();
   static void vectorlEvhtmcLcLMuongR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLMuongR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLMuongR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLMuongR(void *p);
   static void deleteArray_vectorlEvhtmcLcLMuongR(void *p);
   static void destruct_vectorlEvhtmcLcLMuongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Muon>*)
   {
      vector<vhtm::Muon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Muon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Muon>", -2, "vector", 214,
                  typeid(vector<vhtm::Muon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLMuongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Muon>) );
      instance.SetNew(&new_vectorlEvhtmcLcLMuongR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLMuongR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLMuongR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLMuongR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLMuongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Muon> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Muon>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLMuongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Muon>*)0x0)->GetClass();
      vectorlEvhtmcLcLMuongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLMuongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLMuongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Muon> : new vector<vhtm::Muon>;
   }
   static void *newArray_vectorlEvhtmcLcLMuongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Muon>[nElements] : new vector<vhtm::Muon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLMuongR(void *p) {
      delete ((vector<vhtm::Muon>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLMuongR(void *p) {
      delete [] ((vector<vhtm::Muon>*)p);
   }
   static void destruct_vectorlEvhtmcLcLMuongR(void *p) {
      typedef vector<vhtm::Muon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Muon>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLMETgR_Dictionary();
   static void vectorlEvhtmcLcLMETgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLMETgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLMETgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLMETgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLMETgR(void *p);
   static void destruct_vectorlEvhtmcLcLMETgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::MET>*)
   {
      vector<vhtm::MET> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::MET>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::MET>", -2, "vector", 214,
                  typeid(vector<vhtm::MET>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLMETgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::MET>) );
      instance.SetNew(&new_vectorlEvhtmcLcLMETgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLMETgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLMETgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLMETgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLMETgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::MET> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::MET>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLMETgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::MET>*)0x0)->GetClass();
      vectorlEvhtmcLcLMETgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLMETgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLMETgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::MET> : new vector<vhtm::MET>;
   }
   static void *newArray_vectorlEvhtmcLcLMETgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::MET>[nElements] : new vector<vhtm::MET>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLMETgR(void *p) {
      delete ((vector<vhtm::MET>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLMETgR(void *p) {
      delete [] ((vector<vhtm::MET>*)p);
   }
   static void destruct_vectorlEvhtmcLcLMETgR(void *p) {
      typedef vector<vhtm::MET> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::MET>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLJetgR_Dictionary();
   static void vectorlEvhtmcLcLJetgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLJetgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLJetgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLJetgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLJetgR(void *p);
   static void destruct_vectorlEvhtmcLcLJetgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Jet>*)
   {
      vector<vhtm::Jet> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Jet>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Jet>", -2, "vector", 214,
                  typeid(vector<vhtm::Jet>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLJetgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Jet>) );
      instance.SetNew(&new_vectorlEvhtmcLcLJetgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLJetgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLJetgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLJetgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLJetgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Jet> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Jet>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLJetgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Jet>*)0x0)->GetClass();
      vectorlEvhtmcLcLJetgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLJetgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLJetgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Jet> : new vector<vhtm::Jet>;
   }
   static void *newArray_vectorlEvhtmcLcLJetgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Jet>[nElements] : new vector<vhtm::Jet>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLJetgR(void *p) {
      delete ((vector<vhtm::Jet>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLJetgR(void *p) {
      delete [] ((vector<vhtm::Jet>*)p);
   }
   static void destruct_vectorlEvhtmcLcLJetgR(void *p) {
      typedef vector<vhtm::Jet> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Jet>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLGenParticlegR_Dictionary();
   static void vectorlEvhtmcLcLGenParticlegR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLGenParticlegR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLGenParticlegR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLGenParticlegR(void *p);
   static void deleteArray_vectorlEvhtmcLcLGenParticlegR(void *p);
   static void destruct_vectorlEvhtmcLcLGenParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::GenParticle>*)
   {
      vector<vhtm::GenParticle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::GenParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::GenParticle>", -2, "vector", 214,
                  typeid(vector<vhtm::GenParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLGenParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::GenParticle>) );
      instance.SetNew(&new_vectorlEvhtmcLcLGenParticlegR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLGenParticlegR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLGenParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLGenParticlegR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLGenParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::GenParticle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::GenParticle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLGenParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::GenParticle>*)0x0)->GetClass();
      vectorlEvhtmcLcLGenParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLGenParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLGenParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenParticle> : new vector<vhtm::GenParticle>;
   }
   static void *newArray_vectorlEvhtmcLcLGenParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenParticle>[nElements] : new vector<vhtm::GenParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLGenParticlegR(void *p) {
      delete ((vector<vhtm::GenParticle>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLGenParticlegR(void *p) {
      delete [] ((vector<vhtm::GenParticle>*)p);
   }
   static void destruct_vectorlEvhtmcLcLGenParticlegR(void *p) {
      typedef vector<vhtm::GenParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::GenParticle>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLGenMETgR_Dictionary();
   static void vectorlEvhtmcLcLGenMETgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLGenMETgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLGenMETgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLGenMETgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLGenMETgR(void *p);
   static void destruct_vectorlEvhtmcLcLGenMETgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::GenMET>*)
   {
      vector<vhtm::GenMET> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::GenMET>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::GenMET>", -2, "vector", 214,
                  typeid(vector<vhtm::GenMET>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLGenMETgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::GenMET>) );
      instance.SetNew(&new_vectorlEvhtmcLcLGenMETgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLGenMETgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLGenMETgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLGenMETgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLGenMETgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::GenMET> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::GenMET>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLGenMETgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::GenMET>*)0x0)->GetClass();
      vectorlEvhtmcLcLGenMETgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLGenMETgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLGenMETgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenMET> : new vector<vhtm::GenMET>;
   }
   static void *newArray_vectorlEvhtmcLcLGenMETgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenMET>[nElements] : new vector<vhtm::GenMET>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLGenMETgR(void *p) {
      delete ((vector<vhtm::GenMET>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLGenMETgR(void *p) {
      delete [] ((vector<vhtm::GenMET>*)p);
   }
   static void destruct_vectorlEvhtmcLcLGenMETgR(void *p) {
      typedef vector<vhtm::GenMET> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::GenMET>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLGenJetgR_Dictionary();
   static void vectorlEvhtmcLcLGenJetgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLGenJetgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLGenJetgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLGenJetgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLGenJetgR(void *p);
   static void destruct_vectorlEvhtmcLcLGenJetgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::GenJet>*)
   {
      vector<vhtm::GenJet> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::GenJet>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::GenJet>", -2, "vector", 214,
                  typeid(vector<vhtm::GenJet>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLGenJetgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::GenJet>) );
      instance.SetNew(&new_vectorlEvhtmcLcLGenJetgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLGenJetgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLGenJetgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLGenJetgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLGenJetgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::GenJet> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::GenJet>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLGenJetgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::GenJet>*)0x0)->GetClass();
      vectorlEvhtmcLcLGenJetgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLGenJetgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLGenJetgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenJet> : new vector<vhtm::GenJet>;
   }
   static void *newArray_vectorlEvhtmcLcLGenJetgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenJet>[nElements] : new vector<vhtm::GenJet>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLGenJetgR(void *p) {
      delete ((vector<vhtm::GenJet>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLGenJetgR(void *p) {
      delete [] ((vector<vhtm::GenJet>*)p);
   }
   static void destruct_vectorlEvhtmcLcLGenJetgR(void *p) {
      typedef vector<vhtm::GenJet> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::GenJet>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLGenEventgR_Dictionary();
   static void vectorlEvhtmcLcLGenEventgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLGenEventgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLGenEventgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLGenEventgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLGenEventgR(void *p);
   static void destruct_vectorlEvhtmcLcLGenEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::GenEvent>*)
   {
      vector<vhtm::GenEvent> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::GenEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::GenEvent>", -2, "vector", 214,
                  typeid(vector<vhtm::GenEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLGenEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::GenEvent>) );
      instance.SetNew(&new_vectorlEvhtmcLcLGenEventgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLGenEventgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLGenEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLGenEventgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLGenEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::GenEvent> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::GenEvent>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLGenEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::GenEvent>*)0x0)->GetClass();
      vectorlEvhtmcLcLGenEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLGenEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLGenEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenEvent> : new vector<vhtm::GenEvent>;
   }
   static void *newArray_vectorlEvhtmcLcLGenEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::GenEvent>[nElements] : new vector<vhtm::GenEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLGenEventgR(void *p) {
      delete ((vector<vhtm::GenEvent>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLGenEventgR(void *p) {
      delete [] ((vector<vhtm::GenEvent>*)p);
   }
   static void destruct_vectorlEvhtmcLcLGenEventgR(void *p) {
      typedef vector<vhtm::GenEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::GenEvent>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLEventgR_Dictionary();
   static void vectorlEvhtmcLcLEventgR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLEventgR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLEventgR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLEventgR(void *p);
   static void deleteArray_vectorlEvhtmcLcLEventgR(void *p);
   static void destruct_vectorlEvhtmcLcLEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Event>*)
   {
      vector<vhtm::Event> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Event>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Event>", -2, "vector", 214,
                  typeid(vector<vhtm::Event>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Event>) );
      instance.SetNew(&new_vectorlEvhtmcLcLEventgR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLEventgR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLEventgR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Event> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Event>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Event>*)0x0)->GetClass();
      vectorlEvhtmcLcLEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Event> : new vector<vhtm::Event>;
   }
   static void *newArray_vectorlEvhtmcLcLEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Event>[nElements] : new vector<vhtm::Event>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLEventgR(void *p) {
      delete ((vector<vhtm::Event>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLEventgR(void *p) {
      delete [] ((vector<vhtm::Event>*)p);
   }
   static void destruct_vectorlEvhtmcLcLEventgR(void *p) {
      typedef vector<vhtm::Event> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Event>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLElectrongR_Dictionary();
   static void vectorlEvhtmcLcLElectrongR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLElectrongR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLElectrongR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLElectrongR(void *p);
   static void deleteArray_vectorlEvhtmcLcLElectrongR(void *p);
   static void destruct_vectorlEvhtmcLcLElectrongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Electron>*)
   {
      vector<vhtm::Electron> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Electron>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Electron>", -2, "vector", 214,
                  typeid(vector<vhtm::Electron>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLElectrongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Electron>) );
      instance.SetNew(&new_vectorlEvhtmcLcLElectrongR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLElectrongR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLElectrongR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLElectrongR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLElectrongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Electron> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Electron>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLElectrongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Electron>*)0x0)->GetClass();
      vectorlEvhtmcLcLElectrongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLElectrongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLElectrongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Electron> : new vector<vhtm::Electron>;
   }
   static void *newArray_vectorlEvhtmcLcLElectrongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Electron>[nElements] : new vector<vhtm::Electron>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLElectrongR(void *p) {
      delete ((vector<vhtm::Electron>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLElectrongR(void *p) {
      delete [] ((vector<vhtm::Electron>*)p);
   }
   static void destruct_vectorlEvhtmcLcLElectrongR(void *p) {
      typedef vector<vhtm::Electron> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Electron>

namespace ROOT {
   static TClass *vectorlEvhtmcLcLCandidategR_Dictionary();
   static void vectorlEvhtmcLcLCandidategR_TClassManip(TClass*);
   static void *new_vectorlEvhtmcLcLCandidategR(void *p = 0);
   static void *newArray_vectorlEvhtmcLcLCandidategR(Long_t size, void *p);
   static void delete_vectorlEvhtmcLcLCandidategR(void *p);
   static void deleteArray_vectorlEvhtmcLcLCandidategR(void *p);
   static void destruct_vectorlEvhtmcLcLCandidategR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vhtm::Candidate>*)
   {
      vector<vhtm::Candidate> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vhtm::Candidate>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vhtm::Candidate>", -2, "vector", 214,
                  typeid(vector<vhtm::Candidate>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvhtmcLcLCandidategR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vhtm::Candidate>) );
      instance.SetNew(&new_vectorlEvhtmcLcLCandidategR);
      instance.SetNewArray(&newArray_vectorlEvhtmcLcLCandidategR);
      instance.SetDelete(&delete_vectorlEvhtmcLcLCandidategR);
      instance.SetDeleteArray(&deleteArray_vectorlEvhtmcLcLCandidategR);
      instance.SetDestructor(&destruct_vectorlEvhtmcLcLCandidategR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vhtm::Candidate> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vhtm::Candidate>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvhtmcLcLCandidategR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vhtm::Candidate>*)0x0)->GetClass();
      vectorlEvhtmcLcLCandidategR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvhtmcLcLCandidategR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvhtmcLcLCandidategR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Candidate> : new vector<vhtm::Candidate>;
   }
   static void *newArray_vectorlEvhtmcLcLCandidategR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vhtm::Candidate>[nElements] : new vector<vhtm::Candidate>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvhtmcLcLCandidategR(void *p) {
      delete ((vector<vhtm::Candidate>*)p);
   }
   static void deleteArray_vectorlEvhtmcLcLCandidategR(void *p) {
      delete [] ((vector<vhtm::Candidate>*)p);
   }
   static void destruct_vectorlEvhtmcLcLCandidategR(void *p) {
      typedef vector<vhtm::Candidate> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vhtm::Candidate>

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 214,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 214,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 214,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETLorentzVectorgR_Dictionary();
   static void vectorlETLorentzVectorgR_TClassManip(TClass*);
   static void *new_vectorlETLorentzVectorgR(void *p = 0);
   static void *newArray_vectorlETLorentzVectorgR(Long_t size, void *p);
   static void delete_vectorlETLorentzVectorgR(void *p);
   static void deleteArray_vectorlETLorentzVectorgR(void *p);
   static void destruct_vectorlETLorentzVectorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TLorentzVector>*)
   {
      vector<TLorentzVector> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TLorentzVector>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TLorentzVector>", -2, "vector", 214,
                  typeid(vector<TLorentzVector>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETLorentzVectorgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TLorentzVector>) );
      instance.SetNew(&new_vectorlETLorentzVectorgR);
      instance.SetNewArray(&newArray_vectorlETLorentzVectorgR);
      instance.SetDelete(&delete_vectorlETLorentzVectorgR);
      instance.SetDeleteArray(&deleteArray_vectorlETLorentzVectorgR);
      instance.SetDestructor(&destruct_vectorlETLorentzVectorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TLorentzVector> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TLorentzVector>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETLorentzVectorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TLorentzVector>*)0x0)->GetClass();
      vectorlETLorentzVectorgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETLorentzVectorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETLorentzVectorgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLorentzVector> : new vector<TLorentzVector>;
   }
   static void *newArray_vectorlETLorentzVectorgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLorentzVector>[nElements] : new vector<TLorentzVector>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETLorentzVectorgR(void *p) {
      delete ((vector<TLorentzVector>*)p);
   }
   static void deleteArray_vectorlETLorentzVectorgR(void *p) {
      delete [] ((vector<TLorentzVector>*)p);
   }
   static void destruct_vectorlETLorentzVectorgR(void *p) {
      typedef vector<TLorentzVector> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TLorentzVector>

namespace ROOT {
   static TClass *maplEstringcOvectorlEdoublegRsPgR_Dictionary();
   static void maplEstringcOvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlEdoublegRsPgR(void *p = 0);
   static void *newArray_maplEstringcOvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlEdoublegRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlEdoublegRsPgR(void *p);
   static void destruct_maplEstringcOvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<double> >*)
   {
      map<string,vector<double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<double> >", -2, "map", 96,
                  typeid(map<string,vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlEdoublegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,vector<double> >) );
      instance.SetNew(&new_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<double> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,vector<double> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,vector<double> >*)0x0)->GetClass();
      maplEstringcOvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<double> > : new map<string,vector<double> >;
   }
   static void *newArray_maplEstringcOvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<double> >[nElements] : new map<string,vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      delete ((map<string,vector<double> >*)p);
   }
   static void deleteArray_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      delete [] ((map<string,vector<double> >*)p);
   }
   static void destruct_maplEstringcOvectorlEdoublegRsPgR(void *p) {
      typedef map<string,vector<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<double> >

namespace ROOT {
   static TClass *maplEstringcOunsignedsPintgR_Dictionary();
   static void maplEstringcOunsignedsPintgR_TClassManip(TClass*);
   static void *new_maplEstringcOunsignedsPintgR(void *p = 0);
   static void *newArray_maplEstringcOunsignedsPintgR(Long_t size, void *p);
   static void delete_maplEstringcOunsignedsPintgR(void *p);
   static void deleteArray_maplEstringcOunsignedsPintgR(void *p);
   static void destruct_maplEstringcOunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,unsigned int>*)
   {
      map<string,unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,unsigned int>", -2, "map", 96,
                  typeid(map<string,unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOunsignedsPintgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,unsigned int>) );
      instance.SetNew(&new_maplEstringcOunsignedsPintgR);
      instance.SetNewArray(&newArray_maplEstringcOunsignedsPintgR);
      instance.SetDelete(&delete_maplEstringcOunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOunsignedsPintgR);
      instance.SetDestructor(&destruct_maplEstringcOunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,unsigned int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,unsigned int>*)0x0)->GetClass();
      maplEstringcOunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,unsigned int> : new map<string,unsigned int>;
   }
   static void *newArray_maplEstringcOunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,unsigned int>[nElements] : new map<string,unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOunsignedsPintgR(void *p) {
      delete ((map<string,unsigned int>*)p);
   }
   static void deleteArray_maplEstringcOunsignedsPintgR(void *p) {
      delete [] ((map<string,unsigned int>*)p);
   }
   static void destruct_maplEstringcOunsignedsPintgR(void *p) {
      typedef map<string,unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,unsigned int>

namespace ROOT {
   static TClass *maplEstringcOfloatgR_Dictionary();
   static void maplEstringcOfloatgR_TClassManip(TClass*);
   static void *new_maplEstringcOfloatgR(void *p = 0);
   static void *newArray_maplEstringcOfloatgR(Long_t size, void *p);
   static void delete_maplEstringcOfloatgR(void *p);
   static void deleteArray_maplEstringcOfloatgR(void *p);
   static void destruct_maplEstringcOfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,float>*)
   {
      map<string,float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,float>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,float>", -2, "map", 96,
                  typeid(map<string,float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,float>) );
      instance.SetNew(&new_maplEstringcOfloatgR);
      instance.SetNewArray(&newArray_maplEstringcOfloatgR);
      instance.SetDelete(&delete_maplEstringcOfloatgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOfloatgR);
      instance.SetDestructor(&destruct_maplEstringcOfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,float>*)0x0)->GetClass();
      maplEstringcOfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,float> : new map<string,float>;
   }
   static void *newArray_maplEstringcOfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,float>[nElements] : new map<string,float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOfloatgR(void *p) {
      delete ((map<string,float>*)p);
   }
   static void deleteArray_maplEstringcOfloatgR(void *p) {
      delete [] ((map<string,float>*)p);
   }
   static void destruct_maplEstringcOfloatgR(void *p) {
      typedef map<string,float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,float>

namespace ROOT {
   static TClass *maplEdoublecOvectorlEdoublegRsPgR_Dictionary();
   static void maplEdoublecOvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_maplEdoublecOvectorlEdoublegRsPgR(void *p = 0);
   static void *newArray_maplEdoublecOvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_maplEdoublecOvectorlEdoublegRsPgR(void *p);
   static void deleteArray_maplEdoublecOvectorlEdoublegRsPgR(void *p);
   static void destruct_maplEdoublecOvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<double,vector<double> >*)
   {
      map<double,vector<double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<double,vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<double,vector<double> >", -2, "map", 96,
                  typeid(map<double,vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEdoublecOvectorlEdoublegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<double,vector<double> >) );
      instance.SetNew(&new_maplEdoublecOvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_maplEdoublecOvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_maplEdoublecOvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEdoublecOvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_maplEdoublecOvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<double,vector<double> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<double,vector<double> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEdoublecOvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<double,vector<double> >*)0x0)->GetClass();
      maplEdoublecOvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEdoublecOvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEdoublecOvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<double,vector<double> > : new map<double,vector<double> >;
   }
   static void *newArray_maplEdoublecOvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<double,vector<double> >[nElements] : new map<double,vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEdoublecOvectorlEdoublegRsPgR(void *p) {
      delete ((map<double,vector<double> >*)p);
   }
   static void deleteArray_maplEdoublecOvectorlEdoublegRsPgR(void *p) {
      delete [] ((map<double,vector<double> >*)p);
   }
   static void destruct_maplEdoublecOvectorlEdoublegRsPgR(void *p) {
      typedef map<double,vector<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<double,vector<double> >

namespace {
  void TriggerDictionaryInitialization_Dict_Impl() {
    static const char* headers[] = {
"interface/PhysicsObjects.h",
0
    };
    static const char* includePaths[] = {
"./interface",
"./",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08-elfike/include",
"/afs/cern.ch/work/g/gsaha/public/MultileptonJetMetAnalysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$bits/char_traits.h")))  __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace std{template <typename > class __attribute__((annotate("$clingAutoload$bits/memoryfwd.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Muon;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Candidate;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  TriggerObject;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Photon;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Vertex;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Jet;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Tau;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  MET;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  GenMET;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  GenJet;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  GenParticle;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Electron;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  GenEvent;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  Event;}
namespace vhtm{class __attribute__((annotate("$clingAutoload$interface/PhysicsObjects.h")))  PackedPFCandidate;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "interface/PhysicsObjects.h"
#include "interface/PhysicsObjects.h"

#ifdef __CINT__

#include <vector>
#include <map>
#include <string>
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ class vhtm::Candidate+;
#pragma link C++ class std::vector<vhtm::Candidate>+;
#pragma link C++ class vhtm::PackedPFCandidate+;
#pragma link C++ class std::vector<vhtm::PackedPFCandidate>+;
#pragma link C++ class vhtm::Event+;
#pragma link C++ class std::vector<vhtm::Event>+;
#pragma link C++ class vhtm::GenEvent+;
#pragma link C++ class std::vector<vhtm::GenEvent>+;
#pragma link C++ class vhtm::Electron+;
#pragma link C++ class std::vector<vhtm::Electron>+;
#pragma link C++ class vhtm::GenParticle+;
#pragma link C++ class std::vector<vhtm::GenParticle>+;
#pragma link C++ class vhtm::GenJet+;
#pragma link C++ class std::vector<vhtm::GenJet>+;
#pragma link C++ class vhtm::GenMET+;
#pragma link C++ class std::vector<vhtm::GenMET>+;
#pragma link C++ class vhtm::MET+;
#pragma link C++ class std::vector<vhtm::MET>+;
#pragma link C++ class vhtm::Tau+;
#pragma link C++ class std::vector<vhtm::Tau>+;
#pragma link C++ class vhtm::Muon+;
#pragma link C++ class std::vector<vhtm::Muon>+;
#pragma link C++ class std::pair<vhtm::Muon,vhtm::Muon>+;
#pragma link C++ class vhtm::Jet+;
#pragma link C++ class std::vector<vhtm::Jet>+;
#pragma link C++ class vhtm::Vertex+;
#pragma link C++ class std::vector<vhtm::Vertex>+;
#pragma link C++ class vhtm::Photon+;
#pragma link C++ class std::vector<vhtm::Photon>+;
#pragma link C++ class vhtm::TriggerObject+;
#pragma link C++ class std::vector<vhtm::TriggerObject>+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class std::map<std::string,unsigned int>+;
#pragma link C++ class std::pair<std::string,unsigned int>+;
#pragma link C++ class std::map<double,std::vector<double> >+;
#pragma link C++ class std::map<std::string,std::vector<double> >+;
#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"vhtm::Candidate", payloadCode, "@",
"vhtm::Electron", payloadCode, "@",
"vhtm::Event", payloadCode, "@",
"vhtm::GenEvent", payloadCode, "@",
"vhtm::GenJet", payloadCode, "@",
"vhtm::GenMET", payloadCode, "@",
"vhtm::GenParticle", payloadCode, "@",
"vhtm::Jet", payloadCode, "@",
"vhtm::MET", payloadCode, "@",
"vhtm::Muon", payloadCode, "@",
"vhtm::PackedPFCandidate", payloadCode, "@",
"vhtm::Photon", payloadCode, "@",
"vhtm::Tau", payloadCode, "@",
"vhtm::TriggerObject", payloadCode, "@",
"vhtm::Vertex", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Dict() {
  TriggerDictionaryInitialization_Dict_Impl();
}
