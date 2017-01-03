#ifndef _TTHBBLEPTONICREADER_H_
#define _TTHBBLEPTONICREADER_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <map>
#include <vector>

#include "TString.h"
#include "TTree.h"
#include "TDataType.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TLeafI.h"

#include "Report.h"
#include "tth.h"

/*
#include "TTHbbLeptonic/MVAVariables.h"
#include "TTHbbLeptonic/PairedSystem.h"
#include "TopEvent/Event.h"
*/

using namespace std;
using namespace tth;

static Bool_t youHaveBeenWarned_IsolationVars = kFALSE;

class
TTHbbLeptonicReader : public TreeReader
{

private:

  TString selectionTag;
  TString baseCategory;

  Double_t weight;
  
  Double_t minPtEl;
  Double_t minPtMu;
  Double_t minPtJet;
  Double_t minPtLeadingEl;
  Double_t minPtLeadingMu;
  UInt_t minNjets;
  UInt_t minNbtags;
  Double_t maxEtaEl;
  Double_t maxEtaMu;
  Double_t maxYJet;
  
  TString qualityEl;
  TString qualityMu;
  TString qualityLeadingEl;
  TString qualityLeadingMu;

  TString isoEl;
  TString isoMu;
  TString isoLeadingEl;
  TString isoLeadingMu;

  Bool_t useTriggerMatching;
  Bool_t useQualityIsoCuts;

  Int_t passLjet; // 1 => Require events to also pass L+jets SR; -1 => Require to fail L+jets; other => No extra selection
  Int_t truthMatch; // 1 => Require truth-matched leptons; -1 => Require fake leptons; other => No extra selection
  Bool_t allLeptonsMatchTruth;

  //top::Event m_event;
  
  std::vector<TLorentzVector> good_jets;
  std::vector<TLorentzVector> good_jets_btagged;
  std::vector<TLorentzVector> good_el;
  std::vector<TLorentzVector> good_mu;
  std::vector<std::size_t> good_el_id;
  std::vector<std::size_t> good_el_id_forTriggerMatching;
  std::vector<std::size_t> good_mu_id;
  std::vector<std::size_t> good_mu_id_forTriggerMatching;

  std::vector<TLorentzVector> good_el_SL;
  std::vector<TLorentzVector> good_mu_SL;
  std::vector<std::size_t> good_el_id_SL;
  std::vector<std::size_t> good_el_id_forTriggerMatching_SL;
  std::vector<std::size_t> good_mu_id_SL;
  std::vector<std::size_t> good_mu_id_forTriggerMatching_SL;
  
  std::vector<Trig> eTriggers  = { Trig("HLT_e24_lhmedium_L1EM20VH",1,0) , Trig("HLT_e60_lhmedium",1,0) , Trig("HLT_e120_lhloose",1,0) };
  std::vector<Trig> uTriggers  = { Trig("HLT_mu20_iloose_L1MU15",0,1) , Trig("HLT_mu40",0,1) , Trig("HLT_mu60_0eta105_msonly",0,1) };
  std::vector<Trig> eeTriggers = { Trig("HLT_2e12_lhloose_L12EM10VH",2,0) };
  std::vector<Trig> uuTriggers = { Trig("HLT_2mu10",0,2) , Trig("HLT_mu18_mu8noL1",0,1) };
  std::vector<Trig> euTriggers = { Trig("HLT_e17_loose_mu14",1,1) , Trig("HLT_e24_medium_L1EM20VHI_mu8noL1",1,0) , Trig("HLT_e7_medium_mu24",1,1) };

  // These bools get updated for each event depending on whether the above triggers fire. If useTriggerMatching
  // is set to '0' then trigger matching is not required.
  Bool_t eTriggerPass;
  Bool_t uTriggerPass;
  Bool_t eeTriggerPass;
  Bool_t uuTriggerPass;
  Bool_t euTriggerPass;
  Bool_t eTriggerPassSL;
  Bool_t uTriggerPassSL;
  
  Bool_t triggerFired( Trig *t ) { return (map_char[t->name] == 1); }

  Bool_t electronTriggerMatches( Trig *t , std::vector<std::size_t> *vec_elids ) {
    if( !useTriggerMatching ) return kTRUE;
    if( t->nElectronsL1 == 0 ) return kTRUE;
    Int_t nMatches = 0;
    for( std::size_t i = 0 ; i < vec_elids->size() ; ++i ) {
      //cout << "trig = " << t << " = " << t->elMatchName << endl;
      //cout << "size = " << map_vector_char[t->elMatchName]->size() << endl;
      //cout << "ID   = " << vec_elids->at(i) << endl;
      if( map_vector_char[t->elMatchName]->at(vec_elids->at(i)) == 1 ) nMatches++;
    }
    return nMatches >= t->nElectronsL1;
  }

  Bool_t muonTriggerMatches( Trig *t , std::vector<std::size_t> *vec_muids ) {
    if( !useTriggerMatching ) return kTRUE;
    if( t->nMuonsL1 == 0 ) return kTRUE;
    Int_t nMatches = 0;
    for( std::size_t i = 0 ; i < vec_muids->size() ; ++i ) {
      if( map_vector_char[t->muMatchName]->at(vec_muids->at(i)) == 1 ) nMatches++;
    }
    return nMatches >= t->nMuonsL1;
  }

  Bool_t checkTrigger( Trig *t , std::vector<std::size_t> *vec_elids , std::vector<std::size_t> *vec_muids ) {
    if( ! triggerFired(t) ) return kFALSE;
    if( ! electronTriggerMatches(t,vec_elids) ) return kFALSE;
    if( ! muonTriggerMatches(t,vec_muids) ) return kFALSE;
    return kTRUE;
  }

  Bool_t checkTriggers( std::vector<std::size_t> *vec_elids , std::vector<std::size_t> *vec_muids , std::vector<Trig> *vec_trigs ) {
    for( std::size_t itrig = 0 ; itrig < vec_trigs->size() ; ++itrig ) {
      if( checkTrigger( &(vec_trigs->at(itrig)) , vec_elids , vec_muids ) ) return kTRUE;
    }
    return kFALSE;
  }

  // Checks to see if leading leptons pass quality/isolation cuts
  Bool_t goodLeadingEl() {
    if( good_el.size() == 0 ) return kFALSE;
    if( useQualityIsoCuts ) {
      if( map_vector_int[qualityLeadingEl]->at(good_el_id[0]) != 1 ) return kFALSE;
      if( map_vector_char[isoLeadingEl]->at(good_el_id[0]) != 1 ) return kFALSE;
    }
    if( good_el[0].Pt() < minPtLeadingEl ) return kFALSE;
    return kTRUE;
  }
  Bool_t goodLeadingMu() {
    if( good_mu.size() == 0 ) return kFALSE;
    if( useQualityIsoCuts ) {
      if( map_vector_char[qualityLeadingMu]->at(good_mu_id[0]) != 1 ) return kFALSE;
      if( map_vector_char[isoLeadingMu]->at(good_mu_id[0]) != 1 ) return kFALSE;
    }
    if( good_mu[0].Pt() < minPtLeadingMu ) return kFALSE;
    return kTRUE;
  }

  // Generic selections
  Bool_t sel_ee_withz() { return goodLeadingEl() && good_el.size()==2 && good_mu.size()==0 && map_float["good_Mee"] > 15000; }
  Bool_t sel_uu_withz() { return goodLeadingMu() && good_el.size()==0 && good_mu.size()==2 && map_float["good_Muu"] > 15000; }
  Bool_t sel_ee() { return sel_ee_withz() && (map_float["good_Mee"] < 83000 || map_float["good_Mee"] > 99000); }
  Bool_t sel_uu() { return sel_uu_withz() && (map_float["good_Muu"] < 83000 || map_float["good_Muu"] > 99000); }
  Bool_t sel_eu() { return (goodLeadingEl() || goodLeadingMu()) && good_el.size()==1 && good_mu.size()==1 && map_float["good_Ht"] > 130000; }

  Bool_t goodLeadingEl_SL() {
    if( good_el_SL.size()==0 ) return kFALSE;
    if( good_el_SL[0].Pt() < 30000 ) return kFALSE;
    return kTRUE;
  }

  Bool_t goodLeadingMu_SL() {
    if( good_mu_SL.size()==0 ) return kFALSE;
    if( good_mu_SL[0].Pt() < 30000 ) return kFALSE;
    return kTRUE;
  }
  
  Bool_t sel_SL_SR() { return (good_jets.size()==5 && good_jets_btagged.size()>=4) || (good_jets.size()>=6 && good_jets_btagged.size()>=3); }
  Bool_t sel_e() { return goodLeadingEl_SL() && good_el_SL.size()==1 && good_mu_SL.size()==0 && eTriggerPassSL && sel_SL_SR(); }
  Bool_t sel_u() { return goodLeadingMu_SL() && good_el_SL.size()==0 && good_mu_SL.size()==1 && uTriggerPassSL && sel_SL_SR(); }
  Bool_t checkSL() {
    if( passLjet==0 ) return kTRUE;
    if( passLjet==1 ) return sel_e() || sel_u();
    if( passLjet==-1 ) return !sel_e() && !sel_u();
    assert( false );
  }
  
public:

  TTHbbLeptonicReader()
    : TreeReader()
    , selectionTag( "none" )
    , baseCategory( "none" )
    , minPtEl( 0 )
    , minPtMu( 0 )
    , minPtJet( 0 )
    , minPtLeadingEl( 30000 )
    , minPtLeadingMu( 30000 )
    , minNjets( 0 )
    , minNbtags( 0 )
    , maxEtaEl( 2.47 )
    , maxEtaMu( 2.5 )
    , maxYJet( 2.5 )
    , useTriggerMatching( kTRUE )
    , useQualityIsoCuts( kTRUE )
    , qualityEl( "el_LHMedium" )
    , qualityMu( "mu_Medium" )
    , qualityLeadingEl( "el_LHMedium" )
    , qualityLeadingMu( "mu_Medium" )
    , isoEl( "el_isoLoose" )
    , isoMu( "mu_isoLoose" )
    , isoLeadingEl( "el_isoLoose" )
    , isoLeadingMu( "mu_isoLoose" )
  {}
  ~TTHbbLeptonicReader() {}

  void setSelectionTag( TString t ) { selectionTag = t; }
  void setBaseCategory( TString t ) { baseCategory = t; }

  void setMinPtEl( Double_t v ) { minPtEl = v; }
  void setMinPtMu( Double_t v ) { minPtMu = v; }
  void setMinPtLeadingLep( Double_t v ) { minPtLeadingEl = v; minPtLeadingMu = v; }
  void setMinPtJet( Double_t v ) { minPtJet = v; }
  void setMinNjets( UInt_t v ) { minNjets = v; }
  void setMinNbtags( UInt_t v ) { minNbtags = v; }
  void setMaxEtaEl( Double_t v ) { maxEtaEl = v; }
  void setMaxEtaMu( Double_t v ) { maxEtaMu = v; }
  void setMaxYJet( Double_t v ) { maxYJet = v; }
  void setUseTriggerMatching( Bool_t v ) { useTriggerMatching = v; }

  void setQualityLep( TString v ) { qualityEl = "el_LH"+v; qualityMu = "mu_"+v; }
  void setQualityLeadingLep( TString v ) { qualityLeadingEl = "el_LH"+v; qualityLeadingMu = "mu_"+v; }
  void setIsoLep( TString v ) { isoEl = "el_iso"+v; isoMu = "mu_iso"+v; }
  void setIsoLeadingLep( TString v ) { isoLeadingEl = "el_iso"+v; isoLeadingMu = "mu_iso"+v; }

  void setPassLjet( Int_t v ) { passLjet = v; }
  void setTruthMatch( Int_t v ) { truthMatch = v; }
  
  TString getSelectionTag( TString t ) { return selectionTag; }


  void setupOutputTree( const TString &tname ) {

    // Initialize output TTree
    TTree *tmp_tree = new TTree( tname , tname );
    tmp_tree->SetDirectory( 0 );
    tmp_tree->Branch( "weight" , &weight );
    
    tmp_tree->Branch( "H1_all" , &map_float["H1_all"] );
    tmp_tree->Branch( "H2_all" , &map_float["H2_all"] );
    tmp_tree->Branch( "H3_all" , &map_float["H3_all"] );
    tmp_tree->Branch( "H4_all" , &map_float["H4_all"] );
    tmp_tree->Branch( "H5_all" , &map_float["H5_all"] );
    tmp_tree->Branch( "H1transverse_all" , &map_float["H1transverse_all"] );
    tmp_tree->Branch( "H2transverse_all" , &map_float["H2transverse_all"] );
    tmp_tree->Branch( "H3transverse_all" , &map_float["H3transverse_all"] );
    tmp_tree->Branch( "H4transverse_all" , &map_float["H4transverse_all"] );
    tmp_tree->Branch( "H5transverse_all" , &map_float["H5transverse_all"] );

    tmp_tree->Branch( "Thrust_all" , &map_float["Thrust_all"] );
    tmp_tree->Branch( "ThrustAxisX_all" , &map_float["ThrustAxisX_all"] );
    tmp_tree->Branch( "ThrustAxisY_all" , &map_float["ThrustAxisY_all"] );
    tmp_tree->Branch( "ThrustAxisZ_all" , &map_float["ThrustAxisZ_all"] );

    for( variableIter iv = mva_variable_names.begin() , fv = mva_variable_names.end() ; iv != fv ; ++iv ) {
      for( pairingIter ip = mva_pairing_names.begin() , fp = mva_pairing_names.end() ; ip != fp ; ++ip ) {
	TString itag = ip->second + "_" + iv->second;
	tmp_tree->Branch( "M"+itag , &map_float["M"+itag] );
	tmp_tree->Branch( "Pt"+itag , &map_float["Pt"+itag] );
	tmp_tree->Branch( "PtSum"+itag , &map_float["PtSum"+itag] );
	tmp_tree->Branch( "dR"+itag , &map_float["dR"+itag] );
	tmp_tree->Branch( "dPhi"+itag , &map_float["dPhi"+itag] );
	tmp_tree->Branch( "dEta"+itag , &map_float["dEta"+itag] );
      }
      tmp_tree->Branch( "Mjjj_"+iv->second , &map_float["Mjjj_"+iv->second] );
      tmp_tree->Branch( "Ptjjj_"+iv->second , &map_float["Ptjjj_"+iv->second] );
    }
    
    tmp_tree->Branch( "DileptonMass" , &map_float["DileptonMass"] );
    tmp_tree->Branch( "DileptonPt" , &map_float["DileptonPt"] );
    tmp_tree->Branch( "DileptonSumPt" , &map_float["DileptonSumPt"] );
    tmp_tree->Branch( "DileptondR" , &map_float["DileptondR"] );
    tmp_tree->Branch( "DileptondPhi" , &map_float["DileptondPhi"] );
    tmp_tree->Branch( "DileptondEta" , &map_float["DileptondEta"] );
    
    for( collectionIter ic = mva_collection_names.begin() , fc = mva_collection_names.end() ; ic != fc ; ++ic ) {
      tmp_tree->Branch( "Aplanarity_"+ic->second , &map_float["Aplanarity_"+ic->second] );
      tmp_tree->Branch( "Aplanority_"+ic->second , &map_float["Aplanority_"+ic->second] );
      tmp_tree->Branch( "Sphericity_"+ic->second , &map_float["Sphericity_"+ic->second] );
      tmp_tree->Branch( "Spherocity_"+ic->second , &map_float["Spherocity_"+ic->second] );
      tmp_tree->Branch( "SphericityT_"+ic->second , &map_float["SphericityT_"+ic->second] );
      tmp_tree->Branch( "Planarity_"+ic->second , &map_float["Planarity_"+ic->second] );
      tmp_tree->Branch( "Variable_C_"+ic->second , &map_float["Variable_C_"+ic->second] );
      tmp_tree->Branch( "Variable_D_"+ic->second , &map_float["Variable_D_"+ic->second] );
      tmp_tree->Branch( "Circularity_"+ic->second , &map_float["Circularity_"+ic->second] );
      tmp_tree->Branch( "PlanarFlow_"+ic->second , &map_float["PlanarFlow_"+ic->second] );
    }
    
    tmp_tree->Branch( "dRlepbb_MindR" , &map_float["dRlepbb_MindR"] );

    outputTrees.push_back( tmp_tree );

  }


  void fillOutputTree( const std::size_t &ifile , const Long64_t &ievent , const Double_t &w ) {
    outputTrees.back()->Fill();
  }

  
  void saveOutputTrees( const TString &fname ) {
    for( std::size_t i = 0 ; i < outputTrees.size() ; i++ ) {
      TString fullfname = TString::Format("output/%s_%s.root",fname.Data(),outputTrees[i]->GetName());
      report::debug( "Saving output tree %s to file %s" , outputTrees[i]->GetName() , fullfname.Data() );
      TFile *fout = new TFile( fullfname , "recreate" );
      outputTrees[i]->SetDirectory( fout );
      outputTrees[i]->Write( "tthAnaOutput" );
      delete outputTrees[i];
      fout->Close();
      delete fout;
    }
  }
  
  
  void setTree( TTree *_tree ) {

    setTree_Generic( _tree );
    
    // Additional observables calcuated in getEntry function
    setBranch_uint( "good_njets" , kFALSE );
    setBranch_uint( "good_nbtags" , kFALSE );
    setBranch_float( "good_Ht" , kFALSE );
    setBranch_float( "good_Mee" , kFALSE );
    setBranch_float( "good_ptee" , kFALSE );
    setBranch_float( "good_Muu" , kFALSE );
    setBranch_float( "good_ptuu" , kFALSE );
    setBranch_float( "good_Meu" , kFALSE );
    setBranch_float( "good_pteu" , kFALSE );
    setBranch_vector_float( "good_jet_pt" , kFALSE );
    setBranch_vector_float( "good_jet_eta" , kFALSE );
    setBranch_vector_float( "good_jet_phi" , kFALSE );
    setBranch_vector_float( "good_jet_y" , kFALSE );
    setBranch_vector_float( "good_el_pt" , kFALSE );
    setBranch_vector_float( "good_el_eta" , kFALSE );
    setBranch_vector_float( "good_el_phi" , kFALSE );
    setBranch_vector_float( "good_mu_pt" , kFALSE );
    setBranch_vector_float( "good_mu_eta" , kFALSE );
    setBranch_vector_float( "good_mu_phi" , kFALSE );
    setBranch_vector_float( "fake_lep_origin" , kFALSE );
    setBranch_vector_float( "fake_lep_originbkg" , kFALSE );
    setBranch_vector_float( "fake_lep_type" , kFALSE );
    setBranch_vector_float( "fake_lep_pdg" , kFALSE );

    // fox wolfram moments
    setBranch_float( "H1_all" , kFALSE );
    setBranch_float( "H2_all" , kFALSE );
    setBranch_float( "H3_all" , kFALSE );
    setBranch_float( "H4_all" , kFALSE );
    setBranch_float( "H5_all" , kFALSE );
    setBranch_float( "H1transverse_all" , kFALSE );
    setBranch_float( "H2transverse_all" , kFALSE );
    setBranch_float( "H3transverse_all" , kFALSE );
    setBranch_float( "H4transverse_all" , kFALSE );
    setBranch_float( "H5transverse_all" , kFALSE );

    // Thrust
    setBranch_float( "Thrust_all" , kFALSE );
    setBranch_float( "ThrustAxisX_all" , kFALSE );
    setBranch_float( "ThrustAxisY_all" , kFALSE );
    setBranch_float( "ThrustAxisZ_all" , kFALSE );

    // Two-object systems
    for( variableIter iv = mva_variable_names.begin() , fv = mva_variable_names.end() ; iv != fv ; ++iv ) {
      for( pairingIter ip = mva_pairing_names.begin() , fp = mva_pairing_names.end() ; ip != fp ; ++ip ) {
	TString itag = ip->second + "_" + iv->second;
	setBranch_float( "M"+itag , kFALSE );
	setBranch_float( "Pt"+itag , kFALSE );
	setBranch_float( "PtSum"+itag , kFALSE );
	setBranch_float( "dR"+itag , kFALSE );
	setBranch_float( "dPhi"+itag , kFALSE );
	setBranch_float( "dEta"+itag , kFALSE );
      }
      setBranch_float( "Mjjj_"+iv->second , kFALSE );
      setBranch_float( "Ptjjj_"+iv->second , kFALSE );
    }

    // DIL
    setBranch_float( "DileptonMass" , kFALSE );
    setBranch_float( "DileptonPt" , kFALSE );
    setBranch_float( "DileptonSumPt" , kFALSE );
    setBranch_float( "DileptondR" , kFALSE );
    setBranch_float( "DileptondPhi" , kFALSE );
    setBranch_float( "DileptondEta" , kFALSE );

    // Collection features
    for( collectionIter ic = mva_collection_names.begin() , fc = mva_collection_names.end() ; ic != fc ; ++ic ) {
      setBranch_float( "Aplanarity_"+ic->second , kFALSE );
      setBranch_float( "Aplanority_"+ic->second , kFALSE );
      setBranch_float( "Sphericity_"+ic->second , kFALSE );
      setBranch_float( "Spherocity_"+ic->second , kFALSE );
      setBranch_float( "SphericityT_"+ic->second , kFALSE );
      setBranch_float( "Planarity_"+ic->second , kFALSE );
      setBranch_float( "Variable_C_"+ic->second , kFALSE );
      setBranch_float( "Variable_D_"+ic->second , kFALSE );
      setBranch_float( "Circularity_"+ic->second , kFALSE );
      setBranch_float( "PlanarFlow_"+ic->second , kFALSE );
    }

    // L+jets
    setBranch_float( "dRlepbb_MindR" , kFALSE );

    // check if the tree has the branches needed for doing lepton quality & isolation cuts
    if( tree->GetEntries() > 0 ) {
      tree->GetEntry(0);
      useQualityIsoCuts = ( map_vector_int[qualityEl] && map_vector_int["el_LHMedium"] &&
			    map_vector_char[isoEl] && map_vector_char["el_isoGradient"] &&
			    map_vector_char[qualityMu] && map_vector_char["mu_Medium"] &&
			    map_vector_char[isoMu] && map_vector_char["mu_isoGradient"] );
      if( !useQualityIsoCuts && !youHaveBeenWarned_IsolationVars ) {
	report::warn( "Could not find lepton isolation/quality branches in the TTree." );
	report::blank( "This isn't critical, but it means no additional cuts will be applied to lepton isolation/quality" );
	report::blank( "in at least some of your samples (this one plus any other one with the missing variables). It" );
	report::blank( "also means that the L+jets signal region selection checks used for some studies will not be done" );
	report::blank( "correctly. Note that you may run into agreement problems if isolation/quality variables are used" );
	report::blank( "in some but not all of your samples, but this shouldn't happen when they're all generated from the" );
	report::blank( "same cuts file." );
	youHaveBeenWarned_IsolationVars = kTRUE;
      }
    }
    
  }


  
  Bool_t getEntry( const Long64_t &ientry ) {

    tree->GetEntry(ientry);
    
    if( baseCategory != TString("none") )
      if( map_int[baseCategory] == 0 )
	return kFALSE;
    
    //
    // Initialize event variables & vectors
    //

    //m_event.clear();
    good_jets.clear();
    good_jets_btagged.clear();
    good_el.clear();
    good_mu.clear();
    good_el_id.clear();
    good_el_id_forTriggerMatching.clear();
    good_mu_id.clear();
    good_mu_id_forTriggerMatching.clear();
    map_vector_float["good_jet_pt"]->clear();
    map_vector_float["good_jet_eta"]->clear();
    map_vector_float["good_jet_phi"]->clear();
    map_vector_float["good_jet_y"]->clear();
    map_vector_float["good_el_pt"]->clear();
    map_vector_float["good_el_eta"]->clear();
    map_vector_float["good_el_phi"]->clear();
    map_vector_float["good_mu_pt"]->clear();
    map_vector_float["good_mu_eta"]->clear();
    map_vector_float["good_mu_phi"]->clear();
    map_vector_float["fake_lep_pdg"]->clear();
    map_vector_float["fake_lep_origin"]->clear();
    map_vector_float["fake_lep_originbkg"]->clear();
    map_vector_float["fake_lep_type"]->clear();

    good_el_SL.clear();
    good_mu_SL.clear();
    good_el_id_SL.clear();
    good_el_id_forTriggerMatching_SL.clear();
    good_mu_id_SL.clear();
    good_mu_id_forTriggerMatching_SL.clear();

    allLeptonsMatchTruth = kTRUE;

    //
    // Find good jets
    //
    
    for( std::size_t i = 0 ; i < map_vector_float["jet_pt"]->size() ; i++ ) {

      Double_t jet_pt = map_vector_float["jet_pt"]->at(i);
      if( jet_pt < minPtJet ) continue;

      TLorentzVector jet;
      jet.SetPtEtaPhiE( jet_pt , map_vector_float["jet_eta"]->at(i) , map_vector_float["jet_phi"]->at(i) , map_vector_float["jet_e"]->at(i) );
      if( TMath::Abs(jet.Rapidity()) > maxYJet ) continue;

      map_vector_float["good_jet_pt"]->push_back( jet.Pt() );
      map_vector_float["good_jet_eta"]->push_back( jet.Eta() );
      map_vector_float["good_jet_phi"]->push_back( jet.Phi() );
      map_vector_float["good_jet_y"]->push_back( jet.Rapidity() );
      
      map_float["good_Ht"] += jet_pt;

      good_jets.push_back( jet );
      if( map_vector_float["jet_mv2c20"]->at(i) > -0.4434 )
	good_jets_btagged.push_back( jet );
      //m_event.m_jets.push_back( new xAOD::Jet( jet , map_vector_float["jet_mv2c20"]->at(i) > -0.4434 ? 10 : -10 ) );
      
    }

    if( good_jets.size() < minNjets ) return kFALSE;
    if( good_jets_btagged.size() < minNbtags ) return kFALSE;

    
    //
    // Find good electrons
    //

    for( std::size_t i = 0 ; i < map_vector_float["el_pt"]->size() ; i++ ) {
	
      // apply quality/isolation cuts, but only if the variables exist in the ntuple
      if( useQualityIsoCuts ) {
	if( map_vector_int[qualityEl]->at(i) != 1 ) continue;
	if( map_vector_char[isoEl]->at(i) != 1 ) continue;
      }

      Double_t el_pt = map_vector_float["el_pt"]->at(i);
      if( el_pt < minPtEl ) continue;

      Double_t el_eta = map_vector_float["el_eta"]->at(i);
      if( TMath::Abs(el_eta) > maxEtaEl ) continue;
      
      TLorentzVector el;
      el.SetPtEtaPhiE( el_pt , el_eta , map_vector_float["el_phi"]->at(i) , map_vector_float["el_e"]->at(i) );

      map_vector_float["good_el_pt"]->push_back( el.Pt() );
      map_vector_float["good_el_eta"]->push_back( el.Eta() );
      map_vector_float["good_el_phi"]->push_back( el.Phi() );
      
      good_el.push_back( el );
      good_el_id.push_back( i );
      //m_event.m_electrons.push_back( new xAOD::Electron(el) );
      if( TMath::Abs(el_eta) <= 1.37 || TMath::Abs(el_eta) >= 1.52 ) {
	// only consider electrons for trigger matching if they are not in the crack region
	good_el_id_forTriggerMatching.push_back( i );
      }
      map_float["good_Ht"] += el_pt;

      // Check if this lepton would satisfy the nominal L+jets selection
      if( useQualityIsoCuts ) {
	if( map_vector_int["el_LHMedium"]->at(i)==1 && map_vector_char["el_isoGradient"]->at(i)==1 ) {
	  good_el_SL.push_back( el );
	  good_el_id_SL.push_back( i );
	  if( TMath::Abs(el_eta)<=1.37 || TMath::Abs(el_eta)>=1.52 ) {
	    good_el_id_forTriggerMatching_SL.push_back( i );
	  }
	}
      }

      // truth-matching to determine if lepton is "real" or "fake"
      if( map_vector_int["el_true_pdg"] ) {
	Bool_t isEl = map_vector_int["el_true_pdg"]->at(i)==11 || map_vector_int["el_true_pdg"]->at(i)==-11;
	Bool_t isIso = map_vector_int["el_true_type"]->at(i)==2 && (map_vector_int["el_true_origin"]->at(i)==10 || map_vector_int["el_true_origin"]->at(i)==12 || map_vector_int["el_true_origin"]->at(i)==13);
	Bool_t isPhotonConv = map_vector_int["el_true_type"]->at(i)==4 && map_vector_int["el_true_origin"]->at(i)==5 && (map_vector_int["el_true_originbkg"]->at(i)==12 || map_vector_int["el_true_originbkg"]->at(i)==13);
	Bool_t isReal = isEl && (isIso || isPhotonConv);
	allLeptonsMatchTruth = allLeptonsMatchTruth && isReal;
	if( !isReal ) {
	  map_vector_float["fake_lep_pdg"]->push_back( float(map_vector_int["el_true_pdg"]->at(i)) );
	  map_vector_float["fake_lep_type"]->push_back( float(map_vector_int["el_true_type"]->at(i)) );
	  map_vector_float["fake_lep_origin"]->push_back( float(map_vector_int["el_true_origin"]->at(i)) );
	  map_vector_float["fake_lep_originbkg"]->push_back( float(map_vector_int["el_true_originbkg"]->at(i)) );
	}
      }
      
    }
    
    //
    // Extra muon cuts
    //
    for( std::size_t i = 0 ; i < map_vector_float["mu_pt"]->size() ; i++ ) {

      if( useQualityIsoCuts ) {
	if( map_vector_char[qualityMu]->at(i) != 1 ) continue;
	if( map_vector_char[isoMu]->at(i) != 1 ) continue;
      }
      
      Double_t mu_pt = map_vector_float["mu_pt"]->at(i);
      if( mu_pt < minPtMu ) continue;

      Double_t mu_eta = map_vector_float["mu_eta"]->at(i);
      if( TMath::Abs(mu_eta) > maxEtaMu ) continue;

      TLorentzVector mu;
      mu.SetPtEtaPhiE( mu_pt , mu_eta , map_vector_float["mu_phi"]->at(i) , map_vector_float["mu_e"]->at(i) );

      map_vector_float["good_mu_pt"]->push_back( mu.Pt() );
      map_vector_float["good_mu_eta"]->push_back( mu.Eta() );
      map_vector_float["good_mu_phi"]->push_back( mu.Phi() );
      
      good_mu.push_back( mu );
      good_mu_id.push_back( i );
      m_event.m_muons.push_back( new xAOD::Muon(mu) );
      if( TMath::Abs(mu_eta) <= 2.4 ) {
	// only consider muons for trigger matching if they have pseudo-rapidity in L1 trigger range
	good_mu_id_forTriggerMatching.push_back( i );
      }
      map_float["good_Ht"] += mu_pt;

      // check if this lepton would satisfy the nominal L+jets selection
      if( useQualityIsoCuts ) {
	if( map_vector_char["mu_Medium"]->at(i)==1 && map_vector_char["mu_isoGradient"]->at(i)==1 ) {
	  good_mu_SL.push_back( mu );
	  good_mu_id_SL.push_back( i );
	  if( TMath::Abs(mu_eta) <= 2.4 ) {
	    good_mu_id_forTriggerMatching_SL.push_back( i );
	  }
	}
      }

      // truth matching to determine if lepton is "real" or "fake" 
      if( map_vector_int["mu_true_pdg"] ) {
	Bool_t isMu = map_vector_int["mu_true_pdg"]->at(i)==13 || map_vector_int["mu_true_pdg"]->at(i)==-13;
	Bool_t isIso = map_vector_int["mu_true_type"]->at(i)==6 && (map_vector_int["mu_true_origin"]->at(i)==10 || map_vector_int["mu_true_origin"]->at(i)==12 || map_vector_int["mu_true_origin"]->at(i)==13);
	Bool_t isReal = isMu && isIso;
	allLeptonsMatchTruth = allLeptonsMatchTruth && isReal;
	if( !isReal ) {
	  map_vector_float["fake_lep_pdg"]->push_back( float(map_vector_int["mu_true_pdg"]->at(i)) );	  
	  map_vector_float["fake_lep_type"]->push_back( float(map_vector_int["mu_true_type"]->at(i)) );
	  map_vector_float["fake_lep_origin"]->push_back( float(map_vector_int["mu_true_origin"]->at(i)) );
	}
      }
      
    }

    m_event.m_met->setP4( map_float["met_met"] , 0.0 , map_float["met_phi"] , 0.0 );
    
    // return here if we didn't find any good leptons
    if( good_mu.size()==0 && good_el.size()==0 ) return kFALSE;

    // return here if...
    //  - leptons are all real and user asks for at least one fake
    //  - leptons are not all real and user asks for no fakes
    if( allLeptonsMatchTruth && truthMatch==-1 ) return kFALSE;
    if( !allLeptonsMatchTruth && truthMatch==1 ) return kFALSE;

    // define some extra variables that we might want to plot in DataMCPlots.cpp
    map_uint["good_njets"] = good_jets.size();
    map_uint["good_nbtags"] = good_jets_btagged.size();
    if( good_el.size() >= 2 ) {
      map_float["good_Mee"] = (good_el[0]+good_el[1]).M();
      map_float["good_ptee"] = (good_el[0]+good_el[1]).Pt();
    }
    if( good_mu.size() >= 2 ) {
      map_float["good_Muu"] = (good_mu[0]+good_mu[1]).M();
      map_float["good_ptuu"] = (good_mu[0]+good_mu[1]).Pt();
    }
    if( good_el.size() >= 1 && good_mu.size() >= 1 ) {
      map_float["good_Meu"] = (good_mu[0]+good_el[0]).M();
      map_float["good_pteu"] = (good_mu[0]+good_el[0]).Pt();
    }
    
    // These functions check that triggers fired and also perform trigger matching
    // if the user has trigger matching turned on (on by default)
    eTriggerPass  = checkTriggers( &good_el_id_forTriggerMatching , &good_mu_id_forTriggerMatching , &eTriggers );
    uTriggerPass  = checkTriggers( &good_el_id_forTriggerMatching , &good_mu_id_forTriggerMatching , &uTriggers );
    eeTriggerPass = checkTriggers( &good_el_id_forTriggerMatching , &good_mu_id_forTriggerMatching , &eeTriggers );
    uuTriggerPass = checkTriggers( &good_el_id_forTriggerMatching , &good_mu_id_forTriggerMatching , &uuTriggers );
    euTriggerPass = checkTriggers( &good_el_id_forTriggerMatching , &good_mu_id_forTriggerMatching , &euTriggers );
    
    // Also check to see if the nominal L+jet triggers have fired
    eTriggerPassSL = checkTriggers( &good_el_id_forTriggerMatching_SL , &good_mu_id_forTriggerMatching_SL , &eTriggers );
    uTriggerPassSL = checkTriggers( &good_el_id_forTriggerMatching_SL , &good_mu_id_forTriggerMatching_SL , &uTriggers );

    // If passLjet!=0, then the checkSL function will return true or false depending on whether the event passes L+jet SR
    // cuts in addition to our own cuts
    if( ! checkSL() ) return kFALSE;

    //
    // We've made it this far in the selection, so lets calculate all the MV training variables...
    //

    /*
    
    MVAVariables *m_mva = new MVAVariables();
    m_mva->initialise( m_event );

    TLorentzVector vleadingLep;
    const xAOD::IParticle *leadingLep = m_mva->getLeadingPtLepton();
    vleadingLep.SetPtEtaPhiE( leadingLep->pt() , leadingLep->eta() , leadingLep->phi() , leadingLep->e() );
    PairedSystem *ps_lepbb_MindR = new PairedSystem( m_mva->getEntry(pairing::bb,variable::MindR) , vleadingLep );

    // Fox Wolfram Moments
    map_float["H1_all"] = m_mva->FirstFoxWolframMoment(collection::all);
    map_float["H2_all"] = m_mva->SecondFoxWolframMoment(collection::all);
    map_float["H3_all"] = m_mva->ThirdFoxWolframMoment(collection::all);
    map_float["H4_all"] = m_mva->FourthFoxWolframMoment(collection::all);
    map_float["H5_all"] = m_mva->FifthFoxWolframMoment(collection::all);
    map_float["H1transverse_all"] = m_mva->FirstFoxWolframTransverseMoment(collection::all);
    map_float["H2transverse_all"] = m_mva->SecondFoxWolframTransverseMoment(collection::all);
    map_float["H3transverse_all"] = m_mva->ThirdFoxWolframTransverseMoment(collection::all);
    map_float["H4transverse_all"] = m_mva->FourthFoxWolframTransverseMoment(collection::all);
    map_float["H5transverse_all"] = m_mva->FifthFoxWolframTransverseMoment(collection::all);

    // Thrust
    map_float["Thrust_all"]	 = m_mva->getThrust(collection::all);
    map_float["ThrustAxisX_all"] = m_mva->getThrustAxis(collection::all).X();
    map_float["ThrustAxisY_all"] = m_mva->getThrustAxis(collection::all).Y();
    map_float["ThrustAxisZ_all"] = m_mva->getThrustAxis(collection::all).Z();

    // Composite objects
    for( variableIter iv = mva_variable_names.begin() , fv = mva_variable_names.end() ; iv != fv ; ++iv ) {
      for( pairingIter ip = mva_pairing_names.begin() , fp = mva_pairing_names.end() ; ip != fp ; ++ip ) {
	TString itag = ip->second + "_" + iv->second;
	map_float["M"+itag] = m_mva->MassofPair( ip->first , iv->first );
	map_float["Pt"+itag] = m_mva->PtofPair( ip->first , iv->first );
	map_float["PtSum"+itag] = m_mva->PtSumofPair( ip->first , iv->first );
	map_float["dR"+itag] = m_mva->deltaRofPair( ip->first , iv->first );
	map_float["dPhi"+itag] = m_mva->deltaPhiofPair( ip->first , iv->first );
	map_float["dEta"+itag] = m_mva->deltaEtaofPair( ip->first , iv->first );
      }
      map_float["Mjjj_"+iv->second] = m_mva->MassofJetTriplet( iv->first );
      map_float["Ptjjj_"+iv->second] = m_mva->PtofJetTriplet( iv->first );
    }

    // DIL
    if( good_mu.size() + good_el.size() > 1 ) {
      map_float["DileptonMass"]	 = m_mva->DileptonMass();
      map_float["DileptonPt"]	 = m_mva->DileptonPt();
      map_float["DileptonSumPt"] = m_mva->DileptonSumPt();
      map_float["DileptondR"]	 = m_mva->DileptondR();
      map_float["DileptondPhi"]	 = m_mva->DileptondPhi();
      map_float["DileptondEta"]	 = m_mva->DileptondEta();
    } else {
      map_float["DileptonMass"]	 = -1;
      map_float["DileptonPt"]	 = -1;
      map_float["DileptonSumPt"] = -1;
      map_float["DileptondR"]	 = -1;
      map_float["DileptondPhi"]	 = -5;
      map_float["DileptondEta"]	 = -10;
    }

    // Collection features
    for( collectionIter ic = mva_collection_names.begin() , fc = mva_collection_names.end() ; ic != fc ; ++ic ) {
      map_float["Aplanarity_"+ic->second] = m_mva->Aplanarity( ic->first );
      map_float["Aplanority_"+ic->second] = m_mva->Aplanority( ic->first );
      map_float["Sphericity_"+ic->second] = m_mva->Sphericity( ic->first );
      map_float["Spherocity_"+ic->second] = m_mva->Spherocity( ic->first );
      map_float["SphericityT_"+ic->second] = m_mva->SphericityT( ic->first );
      map_float["Planarity_"+ic->second] = m_mva->Planarity( ic->first );
      map_float["Variable_C_"+ic->second] = m_mva->Variable_C( ic->first );
      map_float["Variable_D_"+ic->second] = m_mva->Variable_D( ic->first );
      map_float["Circularity_"+ic->second] = m_mva->Circularity( ic->first );
      map_float["PlanarFlow_"+ic->second] = m_mva->PlanarFlow( ic->first );
    }

    // L+jets
    map_float["dRlepbb_MindR"] = ps_lepbb_MindR->DeltaR();
    
    delete ps_lepbb_MindR;
    delete m_mva;

    */
    
    return kTRUE;
    
  }

  
  Bool_t passesSelection() {

    Bool_t triggers_pass	   = ( eTriggerPass || uTriggerPass || eeTriggerPass || uuTriggerPass );
    Bool_t singlelep_triggers_pass = ( eTriggerPass || uTriggerPass );
    
    if( selectionTag == TString("none") ) return kTRUE;
    
    // Nominal selections with m(ll) away from the Z mass
    else if( selectionTag == TString("ee_nominal") ) return( triggers_pass && sel_ee() );
    else if( selectionTag == TString("uu_nominal") ) return( triggers_pass && sel_uu() );
    else if( selectionTag == TString("eu_nominal") ) return( triggers_pass && sel_eu() );
    else if( selectionTag == TString("ee_singleleptrig") ) return( singlelep_triggers_pass && sel_ee() );
    else if( selectionTag == TString("uu_singleleptrig") ) return( singlelep_triggers_pass && sel_uu() );
    else if( selectionTag == TString("eu_singleleptrig") ) return( singlelep_triggers_pass && sel_eu() );

    // In these selections we use a higher cut on the DIL mass: m(ll) > 60 GeV.
    // This is helpful for data-MC comparisons because we don't have low mass Z+jets samples
    else if( selectionTag == TString("ee_highmll") ) return( triggers_pass && sel_ee() && map_float["good_Mee"] > 60000 );
    else if( selectionTag == TString("uu_highmll") ) return( triggers_pass && sel_uu() && map_float["good_Muu"] > 60000 );

    // Selections that include the Z window
    else if( selectionTag == TString("ee_withz") ) return( triggers_pass && sel_ee_withz() );
    else if( selectionTag == TString("uu_withz") ) return( triggers_pass && sel_uu_withz() );

    else {
      report::error( "Missing definition for selection %s" , selectionTag.Data() );
      assert( false );
    }

  }

  static TString getSelectionTitle( const TString &selTag ) {
    if( selTag == TString("none") ) return "No selection";
    else if( selTag.BeginsWith("ee") ) return "ee";
    else if( selTag.BeginsWith("uu") ) return "#mu#mu";
    else if( selTag.BeginsWith("eu") ) return "e#mu";
    else return selTag;
  }

  // pin the jet requirements on the tail end of the selection title
  static TString getDressedSelectionTitle( const TString &selTag , const UInt_t &nj , const UInt_t &nb ) {
    return TString::Format("%s, #geq%uj, #geq%ub",getSelectionTitle(selTag).Data(),nj,nb);
  }

  TString getSelectionTitle() { return getSelectionTitle(selectionTag); }
  TString getDressedSelectionTitle() { return getDressedSelectionTitle(selectionTag,minNjets,minNbtags); }
  
  Double_t leadingLeptonPt() {
    if( good_el.size()==0 && good_mu.size()==0 ) return -666;
    if( good_el.size()==0 ) return good_mu[0].Pt();
    if( good_mu.size()==0 ) return good_el[0].Pt();
    return good_el[0].Pt() > good_mu[0].Pt() ? good_el[0].Pt() : good_mu[0].Pt();
  }

  void debugInfo() {
    report::blank( "" );
    for( std::size_t i : good_mu_id ) {
      report::debug( "Muon %i : true_pdg = %i" , i , map_vector_int["mu_true_pdg"]->at(i) );
    }
    report::blank( "" );
  }
  
};

#endif
