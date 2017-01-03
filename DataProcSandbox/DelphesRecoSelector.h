#ifndef _DELPHESRECOSELECTOR_H_
#define _DELPHESRECOSELECTOR_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <vector>

#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "Report.h"
#include "Particle.h"
#include "TopDecay.h"
#include "DelphesBtagger.h"

// Delphes includes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

class DelphesRecoSelector
{

private:
  // variables

  // User set jet multiplicity requirements
  Int_t minJets;
  Int_t maxJets;
  
  // vectors for keeping track of classified truth particles
  std::vector<RecoParticle*> v_all;
  std::vector<RecoParticle*> v_jets;
  std::vector<RecoParticle*> v_ljets;
  std::vector<RecoParticle*> v_bjets;
  std::vector<RecoParticle*> v_el;
  std::vector<RecoParticle*> v_mu;
  std::vector<RecoParticle*> v_lep;
  std::vector<RecoParticle*> v_met;
  std::vector<RecoParticle*> v_nu;

  // Keep track of whether we found at least one good
  // neutrino z-momentum solution in the particular event
  Bool_t nuMomentumSolved;
  
  // Readers for Delphes trees
  ExRootTreeReader *ex;
  TClonesArray *br_el;
  TClonesArray *br_mu;
  TClonesArray *br_jet;
  TClonesArray *br_met;
  DelphesBtagger *bt;

private:
  // functions

  void cleanup() {
    // Remove all truth particle currently in the record
    for( std::size_t i = 0 ; i < v_all.size() ; i++ ) delete v_all[i];
    v_all.clear();
    v_jets.clear();
    v_ljets.clear();
    v_bjets.clear();
    v_el.clear();
    v_mu.clear();
    v_lep.clear();
    v_met.clear();
    v_nu.clear();
    nuMomentumSolved = kTRUE;
  }

public:

  DelphesRecoSelector( ExRootTreeReader* _ex , DelphesBtagger *_bt , Int_t _minJets = -1 , Int_t _maxJets = -1 )
    : ex( _ex )
    , bt( _bt )
    , minJets( _minJets )
    , maxJets( _maxJets )
  {
    // Set up the Delphes reader branches
    br_el  = ex->UseBranch( "Electron" );
    br_mu  = ex->UseBranch( "Muon" );
    br_jet = ex->UseBranch( "Jet" );
    br_met = ex->UseBranch( "MissingET" );
  }

  ~DelphesRecoSelector() { cleanup(); }


  //
  // Basic return functions
  //
  std::size_t getN() { return v_all.size(); }
  std::size_t getNjets() { return v_jets.size(); }
  std::size_t getNljets() { return v_ljets.size(); }
  std::size_t getNbjets() { return v_bjets.size(); }
  std::size_t getNel() { return v_el.size(); }
  std::size_t getNmu() { return v_mu.size(); }
  std::size_t getNlep() { return v_el.size() + v_mu.size(); }

  std::size_t getNjetsPtAbove( double pt ) {
    std::size_t res = 0;
    for( std::size_t i = 0 ; i < v_jets.size() ; i++ ) {
      if( v_jets[i]->getPt() >= pt ) res++;
    }
    return res;
  }
  
  std::size_t getNbjets( int wp ) {
    std::size_t res = 0;
    for( std::size_t i = 0 ; i < v_bjets.size() ; i++ ) {
      if( v_bjets[i]->getTagLevel() >= wp ) res++;
    }
    return res;
  }

  Bool_t getNuMomentumSolved() { return nuMomentumSolved; }
  
  std::vector<RecoParticle*> getAll() { return v_all; }
  std::vector<RecoParticle*> getJets() { return v_jets; }
  std::vector<RecoParticle*> getLJets() { return v_ljets; }
  std::vector<RecoParticle*> getBJets() { return v_bjets; }
  std::vector<RecoParticle*> getMu() { return v_mu; }
  std::vector<RecoParticle*> getEl() { return v_el; }
  std::vector<RecoParticle*> getLep() { return v_lep; }

  RecoParticle* getParticle( std::size_t i ) { return v_all[i]; }
  RecoParticle* getJet( std::size_t i ) { return v_jets[i]; }
  RecoParticle* getLJet( std::size_t i ) { return v_ljets[i]; }
  RecoParticle* getBJet( std::size_t i ) { return v_bjets[i]; }
  RecoParticle* getMu( std::size_t i ) { return v_mu[i]; }
  RecoParticle* getEl( std::size_t i ) { return v_el[i]; }
  RecoParticle* getLep( std::size_t i ) { return v_lep[i]; }
  RecoParticle* getMet() { return v_met[0]; }
  RecoParticle* getNu( std::size_t i ) { return v_nu[i]; }

  //
  // Functions for constructing observables for combined objects
  // representing full system
  //

  Double_t getHtHad() {
    Double_t res = 0;
    for( std::size_t i = 0 ; i < v_jets.size() ; i++ ) res += v_jets[i]->getPt();
    return res;
  }

  Double_t getHtAll() {
    Double_t res = getHtHad();
    for( std::size_t i = 0 ; i < v_lep.size() ; i++ ) res += v_lep[i]->getPt();
    res += v_met[0]->getPt();
    return res;
  }
    
  
  

  //
  // Real nuts & bolts function used to process each event and
  // create a truth record that we can use to match reconstructed
  // objects
  //
  Bool_t processRecoRecord() {

    // Load event number "iev" in the tree
    //ex->ReadEntry( iev ); This is handled in DelphesTtbar.cpp
    
    // Make sure we clear vectors and reset indicators from previous event
    cleanup();

    //
    // Loop over particles in the record and keep the ones that pass minimum energy/eta requirements
    //

    // Jet selection
    for( std::size_t i = 0 ; i < br_jet->GetEntriesFast() ; ++i ) {

      Jet *jet = (Jet*) br_jet->At(i);

      if( jet->PT < 20 ) continue;
      if( TMath::Abs(jet->Eta) > 2.5 ) continue;

      Int_t tagLevel = bt->getTagLevel( jet->PT , jet->Flavor );

      v_all.push_back( new RecoParticle( jet->PT , jet->Eta , jet->Phi , jet->Mass , RecoParticle::JET , tagLevel ) );
      v_jets.push_back( v_all.back() );
      if( tagLevel > 0 )
	v_bjets.push_back( v_all.back() );
      else
	v_ljets.push_back( v_all.back() );

    }

    // Cut on min Njets
    if( int(v_jets.size()) < minJets ) return kFALSE;
    if( int(v_jets.size()) > maxJets && maxJets > 0 ) return kFALSE;
    
    // Electron selection
    for( std::size_t i = 0 ; i < br_el->GetEntriesFast() ; ++i ) {

      Electron *el = (Electron*) br_el->At(i);

      if( el->PT < 20 ) continue;
      if( TMath::Abs(el->Eta) > 2.5 ) continue;

      v_all.push_back( new RecoParticle( el->PT , el->Eta , el->Phi , 0.0 , RecoParticle::EL ) );
      v_el.push_back( v_all.back() );
      v_lep.push_back( v_all.back() );

    }

    // Muon selection
    for( std::size_t i = 0 ; i < br_mu->GetEntriesFast() ; ++i ) {

      Muon *mu = (Muon*) br_mu->At(i);

      if( mu->PT < 20 ) continue;
      if( TMath::Abs(mu->Eta) > 2.5 ) continue;

      v_all.push_back( new RecoParticle( mu->PT , mu->Eta , mu->Phi , 0.0 , RecoParticle::MU ) );
      v_mu.push_back( v_all.back() );
      v_lep.push_back( v_all.back() );

    }

    // Cut on the number of leptons
    //if( v_lep.size() != 1 ) return kFALSE;

    // Load the MET
    MissingET *met = (MissingET*) br_met->At(0);
    v_all.push_back( new RecoParticle( met->MET , 0.0 , met->Phi , 0.0 , RecoParticle::MET ) );
    v_met.push_back( v_all.back() );

    // If this is a single lepton event, calculate the z momentum of the neutrino by fixing
    // the W mass of the lepton+nu system. There are two possible solutions.
    if( v_lep.size() == 1 ) {

      Double_t mW  = 80.399;
      Double_t El  = v_lep[0]->getE();
      Double_t plx = v_lep[0]->getPx();
      Double_t ply = v_lep[0]->getPy();
      Double_t plz = v_lep[0]->getPz();
      Double_t pvx = v_met[0]->getPx();
      Double_t pvy = v_met[0]->getPy();

      // Solving quadratic
      Double_t psi = (0.5*mW*mW) + (plx*pvx) + (ply*pvy);
      Double_t A = (El*El) - (plz*plz);
      Double_t B = -2.0 * plz * psi;
      Double_t C = (El*El*pvx*pvx) + (El*El*pvy*pvy) - (psi*psi);

      Double_t pvz1 = 0.0;
      Double_t pvz2 = 0.0;

      if( ( (B*B) - (4.0*A*C) ) < 0 ) {
	// No solutions exist
	// Keep the z-momenta at zero
	//std::cout << std::endl;
	//report::warn( "No solutions found for z-momentum" );
	nuMomentumSolved = kFALSE;
      } else {
	pvz1 = ( -B + TMath::Sqrt( (B*B) - (4.0*A*C) ) ) / ( 2.0*A );
	pvz2 = ( -B - TMath::Sqrt( (B*B) - (4.0*A*C) ) ) / ( 2.0*A );
	// Order the two solutions by decreasing z-momentum magnitude
	if( TMath::Abs(pvz2) > TMath::Abs(pvz1) ) std::swap( pvz1 , pvz2 );
      }

      Double_t Ev1 = TMath::Sqrt( (pvx*pvx) + (pvy*pvy) + (pvz1*pvz1) );
      Double_t Ev2 = TMath::Sqrt( (pvx*pvx) + (pvy*pvy) + (pvz2*pvz2) );

      Double_t eta1 = 0.5 * TMath::Log( (Ev1 + pvz1) / (Ev1 - pvz1) );
      Double_t eta2 = 0.5 * TMath::Log( (Ev2 + pvz2) / (Ev2 - pvz2) );

      v_nu.push_back( new RecoParticle( v_met[0]->getPt() , eta1 , v_met[0]->getPhi() , 0.0 , RecoParticle::NU ) );
      v_nu.push_back( new RecoParticle( v_met[0]->getPt() , eta2 , v_met[0]->getPhi() , 0.0 , RecoParticle::NU ) );

      //std::cout << std::endl;
      //report::debug( "mW = %g , El = %g , plx = %g , ply = %g , plz = %g , pvx = %g , pvy = %g  " , mW , El , plx , ply , plz , pvx , pvy );
      //report::debug( "psi = %g , A = %g , B = %g , C = %g" , psi , A , B , C );
      //report::debug( "(B*B) - (4.0*A*C) = %g" , (B*B) - (4.0*A*C) );
      //report::debug( "Solution #1 : pvz = %10g , eta = %10g , MW = %10g" , pvz1 , eta1 , (v_nu[0]->getV4()+v_lep[0]->getV4()).M() );
      //report::debug( "Solution #2 : pvz = %10g , eta = %10g , MW = %10g" , pvz2 , eta2 , (v_nu[1]->getV4()+v_lep[0]->getV4()).M() );

    }

    // Success!
    return kTRUE;

  }

  // Return true if the reconstruction record passes the user-defined event selection
  Bool_t passesSelection( TString sel ) {
    if( int(v_jets.size()) < minJets ) return kFALSE;
    if( int(v_jets.size()) > maxJets && maxJets > 0 ) return kFALSE;
    if( sel==TString("ljet") ) {
      if( v_lep.size() != 1 ) return kFALSE;
      if( v_jets.size() < 4 ) return kFALSE;
      if( v_bjets.size() < 2 ) return kFALSE;
      //if( met->getPt() < 20 ) return kFALSE;
    } else if( sel==TString("ejet") ) {
      if( v_el.size() != 1 ) return kFALSE;
      if( v_mu.size() != 0 ) return kFALSE;
      if( v_jets.size() < 4 ) return kFALSE;
      if( v_bjets.size() < 2 ) return kFALSE;
      //if( met->getPt() < 20 ) return kFALSE;
    } else if( sel==TString("mjet") ) {
      if( v_el.size() != 0 ) return kFALSE;
      if( v_mu.size() != 1 ) return kFALSE;
      if( v_jets.size() < 4 ) return kFALSE;
      if( v_bjets.size() < 2 ) return kFALSE;
      //if( met->getPt() < 20 ) return kFALSE;
    } else if( sel==TString("dil") ) {
      if( v_lep.size() != 2 ) return kFALSE;
      if( v_jets.size() < 2 ) return kFALSE;
      if( v_bjets.size() < 2 ) return kFALSE;
      //if( met->getPt() < 20 ) return kFALSE;
    }
    return kTRUE;
  }

};

#endif
