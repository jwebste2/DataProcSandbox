#ifndef _TTBARFEATUREEXTRACTOR_H_
#define _TTBARFEATUREEXTRACTOR_H_

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
#include "TClonesArray.h"

#include "Report.h"
#include "Particle.h"
#include "RecoParticleCollection.h"
#include "DelphesRecoSelector.h"
#include "DelphesTruthSelector.h"


class
TtbarFeatureExtractor
{

  //
  // For now I assume that we are targetting the semi-leptonic ttbar decay
  //

private:

  DelphesRecoSelector *rSel;
  DelphesTruthSelector *tSel;
  
  std::map<TString,Double_t> m_features;
  std::vector<TString> v_featureNames;

  std::vector<TString> v_collectionNames;
  std::vector< std::pair<TString,TString> > v_pairNames;
  
  TString outputPath;
  std::ofstream *outputCsv;

  void addFeature( TString s ) { v_featureNames.push_back(s); }

public:

  TtbarFeatureExtractor()
    : rSel( (DelphesRecoSelector*)0 )
    , tSel( (DelphesTruthSelector*)0 )
  {

    // Bookkeeping features
    addFeature( "EventId" );
    addFeature( "CombId" );
    addFeature( "Tag" );

    // Event-wide features
    addFeature( "nJets" );
    for( double pt : { 30 , 40 , 50 } ) addFeature( TString::Format("nJetsPtAbove%i",int(pt)) );
    for( int i = 1 ; i <= 5 ; i++ ) addFeature( TString::Format("nBtags%i",i) );
    addFeature( "pt_met" );
    addFeature( "phi_met" );
    for( int i = 1 ; i <= 4 ; i++ ) {
      addFeature( TString::Format("pt_lep%i",i) );
      addFeature( TString::Format("eta_lep%i",i) );
      addFeature( TString::Format("phi_lep%i",i) );
      addFeature( TString::Format("isMuon_lep%i",i) );
    }
    for( int i = 1 ; i <= 6 ; i++ ) {
      addFeature( TString::Format("pt_jet%i",i) );
      addFeature( TString::Format("eta_jet%i",i) );
      addFeature( TString::Format("phi_jet%i",i) );
      addFeature( TString::Format("m_jet%i",i) );
      addFeature( TString::Format("wp_jet%i",i) );
    }
    for( int i = 1 ; i <= 3 ; i++ ) {
      addFeature( TString::Format("pt_bjet%i",i) );
      addFeature( TString::Format("eta_bjet%i",i) );
      addFeature( TString::Format("phi_bjet%i",i) );
      addFeature( TString::Format("m_bjet%i",i) );
      addFeature( TString::Format("wp_bjet%i",i) );
    }
    addFeature( "HT_all" );
    addFeature( "HT_had" );

    addFeature( "decayT" );
    addFeature( "decayTbar" );
    addFeature( "decayWp" );
    addFeature( "decayWm" );

  }

  ~TtbarFeatureExtractor() {}

  //
  // Basic functions
  //

  void setRecoSelector( DelphesRecoSelector *s ) { rSel = s; }
  void setTruthSelector( DelphesTruthSelector *s ) { tSel = s; }


  //
  // Functions for building CSV file
  //

  void openCsv( TString fname ) {
    outputPath = fname;
    outputCsv  = new std::ofstream( fname.Data() );
    (*outputCsv) << v_featureNames[0];
    for( std::size_t i = 1 ; i < v_featureNames.size() ; i++ ) {
      (*outputCsv) << "," << v_featureNames[i];
    }
    (*outputCsv) << std::endl;
  }
  void closeCsv() { outputCsv->close(); }

  void closeAndZipCsv() {
    outputCsv->close();
    report::info( "Bzipping %s" , outputPath.Data() );
    std::system( TString::Format( "bzip2 %s" , outputPath.Data() ).Data() );
  }
  
  void fill( Int_t eventId , Int_t combId , Int_t tag ) {

    m_features["EventId"] = eventId;
    m_features["CombId"]  = combId;
    m_features["Tag"]     = tag;

    if( combId <= 0 ) {

      // This is the first fill of a new event
      // so we should reset all the event-wide variables that are fixed for all
      // combinations in a given event, e.g. Njets, Nbtags...
      
      m_features["nJets"] = rSel->getNjets();
      for( double pt : { 30 , 40 , 50 } ) m_features[TString::Format("nJetsPtAbove%i",int(pt))] = rSel->getNjetsPtAbove(pt);
      for( int i = 1 ; i <= 5 ; i++ ) m_features[TString::Format("nBtags%i",i)] = rSel->getNbjets(i);

      m_features["pt_met"]  = rSel->getMet()->getPt();
      m_features["phi_met"] = rSel->getMet()->getPhi();
      for( int i = 1 ; i <= 4 ; i++ ) {
	if( rSel->getNlep() < i ) {
	  m_features[TString::Format("pt_lep%i",i)]     = -10;
	  m_features[TString::Format("eta_lep%i",i)]    = -10;
	  m_features[TString::Format("phi_lep%i",i)]    = -10;
	  m_features[TString::Format("isMuon_lep%i",i)] = -10;
	} else {
	  m_features[TString::Format("pt_lep%i",i)]     = rSel->getLep(i-1)->getPt();
	  m_features[TString::Format("eta_lep%i",i)]    = rSel->getLep(i-1)->getEta();
	  m_features[TString::Format("phi_lep%i",i)]    = rSel->getLep(i-1)->getPhi();
	  m_features[TString::Format("isMuon_lep%i",i)] = ( rSel->getLep(i-1)->getType()==RecoParticle::MU ? 1 : 0 );
	}
      }
      for( int i = 1 ; i <= 6 ; i++ ) {
        if( rSel->getNjets() < i ) {
          m_features[TString::Format("pt_jet%i",i)]  = -10;
          m_features[TString::Format("eta_jet%i",i)] = -10;
          m_features[TString::Format("phi_jet%i",i)] = -10;
          m_features[TString::Format("m_jet%i",i)]   = -10;
          m_features[TString::Format("wp_jet%i",i)]  = -10;
        } else {
          m_features[TString::Format("pt_jet%i",i)]  = rSel->getJet(i-1)->getPt();
          m_features[TString::Format("eta_jet%i",i)] = rSel->getJet(i-1)->getEta();
          m_features[TString::Format("phi_jet%i",i)] = rSel->getJet(i-1)->getPhi();
          m_features[TString::Format("m_jet%i",i)]   = rSel->getJet(i-1)->getM();
          m_features[TString::Format("wp_jet%i",i)]  = rSel->getJet(i-1)->getTagLevel();
        }
      }
      for( int i = 1 ; i <= 3 ; i++ ) {
        if( rSel->getNbjets() < i ) {
          m_features[TString::Format("pt_bjet%i",i)]  = -10;
          m_features[TString::Format("eta_bjet%i",i)] = -10;
          m_features[TString::Format("phi_bjet%i",i)] = -10;
          m_features[TString::Format("m_bjet%i",i)]   = -10;
          m_features[TString::Format("wp_bjet%i",i)]  = -10;
        } else {
          m_features[TString::Format("pt_bjet%i",i)]  = rSel->getBJet(i-1)->getPt();
          m_features[TString::Format("eta_bjet%i",i)] = rSel->getBJet(i-1)->getEta();
          m_features[TString::Format("phi_bjet%i",i)] = rSel->getBJet(i-1)->getPhi();
          m_features[TString::Format("m_bjet%i",i)]   = rSel->getBJet(i-1)->getM();
          m_features[TString::Format("wp_bjet%i",i)]  = rSel->getBJet(i-1)->getTagLevel();
        }
      }
      m_features["HT_all"] = rSel->getHtAll();
      m_features["HT_had"] = rSel->getHtHad();

      m_features["decayT"]    = tSel->getDecayT();
      m_features["decayTbar"] = tSel->getDecayTbar();
      m_features["decayWp"]   = tSel->getDecayWp();
      m_features["decayWm"]   = tSel->getDecayWm();
      
    }

  }

  void save() {
    (*outputCsv) << m_features[v_featureNames[0]];
    for( std::size_t i = 1 ; i < v_featureNames.size() ; i++ ) {
      (*outputCsv) << "," << m_features[v_featureNames[i]];
    }
    (*outputCsv) << std::endl;
  }

  void dump() {
    std::cout << std::endl;
    for( std::size_t i = 0 ; i < v_featureNames.size() ; i++ ) {
      std::cout << TString::Format("%25s",v_featureNames[i].Data()) << " = " << m_features[v_featureNames[i]] << std::endl;
    }
    std::cout << std::endl;
  }



};

#endif
