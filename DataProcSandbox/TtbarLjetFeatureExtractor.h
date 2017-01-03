#ifndef _TTBARLJETFEATUREEXTRACTOR_H_
#define _TTBARLJETFEATUREEXTRACTOR_H_

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


class
TtbarLjetFeatureExtractor
{

  //
  // For now I assume that we are targetting the semi-leptonic ttbar decay
  //

private:

  DelphesRecoSelector *rSel;
  RecoParticle *lepTopJet;
  RecoParticle *hadTopJet;
  RecoParticle *hadWJet1;
  RecoParticle *hadWJet2;
  RecoParticle *lep;
  RecoParticle *nuSol1;
  RecoParticle *nuSol2;

  std::map<TString,Double_t> m_features;
  std::vector<TString> v_featureNames;

  std::vector<TString> v_eventCollectionNames;
  std::vector<TString> v_extendedCollectionNames;
  std::vector<TString> v_shapeNames;
  std::vector< std::pair<TString,TString> > v_pairNames;

  TString outputPath;
  std::ofstream *outputCsv;

  void addFeature( TString s ) { v_featureNames.push_back(s); }


public:

  TtbarLjetFeatureExtractor() {

    // Multi-particle systems to use when building event features
    v_eventCollectionNames =  {
      "lepWSol1" ,
      "lepWSol2"
    };
			       
    // Multi-particle systems to use when building extended features
    v_extendedCollectionNames = {
      "hadW" ,
      //"lepWSol1" ,  // This collection is included with event-wide features
      //"lepWSol2" ,  // This collection is included with event-wide features
      "hadTop" ,
      "lepTopSol1" ,
      "lepTopSol2" ,
      "ttbarSol1" ,
      "ttbarSol2"
    };

    // Names of shape parameters to use when creating object or multi-object features
    v_shapeNames = {
      "Aplanarity" ,
      "AplanarityO" ,
      "Sphericity" ,
      "SphericityO" ,
      "SphericityT" ,
      "Planarity" ,
      "VariableC" ,
      "VariableD" ,
      "Circularity" ,
      "PlanarFlow"
    };
    
    // Object pairs to use when constructing deta/dphi/dr features
    v_pairNames = {
      // hadTop separations
      std::make_pair( "hadTop" , "lepTopSol1" ) ,
      std::make_pair( "hadTop" , "lepTopSol2" ) ,
      std::make_pair( "hadTop" , "hadW" ) ,
      std::make_pair( "hadTop" , "lepWSol1" ) ,
      std::make_pair( "hadTop" , "lepWSol2" ) ,
      std::make_pair( "hadTop" , "hadTopJet" ) ,
      std::make_pair( "hadTop" , "hadWJet1" ) ,
      std::make_pair( "hadTop" , "hadWJet2" ) ,
      std::make_pair( "hadTop" , "lepTopJet" ) ,
      std::make_pair( "hadTop" , "lep" ) ,
      std::make_pair( "hadTop" , "nuSol1" ) ,
      std::make_pair( "hadTop" , "nuSol2" ) ,
      // hadW separations
      std::make_pair( "hadW" , "lepTopSol1" ) ,
      std::make_pair( "hadW" , "lepTopSol2" ) ,
      std::make_pair( "hadW" , "lepWSol1" ) ,
      std::make_pair( "hadW" , "lepWSol2" ) ,
      std::make_pair( "hadW" , "hadTopJet" ) ,
      std::make_pair( "hadW" , "hadWJet1" ) ,
      std::make_pair( "hadW" , "hadWJet2" ) ,
      std::make_pair( "hadW" , "lepTopJet" ) ,
      std::make_pair( "hadW" , "lep" ) ,
      std::make_pair( "hadW" , "nuSol1" ) ,
      std::make_pair( "hadW" , "nuSol2" ) ,
      // hadTopJet separations
      std::make_pair( "hadTopJet" , "lepTopSol1" ) ,
      std::make_pair( "hadTopJet" , "lepTopSol2" ) ,
      std::make_pair( "hadTopJet" , "lepWSol1" ) ,
      std::make_pair( "hadTopJet" , "lepWSol2" ) ,
      std::make_pair( "hadTopJet" , "hadWJet1" ) ,
      std::make_pair( "hadTopJet" , "hadWJet2" ) ,
      std::make_pair( "hadTopJet" , "lepTopJet" ) ,
      std::make_pair( "hadTopJet" , "lep" ) ,
      std::make_pair( "hadTopJet" , "nuSol1" ) ,
      std::make_pair( "hadTopJet" , "nuSol2" ) ,
      // hadWJet1 seprations
      std::make_pair( "hadWJet1" , "lepTopSol1" ) ,
      std::make_pair( "hadWJet1" , "lepTopSol2" ) ,
      std::make_pair( "hadWJet1" , "lepWSol1" ) ,
      std::make_pair( "hadWJet1" , "lepWSol2" ) ,
      std::make_pair( "hadWJet1" , "hadWJet2" ) ,
      std::make_pair( "hadWJet1" , "lepTopJet" ) ,
      std::make_pair( "hadWJet1" , "lep" ) ,
      std::make_pair( "hadWJet1" , "nuSol1" ) ,
      std::make_pair( "hadWJet1" , "nuSol2" ) ,      
      // hadWJet2 separations
      std::make_pair( "hadWJet2" , "lepTopSol1" ) ,
      std::make_pair( "hadWJet2" , "lepTopSol2" ) ,
      std::make_pair( "hadWJet2" , "lepWSol1" ) ,
      std::make_pair( "hadWJet2" , "lepWSol2" ) ,
      std::make_pair( "hadWJet2" , "lepTopJet" ) ,
      std::make_pair( "hadWJet2" , "lep" ) ,
      std::make_pair( "hadWJet2" , "nuSol1" ) ,
      std::make_pair( "hadWJet2" , "nuSol2" ) ,      
      // lepTopSol separations
      std::make_pair( "lepTopSol1" , "lepWSol1" ) ,
      std::make_pair( "lepTopSol1" , "lepWSol2" ) ,
      std::make_pair( "lepWSol1" , "lepTopSol2" ) ,
      std::make_pair( "lepWSol2" , "lepTopSol2" ) ,
      std::make_pair( "lepTopJet" , "lepTopSol1" ) ,
      std::make_pair( "lepTopJet" , "lepTopSol2" ) ,
      std::make_pair( "lep" , "lepTopSol1" ) ,
      std::make_pair( "lep" , "lepTopSol2" ) ,
      std::make_pair( "lepTopSol1" , "nuSol1" ) ,
      std::make_pair( "lepTopSol1" , "nuSol2" ) ,
      std::make_pair( "nuSol1" , "lepTopSol2" ) ,
      std::make_pair( "nuSol2" , "lepTopSol2" ) ,
      // lepWSol separations
      std::make_pair( "lepTopJet" , "lepWSol1" ) ,
      std::make_pair( "lepTopJet" , "lepWSol2" ) ,
      std::make_pair( "lep" , "lepWSol1" ) ,
      std::make_pair( "lep" , "lepWSol2" ) ,
      std::make_pair( "lepWSol1" , "nuSol1" ) ,
      std::make_pair( "lepWSol1" , "nuSol2" ) ,
      std::make_pair( "nuSol1" , "lepWSol2" ) ,
      std::make_pair( "nuSol2" , "lepWSol2" ) ,
      // lepTopJet separations
      std::make_pair( "lep" , "lepTopJet" ) ,
      std::make_pair( "lepTopJet" , "nuSol1" ) ,
      std::make_pair( "lepTopJet" , "nuSol2" ) ,
      // lep separations
      std::make_pair( "lep" , "nuSol1" ) 
      // std::make_pair( "lep" , "nuSol2" ) // dR is 100% correlated with dR(lep,nuSol1)
    };
    
    // Bookkeeping features
    addFeature( "EventId" );
    addFeature( "CombId" );

    // Event-wide features
    addFeature( "nJets" );
    for( double pt : { 30 , 40 , 50 } ) addFeature( TString::Format("nJetsPtAbove%i",int(pt)) );
    for( int i = 1 ; i <= 5 ; i++ ) addFeature( TString::Format("nBtags%i",i) );
    addFeature( "leptonIsMuon" );
    addFeature( "nuMomentumSolved" );
    addFeature( "pt_lep" );
    addFeature( "eta_lep" );
    addFeature( "phi_lep" );
    addFeature( "pt_nuSol1" );
    addFeature( "eta_nuSol1" );
    addFeature( "eta_nuSol2" );
    addFeature( "phi_nuSol1" );
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
    for( TString n : v_eventCollectionNames ) {
      addFeature( "eta_"+n );
      if( ! n==TString("lepWSol2") ) {
	addFeature( "m_"+n );
      }
      if( ! n.Contains("Sol2") ) {
	addFeature( "pt_"+n );
	addFeature( "phi_"+n );
	addFeature( "mt_"+n );
	addFeature( "ptsum_"+n );
      }
      for( TString s : v_shapeNames ) {
	// For 2 object systems (i.e. reconstructed W candidates) some of the shape
	// values are 100% correlated:
	//    Sphericity =~= Planarity =~= Circularity, and
	//    VariableC =~= PlanarFlow
	if( s==TString("Planarity") ) continue;
	if( s==TString("Circularity") ) continue;
	if( s==TString("PlanarFlow") ) continue;
	// Furthermore, SphericityT is independent of neutrino pz solution
	if( n.Contains("Sol2") && s==TString("SphericityT") ) continue;
	addFeature( s+"_"+n );
      }
    }

    // Combination-specific basic features
    for( int i = 1 ; i <= 5 ; i++ ) addFeature( TString::Format("nBtags%i_ttbar",i) );
    for( int i = 1 ; i <= 5 ; i++ ) addFeature( TString::Format("nBtags%i_ttbarDecay",i) );
    addFeature( "pt_lepTopJet" );
    addFeature( "eta_lepTopJet" );
    addFeature( "phi_lepTopJet" );
    addFeature( "m_lepTopJet" );
    addFeature( "wp_lepTopJet" );
    addFeature( "pt_hadTopJet" );
    addFeature( "eta_hadTopJet" );
    addFeature( "phi_hadTopJet" );
    addFeature( "m_hadTopJet" );
    addFeature( "wp_hadTopJet" );
    addFeature( "pt_hadWJet1" );
    addFeature( "eta_hadWJet1" );
    addFeature( "phi_hadWJet1" );
    addFeature( "m_hadWJet1" );
    addFeature( "wp_hadWJet1" );
    addFeature( "pt_hadWJet2" );
    addFeature( "eta_hadWJet2" );
    addFeature( "phi_hadWJet2" );
    addFeature( "m_hadWJet2" );
    addFeature( "wp_hadWJet2" );


    // Combination-specific extended features
    for( TString n : v_extendedCollectionNames ) {
      addFeature( "eta_"+n );
      addFeature( "m_"+n );
      if( ! n.Contains("Sol2") ) {
	addFeature( "pt_"+n );
	addFeature( "phi_"+n );
	addFeature( "ptsum_"+n );
      }
      if( n.Contains("had") || (n.Contains("ttbar") && !n.Contains("Sol2")) ) {
	addFeature( "wpsum_"+n );
      }
      for( TString s : v_shapeNames ) {
	// For 2 object systems (i.e. reconstructed W candidates) some of the shape
	// values are 100% correlated:
	//    Sphericity =~= Planarity =~= Circularity, and
	//    VariableC =~= PlanarFlow
	if( n==TString("hadW") ) {
	  if( s==TString("Planarity") ) continue;
	  if( s==TString("Circularity") ) continue;
	  if( s==TString("PlanarFlow") ) continue;
	}
	// Furthermore, SphericityT is independent of neutrino pz solution
	if( n.Contains("Sol2") && s==TString("SphericityT") ) continue;
	addFeature( s+"_"+n );
      }
    }

    // Further extended features
    for( std::pair<TString,TString> p : v_pairNames ) {
      if( ! p.second.Contains("Sol2") ) {
	addFeature( "dphi_"+p.first+"_"+p.second );
      }
      addFeature( "deta_"+p.first+"_"+p.second );
      addFeature( "dr_"+p.first+"_"+p.second );
    }
    
    // Training output
    addFeature( "signal" );
  }

  ~TtbarLjetFeatureExtractor() {}

  //
  // Basic functions
  //

  void setRecoSelector( DelphesRecoSelector *_rSel ) { rSel = _rSel; }
  Int_t getNbtagsTtbarDecay( int wp ) {
    Int_t res = 0;
    if( lepTopJet->getTagLevel() >= wp ) res++;
    if( hadTopJet->getTagLevel() >= wp ) res++;
    return res;
  }
  Int_t getNbtagsTtbarSystem( int wp ) {
    Int_t res = getNbtagsTtbarDecay( wp );
    if( hadWJet1->getTagLevel() >= wp ) res++;
    if( hadWJet2->getTagLevel() >= wp ) res++;
    return res;
  }
  Int_t getWPSumTtbarSystem() {
    Int_t res = 0;
    res += lepTopJet->getTagLevel();
    res += hadTopJet->getTagLevel();
    res += hadWJet1->getTagLevel();
    res += hadWJet2->getTagLevel();
    return res;
  }



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

  void fill( Int_t eventId , Int_t combId , RecoParticle *_lepTopJet , RecoParticle *_hadTopJet , RecoParticle *_hadWJet1 , RecoParticle *_hadWJet2 ) {

    lepTopJet = _lepTopJet;
    hadTopJet = _hadTopJet;
    hadWJet1  = _hadWJet1;
    hadWJet2  = _hadWJet2;
    lep       = rSel->getLep(0);
    nuSol1    = rSel->getNu(0);
    nuSol2    = rSel->getNu(1);

    // Build some collections for extended shape features
    std::map<TString,RecoParticleCollection> m_collections;
    m_collections["hadW"]       = RecoParticleCollection( { hadWJet1 , hadWJet2 } );
    m_collections["lepWTrans"]  = RecoParticleCollection( { lep , rSel->getMet() } );
    m_collections["lepWSol1"]   = RecoParticleCollection( { lep , nuSol1 } );
    m_collections["lepWSol2"]   = RecoParticleCollection( { lep , nuSol2 } );
    m_collections["hadTop"]     = RecoParticleCollection( { hadTopJet , hadWJet1 , hadWJet2 } );
    m_collections["lepTopSol1"] = RecoParticleCollection( { lepTopJet , lep , nuSol1 } );
    m_collections["lepTopSol2"] = RecoParticleCollection( { lepTopJet , lep , nuSol2 } );
    m_collections["ttbarSol1"]  = RecoParticleCollection( { hadTopJet , hadWJet1 , hadWJet2 , lepTopJet , lep , nuSol1 } );
    m_collections["ttbarSol2"]  = RecoParticleCollection( { hadTopJet , hadWJet1 , hadWJet2 , lepTopJet , lep , nuSol2 } );
    
    std::map<TString,TLorentzVector> m_allV4;
    m_allV4["lepTopJet"]  = lepTopJet->getV4();
    m_allV4["hadTopJet"]  = hadTopJet->getV4();
    m_allV4["hadWJet1"]	  = hadWJet1->getV4();
    m_allV4["hadWJet2"]	  = hadWJet2->getV4();
    m_allV4["lep"]	  = lep->getV4();
    m_allV4["nuSol1"]	  = nuSol1->getV4();
    m_allV4["nuSol2"]	  = nuSol2->getV4();
    m_allV4["hadW"]	  = m_collections["hadW"].getV4();
    m_allV4["lepWSol1"]	  = m_collections["lepWSol1"].getV4();
    m_allV4["lepWSol2"]	  = m_collections["lepWSol2"].getV4();
    m_allV4["hadTop"]	  = m_collections["hadTop"].getV4();
    m_allV4["lepTopSol1"] = m_collections["lepTopSol1"].getV4();
    m_allV4["lepTopSol2"] = m_collections["lepTopSol2"].getV4();
    m_allV4["ttbarSol1"]  = m_collections["ttbarSol1"].getV4();
    m_allV4["ttbarSol2"]  = m_collections["ttbarSol2"].getV4();
    
    m_features["EventId"] = eventId;
    m_features["CombId"]  = combId;

    if( combId == 0 ) {

      // This is the first fill of a new event
      // so we should reset all the event-wide variables that are fixed for all
      // combinations in a given event, e.g. Njets, Nbtags...

      m_features["nJets"] = rSel->getNjets();
      for( double pt : { 30 , 40 , 50 } ) m_features[TString::Format("nJetsPtAbove%i",int(pt))] = rSel->getNjetsPtAbove(pt);
      for( int i = 1 ; i <= 5 ; i++ ) m_features[TString::Format("nBtags%i",i)] = rSel->getNbjets(i);

      m_features["leptonIsMuon"] = ( lep->getType()==RecoParticle::MU ? 1 : 0 );
      m_features["nuMomentumSolved"] = ( rSel->getNuMomentumSolved() ? 1 : 0 );
      m_features["pt_lep"]  = lep->getPt();
      m_features["eta_lep"] = lep->getEta();
      m_features["phi_lep"] = lep->getPhi();
      m_features["pt_nuSol1"]  = nuSol1->getPt();
      m_features["eta_nuSol1"] = nuSol1->getEta();
      m_features["eta_nuSol2"] = nuSol2->getEta();
      m_features["phi_nuSol1"] = nuSol1->getPhi();
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
      m_features["HT_all"]     = rSel->getHtAll();
      m_features["HT_had"]     = rSel->getHtHad();

      for( TString n : v_eventCollectionNames ) {
	m_features["eta_"+n]   = m_collections[n].getEta();
	m_features["m_"+n]     = m_collections[n].getM();
	m_features["pt_"+n]    = m_collections[n].getPt();
	m_features["phi_"+n]   = m_collections[n].getPhi();
	m_features["mt_"+n]    = m_collections[n].getMt();
	m_features["ptsum_"+n] = m_collections[n].getPtsum();
	for( TString s : v_shapeNames )
	  m_features[s+"_"+n]  = m_collections[n].getShape(s);
      }

    }

    // basic features
    for( int i = 0 ; i <= 5 ; i++ ) m_features[TString::Format("nBtags%i_ttbar",i)] = Double_t(getNbtagsTtbarSystem(i));
    for( int i = 0 ; i <= 5 ; i++ ) m_features[TString::Format("nBtags%i_ttbarDecay",i)] = Double_t(getNbtagsTtbarDecay(i));
    m_features["pt_lepTopJet"]  = lepTopJet->getPt();
    m_features["eta_lepTopJet"] = lepTopJet->getEta();
    m_features["phi_lepTopJet"] = lepTopJet->getPhi();
    m_features["m_lepTopJet"]   = lepTopJet->getM();
    m_features["wp_lepTopJet"]  = lepTopJet->getTagLevel();
    m_features["pt_hadTopJet"]  = hadTopJet->getPt();
    m_features["eta_hadTopJet"] = hadTopJet->getEta();
    m_features["phi_hadTopJet"] = hadTopJet->getPhi();
    m_features["m_hadTopJet"]   = hadTopJet->getM();
    m_features["wp_hadTopJet"]  = hadTopJet->getTagLevel();
    m_features["pt_hadWJet1"]   = hadWJet1->getPt();
    m_features["eta_hadWJet1"]  = hadWJet1->getEta();
    m_features["phi_hadWJet1"]  = hadWJet1->getPhi();
    m_features["m_hadWJet1"]    = hadWJet1->getM();
    m_features["wp_hadWJet1"]   = hadWJet1->getTagLevel();
    m_features["pt_hadWJet2"]   = hadWJet2->getPt();
    m_features["eta_hadWJet2"]  = hadWJet2->getEta();
    m_features["phi_hadWJet2"]  = hadWJet2->getPhi();
    m_features["m_hadWJet2"]    = hadWJet2->getM();
    m_features["wp_hadWJet2"]   = hadWJet2->getTagLevel();

    // extended features
    for( TString n : v_extendedCollectionNames ) {
      m_features["eta_"+n]   = m_collections[n].getEta();
      m_features["m_"+n]     = m_collections[n].getM();
      m_features["pt_"+n]    = m_collections[n].getPt();
      m_features["phi_"+n]   = m_collections[n].getPhi();
      m_features["wpsum_"+n] = m_collections[n].getWpsum();
      m_features["ptsum_"+n] = m_collections[n].getPtsum();
      for( TString s : v_shapeNames )
	m_features[s+"_"+n]  = m_collections[n].getShape(s);
    }
    
    // further extended features
    for( std::pair<TString,TString> p : v_pairNames ) {
      m_features["dphi_"+p.first+"_"+p.second] = m_allV4[p.first].DeltaPhi( m_allV4[p.second] );
      m_features["deta_"+p.first+"_"+p.second] = m_allV4[p.first].Eta() - m_allV4[p.second].Eta();
      m_features["dr_"+p.first+"_"+p.second]   = m_allV4[p.first].DeltaR( m_allV4[p.second] );;
    }

    // check if the combination is properly matched to a ttbar decay
    Bool_t hadTopMatched = hadTopJet->fromCommonWofSameTop( hadWJet1 , hadWJet2 );
    Bool_t lepTopMatched = lepTopJet->fromCommonTop( lep );
    m_features["signal"] = hadTopMatched && lepTopMatched ? 1 : 0;

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
