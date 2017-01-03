#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>

#include <TMath.h>
#include <TLorentzVector.h>

class Particle
{

protected:
  TLorentzVector v4;

public:
  Particle( Double_t _pt , Double_t _eta , Double_t _phi , Double_t _m ) {
    v4.SetPtEtaPhiM( _pt , _eta , _phi , _m );
  }
  ~Particle() {}
  
  TLorentzVector getV4() { return v4; }
  Double_t getPt() { return v4.Pt(); }
  Double_t getPx() { return v4.Px(); }
  Double_t getPy() { return v4.Py(); }
  Double_t getPz() { return v4.Pz(); }
  Double_t getEta() { return v4.Eta(); }
  Double_t getPhi() { return v4.Phi(); }
  Double_t getM() { return v4.M(); }
  Double_t getE() { return v4.E(); }

  Double_t deltaR( const TLorentzVector &other ) { return v4.DeltaR(other); }

  Double_t deltaR( Particle *other ) { return v4.DeltaR(other->getV4()); }
  Double_t deltaPhi( Particle *other ) { return v4.DeltaPhi(other->getV4()); }
  Double_t deltaEta( Particle *other ) { return v4.Eta() - other->getEta(); }

  // Function returns merged mass of a vector of particles
  static Double_t getMergedMass( std::vector<Particle*> v ) {
    if( v.size()==0 ) return 0.0;
    TLorentzVector l = v[0]->getV4();
    for( std::size_t i = 1 ; i < v.size() ; i++ ) l += v[i]->getV4();
    return l.M();
  }
  
};

class
TruthParticle : public Particle
{

private:
  Int_t parentPDG;
  Int_t PDG;
  
public:

  TruthParticle( Double_t _pt , Double_t _eta , Double_t _phi , Double_t _m , Int_t _PDG , Int_t _parentPDG )
    : Particle( _pt , _eta , _phi , _m )
    , PDG( _PDG )
    , parentPDG( _parentPDG )
  {}
  
  ~TruthParticle() {}

  Int_t getPDG() { return PDG; }
  Int_t getParentPDG() { return parentPDG; }
  
  Bool_t isFromWp() { return( parentPDG == 24 ); }
  Bool_t isFromWm() { return( parentPDG == -24 ); }
  Bool_t isFromW() { return( TMath::Abs(parentPDG) == 24 ); }
  Bool_t isFromTp() { return( parentPDG == 6 ); }
  Bool_t isFromTm() { return( parentPDG == -6 ); }
  Bool_t isFromT() { return( TMath::Abs(parentPDG) == 6 ); }

 
};


class
RecoParticle : public Particle 
{

public:
  enum PType {
    JET ,
    EL ,
    MU ,
    TAU ,
    PHOTON ,
    MET ,
    BJET ,
    LJET ,
    NU ,
    UNDEFINED
  };
  
private:
  PType type;
  Int_t tagLevel;
  TruthParticle *truth_match;

public:
  
  RecoParticle( Double_t _pt , Double_t _eta , Double_t _phi , Double_t _m , PType _type , Int_t _tagLevel = -1 )
    : Particle( _pt , _eta , _phi , _m )
    , type( _type )
    , tagLevel( _tagLevel )
    , truth_match( (TruthParticle*)0 )
  {}

  ~RecoParticle() {}

  PType getType() { return type; }
  TruthParticle* getTruthParticle() { return truth_match; }
  void setTagLevel( const Int_t &v ) { tagLevel = v; }
  Int_t getTagLevel() { return tagLevel; }
  void setTruthParticle( TruthParticle *p ) { truth_match = p; }
  Bool_t isTruthMatched() { return truth_match != 0; }
  Bool_t matchesPDG( Int_t PDG ) {
    Int_t AbsPDG = TMath::Abs(PDG);
    if( type==JET && AbsPDG>=1 && AbsPDG<=5 ) return kTRUE;
    if( type==EL && AbsPDG==11 ) return kTRUE;
    if( type==MU && AbsPDG==13 ) return kTRUE;
    if( type==TAU && AbsPDG==15 ) return kTRUE;
    if( type==PHOTON && PDG==22 ) return kTRUE;
    if( type==BJET && AbsPDG==5 ) return kTRUE;
    if( type==LJET && AbsPDG>=1 && AbsPDG<=4 ) return kTRUE;
    //if( type==NU && (AbsPDG==12 || AbsPDG==14 || AbsPDG==16) ) return kTRUE;
    return kFALSE;
  }

  //
  // Function returns true if the input reco particle is from the same top decay as
  // this particle.
  //
  Bool_t fromCommonTop( RecoParticle *other ) {
    if( ! isTruthMatched() ) return kFALSE;
    if( ! other->isTruthMatched() ) return kFALSE;
    if( getTruthParticle()->isFromWp() || getTruthParticle()->isFromTp() )
      return other->getTruthParticle()->isFromWp() || other->getTruthParticle()->isFromTp();
    if( getTruthParticle()->isFromWm() || getTruthParticle()->isFromTm() )
      return other->getTruthParticle()->isFromWm() || other->getTruthParticle()->isFromTm();
    return kFALSE;
  }

  //
  // Function returns true if the two reconstructed particles are from the same top decay
  // as this particle.
  //
  Bool_t fromCommonTop( RecoParticle *other1 , RecoParticle *other2 ) {
    if( ! isTruthMatched() ) return kFALSE;
    if( ! other1->isTruthMatched() ) return kFALSE;
    if( ! other2->isTruthMatched() ) return kFALSE;
    if( getTruthParticle()->isFromWp() || getTruthParticle()->isFromTp() ) {
      return ( (other1->getTruthParticle()->isFromWp() || other1->getTruthParticle()->isFromTp()) &&
	       (other2->getTruthParticle()->isFromWp() || other2->getTruthParticle()->isFromTp()) );
    }
    if( getTruthParticle()->isFromWm() || getTruthParticle()->isFromTm() ) {
      return ( (other1->getTruthParticle()->isFromWm() || other1->getTruthParticle()->isFromTm()) &&
	       (other2->getTruthParticle()->isFromWm() || other2->getTruthParticle()->isFromTm()) );
    }
    return kFALSE;
  }

  //
  // Function returns true if the input reco particles are from the W decay
  // of the same top decay as this particle.
  //
  Bool_t fromCommonWofSameTop( RecoParticle *other1 , RecoParticle *other2 ) {
    if( ! isTruthMatched() ) return kFALSE;
    if( ! getTruthParticle()->isFromT() ) return kFALSE;
    if( ! other1->isTruthMatched() ) return kFALSE;
    if( ! other2->isTruthMatched() ) return kFALSE;
    if( getTruthParticle()->isFromTp() )
      return ( other1->getTruthParticle()->isFromWp() && other2->getTruthParticle()->isFromWp() );
    if( getTruthParticle()->isFromTm() )
      return ( other1->getTruthParticle()->isFromWm() && other2->getTruthParticle()->isFromWm() );
    return kFALSE;
  }
  
  //
  // Function matches each reconstructed jet/lepton with the nearest truth
  // jet/lepton. Each particle is matched with one truth particle at most
  // by construction.
  //
  static void truthMatch( std::vector<RecoParticle*> reco , std::vector<TruthParticle*> truth ) {
    // Loop over the reco particles
    for( std::size_t ireco = 0 ; ireco < reco.size() ; ireco++ ) {
      RecoParticle *r1 = reco[ireco];

      Double_t bestdR = 999;
      int bestId      = -1;
      
      // For each reco particle, find the nearest truth particle of the same type
      for( std::size_t itruth = 0 ; itruth < truth.size() ; itruth++ ) {

	TruthParticle *t1 = truth[itruth];

	if( ! r1->matchesPDG(t1->getPDG()) ) continue;
	Double_t tmpdR1 = r1->deltaR(t1);
	if( tmpdR1 > 4.0 ) continue;
	if( tmpdR1 < bestdR ) {

	  Bool_t best_possible_match = kTRUE;
	  // We found a new nearest truth candidate
	  // But before counting these we need to make sure it's not an even better match for a different
	  // reconstructed particle
	  
	  for( std::size_t jreco = ireco+1 ; jreco < reco.size() ; jreco++ ) {
	    RecoParticle *r2 = reco[jreco];
	    if( ! r2->matchesPDG(t1->getPDG()) ) continue;
	    Double_t tmpdR2 = r2->deltaR(t1);
	    if( tmpdR2 < tmpdR1 ) {
	      best_possible_match = kFALSE;
	      break;
	    }
	  }

	  if( best_possible_match ) {
	    bestdR = tmpdR1;
	    bestId = int(itruth);
	  }
	  
	}
	  
      }

      // If a match was found, then update the RecoParticle information
      // and be sure to remove the truth particle from the vector so it doesn't get matched twice
      if( bestId >= 0 ) {
	r1->setTruthParticle( truth[bestId] );
	truth.erase( truth.begin() + bestId );
      }
      
    }
  }

};
  
#endif
