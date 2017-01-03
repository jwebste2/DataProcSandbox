#ifndef _RECOPARTICLECOLLECTION_H_
#define _RECOPARTICLECOLLECTION_H_

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
#include "TF1.h"
#include "TRandom3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "Report.h"
#include "Particle.h"




class
RecoParticleCollection
{

  //
  // Class contains a collection of reco particles, can be used to calculate
  // shape information for the collection of particles
  //

private:

  std::vector<RecoParticle*> v_particles;

  TLorentzVector v4;
  Int_t wpsum;
  Double_t ptsum;

  std::vector<Double_t> v_eigenvals;
  std::vector<Double_t> v_eigenvalsT;
  std::vector<Double_t> v_eigenvalsO;

public:

  
  RecoParticleCollection() {}
  RecoParticleCollection( std::vector<RecoParticle*> v )
    : v_particles( v )
  {
    assert( v_particles.size() > 0 );

    // Build a combined 4-vector for the full system
    wpsum = 0;
    ptsum = 0;
    v4	  = TLorentzVector();
    for( std::size_t ip = 0 ; ip < v_particles.size() ; ip++ ) {
      v4    += v_particles[ip]->getV4();
      ptsum += v_particles[ip]->getPt();
      if( v_particles[ip]->getType() == RecoParticle::JET || v_particles[ip]->getType() == RecoParticle::BJET )
	wpsum += v_particles[ip]->getTagLevel();
    }

    // Build momentum tensors to calculate eigenvalues, which will
    // provide shape information
    TMatrixDSym momentumTensor( 3 );
    TMatrixDSym momentumTensorT( 2 );
    TMatrixDSym momentumTensorO( 3 );

    Double_t Sxx=0, Sxy=0, Sxz=0, Syy=0, Syz=0, Szz=0, normal=0;
    Double_t Oxx=0, Oxy=0, Oxz=0, Oyy=0, Oyz=0, Ozz=0, Onormal=0;

    for( RecoParticle *p : v_particles ) {
      Sxx += p->getV4().Px() * p->getV4().Px();
      Sxy += p->getV4().Px() * p->getV4().Py();
      Sxz += p->getV4().Px() * p->getV4().Pz();
      Syy += p->getV4().Py() * p->getV4().Py();
      Syz += p->getV4().Py() * p->getV4().Pz();
      Szz += p->getV4().Pz() * p->getV4().Pz();
      normal += p->getV4().P() * p->getV4().P();

      Oxx += (p->getV4().Px() * p->getV4().Px())/p->getV4().P();
      Oxy += (p->getV4().Px() * p->getV4().Py())/p->getV4().P();
      Oxz += (p->getV4().Px() * p->getV4().Pz())/p->getV4().P();
      Oyy += (p->getV4().Py() * p->getV4().Py())/p->getV4().P();
      Oyz += (p->getV4().Py() * p->getV4().Pz())/p->getV4().P();
      Ozz += (p->getV4().Pz() * p->getV4().Pz())/p->getV4().P();
      Onormal += p->getV4().P();
    }

    momentumTensor( 0 , 0 ) = Sxx/normal;
    momentumTensor( 0 , 1 ) = Sxy/normal;
    momentumTensor( 0 , 2 ) = Sxz/normal;
    momentumTensor( 1 , 0 ) = momentumTensor( 0 , 1 );
    momentumTensor( 1 , 1 ) = Syy/normal;
    momentumTensor( 1 , 2 ) = Syz/normal;
    momentumTensor( 2 , 0 ) = momentumTensor( 0 , 2 );
    momentumTensor( 2 , 1 ) = momentumTensor( 1 , 2 );
    momentumTensor( 2 , 2 ) = Szz/normal;

    momentumTensorT( 0 , 0 ) = momentumTensor( 0 , 0 );
    momentumTensorT( 0 , 1 ) = momentumTensor( 0 , 1 );
    momentumTensorT( 1 , 1 ) = momentumTensor( 1 , 1 );
    momentumTensorT( 1 , 0 ) = momentumTensor( 1 , 0 );

    momentumTensorO( 0 , 0 ) = Oxx/Onormal;
    momentumTensorO( 0 , 1 ) = Oxy/Onormal;
    momentumTensorO( 0 , 2 ) = Oxz/Onormal;
    momentumTensorO( 1 , 0 ) = momentumTensorO( 0 , 1 );
    momentumTensorO( 1 , 1 ) = Oyy/Onormal;
    momentumTensorO( 1 , 2 ) = Oyz/Onormal;
    momentumTensorO( 2 , 0 ) = momentumTensorO( 0 , 2 );
    momentumTensorO( 2 , 1 ) = momentumTensorO( 1 , 2 );
    momentumTensorO( 2 , 2 ) = Ozz/Onormal;

    // Compute eigenvalues

    TMatrixDSymEigen Eigenvalues( momentumTensor );
    TVectorD eigenVec = Eigenvalues.GetEigenValues();
    for( std::size_t i = 0 ; i < 3 ; i++ )
      v_eigenvals.push_back( eigenVec(i) );
    
    TMatrixDSymEigen EigenvaluesT(momentumTensorT);
    TVectorD eigenVecT = EigenvaluesT.GetEigenValues();
    for( std::size_t i = 0 ; i < 2 ; i++ )
      v_eigenvalsT.push_back( eigenVecT(i) );

    TMatrixDSymEigen EigenvaluesO(momentumTensorO);
    TVectorD eigenVecO = EigenvaluesO.GetEigenValues();
    for( std::size_t i = 0 ; i < 3 ; i++ )
      v_eigenvalsO.push_back( eigenVecO(i) );

  }

  ~RecoParticleCollection() {}

  TLorentzVector getV4() { return v4; }
  Double_t getPt() { return v4.Pt(); }
  Double_t getEta() { return v4.Eta(); }
  Double_t getPhi() { return v4.Phi(); }
  Double_t getM() { return v4.M(); }
  Double_t getMt() { return TMath::Sqrt( v4.Et()*v4.Et() - v4.Px()*v4.Px() - v4.Py()*v4.Py() ); }
  Int_t getWpsum() { return wpsum; }
  Double_t getPtsum() { return ptsum; }

  Double_t getAplanarity() { return 1.5 * v_eigenvals[2]; }
  Double_t getAplanarityO() { return 1.5 * v_eigenvalsO[2]; }
  Double_t getSphericity() { return 1.5 * ( v_eigenvals[1] + v_eigenvals[2] ); }
  Double_t getSphericityO() { return 1.5 * ( v_eigenvalsO[1] + v_eigenvalsO[2] ); }
  Double_t getSphericityT() { return 2.0 * v_eigenvalsT[1] / ( v_eigenvalsT[0] + v_eigenvalsT[1] ); }
  Double_t getPlanarity() { return v_eigenvals[1] - v_eigenvals[2]; }
  Double_t getVariableC() {
    return ( 3.0 * (v_eigenvals[0]*v_eigenvals[1]
		    + v_eigenvals[0]*v_eigenvals[2]
		    + v_eigenvals[1]*v_eigenvals[2]) );
  }
  Double_t getVariableD() { return 27.0 * v_eigenvals[0] * v_eigenvals[1] * v_eigenvals[2]; }
  Double_t getCircularity() {
    if( (v_eigenvals[0] + v_eigenvals[1]) == 0 ) return 0.0;
    return 2.0 * v_eigenvals[1] / ( v_eigenvals[0] + v_eigenvals[1] );
  }
  Double_t getPlanarFlow() {
    if( (v_eigenvals[0] + v_eigenvals[1]) * (v_eigenvals[0] + v_eigenvals[1]) == 0 ) return 0.0;
    return 4.0 * (v_eigenvals[0] * v_eigenvals[1]) /
      ((v_eigenvals[0] + v_eigenvals[1]) * (v_eigenvals[0] + v_eigenvals[1]));
  }

  Double_t getShape( TString n ) {
    if( n==TString("Aplanarity") )  return getAplanarity();
    if( n==TString("AplanarityO") ) return getAplanarityO();
    if( n==TString("Sphericity") )  return getSphericity();
    if( n==TString("SphericityO") ) return getSphericityO();
    if( n==TString("SphericityT") ) return getSphericityT();
    if( n==TString("Planarity") )   return getPlanarity();
    if( n==TString("VariableC") )   return getVariableC();
    if( n==TString("VariableD") )   return getVariableD();
    if( n==TString("Circularity") ) return getCircularity();
    if( n==TString("PlanarFlow") )  return getPlanarFlow();
    assert( false );
  }

  
};


#endif

