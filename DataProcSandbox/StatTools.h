#ifndef _STATTOOLS_H_
#define _STATTOOLS_H_

#include <iostream>
#include <vector>
#include <assert.h>

#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TFile.h>
#include <TMath.h>
#include <TIterator.h>
#include <TStopwatch.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooMinimizer.h>
#include <RooWorkspace.h>
#include <RooNLLVar.h>
#include <RooAbsArg.h>

#include "Report.h"

namespace
stats
{

  struct TwoSidedMeasurement {
    Double_t median , lo68 , hi68 , lo95 , hi95;
    Double_t dmedian , dlo68 , dhi68 , dlo95 , dhi95;
    Double_t mad; // median-absolute-deviation
    Double_t mean , dmean;
    Double_t xval; // in case this measurement is associated with some other parameter
    TString getString() {
      return TString::Format( "{ %g [ %g , %g , %g ] %g }" , lo95 , lo68 , median , hi68 , hi95 );
    }
    TString getMedianString() {
      return TString::Format( "%g ± %g" , median , mad );
    }
    TString getMeanString() {
      return TString::Format( "%g ± %g" , mean , dmean );
    }
    Double_t eylo() { return median - lo68; }
    Double_t eyhi() { return hi68 - median; }
    Double_t eylo95() { return median - lo95; }
    Double_t eyhi95() { return hi95 - median; }
  };

  struct OneSidedMeasurement {
    Double_t hi68 , hi95 , hi99;
    Double_t median , mad;
    Double_t mean , dmean;
    Double_t xval; // in case this measurement is associated with some other parameter
    TString getString() {
      return TString::Format( "50%%:%g, 68%%:%g, 95%%:%g, 99%%:%g" , median , hi68 , hi95 , hi99 );
    }
    TString getMedianString() {
      return TString::Format( "%g ± %g" , median , mad );
    }
    TString getMeanString() {
      return TString::Format( "%g ± %g" , mean , dmean );
    }
  };

  static Double_t getChi2Integral( Double_t x , Bool_t one_sided = kFALSE , Int_t ndf = 1 ) {
    if( x <= 0.0 ) return 0.0;
    if( ndf == 1 ) {
      if( one_sided ) return 0.5 + 0.5*TMath::Erf( TMath::Sqrt( 0.5 * x ) );
      return TMath::Erf( TMath::Sqrt( 0.5 * x ) );
    }
    assert( false );
  }

  static Double_t probToSigma( Double_t prob ) {
    assert (prob >=0 && prob <=1);
    prob = TMath::Min( prob , 0.9999999999999999 ); //else it returns 0 instead of something big and negative.
    prob = TMath::Max( prob , 0.0000000000000001 ); //else it returns 0 instead of something big.
    return sqrt(2.0) * TMath::ErfInverse(1-2*prob);
  }

  static Double_t getBinomialUncert( Double_t pval , Int_t ntrials ) {
    // only valid for large ntrials, better for pval close to 50%
    return sqrt( pval * (1.0-pval) / Double_t(ntrials) );
  }

  template< typename valT >
  void calculatePval( Double_t &pval , Double_t &dpval , valT successes , valT trials ) {
    pval = Double_t(successes) / Double_t(trials);
    dpval = getBinomialUncert( pval , Int_t(trials) );
  }

  template< typename valT >
  void calculatePval( Double_t &pval , Double_t &dpval , Double_t &zval , valT successes , valT trials ) {
    calculatePval( pval , dpval , successes , trials );
    zval = probToSigma( pval );
  }

  // calculate zval from pval
  template< typename valT > static valT probToSigma( valT prob ) {
    assert (prob >=0 && prob <=1);
    prob = TMath::Min( prob , 0.9999999999999999 ); //else it returns 0 instead of something big and negative.
    prob = TMath::Max( prob , 0.0000000000000001 ); //else it returns 0 instead of something big.
    return TMath::Sqrt(2.0) * TMath::ErfInverse(1-2*prob);
  }

  template< typename valT > static valT getMedian( std::vector<valT> vec ) {
    std::sort( vec.begin() , vec.end() );
    Int_t iter = Int_t( Double_t(vec.size()) * 0.5 );
    report::printassert( iter >= 0 && iter < vec.size() , "identified median value of vector at iter %i, but vector size = %i" , iter , Int_t(vec.size()) );
    return vec[iter];
  }

  template< typename valT > static valT getMean( std::vector<valT> vec ) {
    valT sum( 0 );
    for( std::size_t i = 0 ; i < vec.size() ; i++ ) sum += vec[i];
    return sum / valT(vec.size());
  }

  template< typename valT > static valT getRMS( std::vector<valT> vec , valT *rmsError = 0 ) {
    valT mean = getMean(vec);
    valT sum( 0 );
    valT N( vec.size() );
    for( valT v : vec ) sum += ( (mean - v)*(mean - v) );
    sum /= N;
    valT rmsVal = TMath::Sqrt( sum );
    if( rmsError ) (*rmsError) = rmsVal / TMath::Sqrt( valT(2) * N );
    return rmsVal;
  }

  // Calculate RMS ignoring outliers that are > 3*actual_RMS from mean
  template< typename valT > static valT getATLASRMS( std::vector<valT> vec , valT *rmsError = 0 , Int_t *nRej = 0 ) {
    return  getRMS( vec , rmsError );
    valT mean = getMean(vec);
    valT realRMS = getRMS(vec);
    valT sum( 0 );
    valT N( 0 );
    if( nRej ) (*nRej) = 0;
    for( valT v : vec ) {
      if( TMath::Abs(mean-v) > 3.0*realRMS ) {
	if( nRej ) (*nRej)++;
      } else {
	sum += ( (mean-v)*(mean-v) );
	N++;
      }
    }
    sum /= N;
    valT rmsVal = TMath::Sqrt( sum );
    if( rmsError ) (*rmsError) = rmsVal / TMath::Sqrt( valT(2) * N );
    return rmsVal;
  }
  
  template< typename valT > void calculateMedian( std::vector<valT> vec , valT &median , valT &mad ) {
    median = getMedian( vec );
    std::vector<valT> mad_vec;
    for( std::size_t i = 0 ; i < vec.size() ; i++ ) mad_vec.push_back( TMath::Abs( median - vec[i] ) );
    mad = getMedian( mad_vec );
  }

  template< typename valT > void calculateMean( std::vector<valT> vec , valT &mean , valT &dmean ) {
    mean = getMean( vec );
    dmean = 0;
    for( std::size_t i = 0 ; i < vec.size() ; i++ ) dmean += TMath::Power( vec[i] - mean , 2 );
    dmean = TMath::Sqrt( dmean / TMath::Power( valT(vec.size()) , 2 ) );
  }

  static TwoSidedMeasurement getTwoSidedMeasurement( std::vector<Double_t> vec ) {
    TwoSidedMeasurement m;
    std::sort( vec.begin() , vec.end() );
    Double_t vecsize = Double_t( vec.size() );
    //m.median = vec[ Int_t( vecsize * 0.5 ) ];
    calculateMedian( vec , m.median , m.mad );
    calculateMean( vec , m.mean , m.dmean );
    //m.median = getMean( vec );
    m.lo68 = vec[ Int_t( vecsize * 0.16 ) ];
    m.hi68 = vec[ Int_t( vecsize * 0.84 ) ];
    m.lo95 = vec[ Int_t( vecsize * 0.025 ) ];
    m.hi95 = vec[ Int_t( vecsize * 0.975 ) ];
    return m;
  }

  static TwoSidedMeasurement getTwoSidedMeasurement( RooDataSet *set , TString var ) {
    assert( set );
    std::vector<Double_t> vec;
    for( Int_t i = 0 ; i < set->numEntries() ; i++ ) vec.push_back( set->get(i)->getRealValue(var) );
    return getTwoSidedMeasurement( vec );
  }

  static TwoSidedMeasurement getTwoSidedEfficiency( Int_t total , Int_t passing ) {
    TwoSidedMeasurement m;
    m.mean   = Double_t(passing)/Double_t(total);
    m.median = Double_t(passing)/Double_t(total);
    m.lo68   = TEfficiency::ClopperPearson(total,passing,0.68,kFALSE);
    m.hi68   = TEfficiency::ClopperPearson(total,passing,0.68,kTRUE);
    m.lo95   = TEfficiency::ClopperPearson(total,passing,0.95,kFALSE);
    m.hi95   = TEfficiency::ClopperPearson(total,passing,0.95,kTRUE);
    return m;
  }

  static void addToGraph( TGraphAsymmErrors *g , TwoSidedMeasurement m , Double_t xval , Double_t xerr = 0.1 ) {
    Int_t ip = g->GetN();
    g->SetPoint( ip , xval , m.mean );
    g->SetPointError( ip , xerr , xerr , m.mean - m.lo68 , m.hi68 - m.mean );
  }

  static void addToGraphs( TGraphAsymmErrors *g68 , TGraphAsymmErrors *g95 , TwoSidedMeasurement m , Double_t xval , Double_t xerr = 0.1 ) {
    assert( g68->GetN() == g95->GetN() );
    Int_t ip = g68->GetN();
    g68->SetPoint( ip , xval , m.mean );
    g68->SetPointError( ip , xerr , xerr , m.mean - m.lo68 , m.hi68 - m.mean );
    g95->SetPoint( ip , xval , m.mean );
    g95->SetPointError( ip , xerr , xerr , m.mean - m.lo95 , m.hi95 - m.mean );
  }

  
  static OneSidedMeasurement getOneSidedMeasurement( std::vector<Double_t> vec ) {
    OneSidedMeasurement m;
    std::sort( vec.begin() , vec.end() );
    Double_t vecsize = Double_t( vec.size() );
    calculateMedian( vec , m.median , m.mad );
    calculateMean( vec , m.mean , m.dmean );
    m.hi68 = vec[ Int_t( vecsize * 0.68 ) ];
    m.hi95 = vec[ Int_t( vecsize * 0.95 ) ];
    m.hi99 = vec[ Int_t( vecsize * 0.99 ) ];
    return m;
  }

  static OneSidedMeasurement getOneSidedMeasurement( RooDataSet *set , TString var ) {
    assert( set );
    std::vector<Double_t> vec;
    for( Int_t i = 0 ; i < set->numEntries() ; i++ ) vec.push_back( set->get(i)->getRealValue(var) );
    return getOneSidedMeasurement( vec );
  }

  template< typename valT > void calculatePvalFromQ( Double_t &pval , Double_t &dpval , std::vector<valT> vec , valT q ) {
    Int_t trials( vec.size() );
    Int_t successes = 0;
    for( std::size_t i = 0 ; i < vec.size() ; i++ ) {
      if( vec[i] >= q ) successes++;
    }
    calculatePval( pval , dpval , successes , trials );
  }

  static TwoSidedMeasurement getQFromConfLevel( std::vector<Double_t> vec , Double_t conflevel ) {
    TwoSidedMeasurement m;
    if( vec.size()==0 ) return m;
    Double_t dconflevel = getBinomialUncert( conflevel , Int_t(vec.size()) );
    //std::cout << "DEBUG : calculated dconflevel = " << dconflevel << " for conflevel = " << conflevel << " and vector size = " << vec.size() << std::endl;
    std::sort( vec.begin() , vec.end() );
    Double_t vecsize = Double_t( vec.size() );
    m.median = vec[ Int_t( vecsize * conflevel ) ];
    m.lo68 = vec[ Int_t( vecsize * (conflevel - dconflevel) ) ];
    m.hi68 = vec[ Int_t( vecsize * (conflevel + dconflevel) ) ];
    m.hi95 = -666.0;
    m.lo95 = -666.0;
    return m;
  }

  static TwoSidedMeasurement getPvalFromQ( TwoSidedMeasurement q , Bool_t one_sided = kTRUE ) {
    // assume H1 is distributed like a Chi-square distribution with 1 DOF
    TwoSidedMeasurement m;
    m.median = 1.0 - getChi2Integral( q.median , one_sided );
    m.lo68 = 1.0 - getChi2Integral( q.hi68 , one_sided );
    m.hi68 = 1.0 - getChi2Integral( q.lo68 , one_sided );
    m.lo95 = 1.0 - getChi2Integral( q.hi95 , one_sided );
    m.hi95 = 1.0 - getChi2Integral( q.lo95 , one_sided );
    m.dmedian = m.dlo68 = m.dhi68 = m.dlo95 = m.dhi95 = 0.0;
    return m;
  }


  static TwoSidedMeasurement getPvalFromQ( std::vector<Double_t> vec , TwoSidedMeasurement q ) {
    Int_t trials( vec.size() );
    Int_t ss_median = 0;
    Int_t ss_lo68 = 0;
    Int_t ss_hi68 = 0;
    Int_t ss_lo95 = 0;
    Int_t ss_hi95 = 0;
    for( std::size_t i = 0 ; i < vec.size() ; i++ ) {
      if( vec[i] >= q.median ) ss_median++;
      if( vec[i] >= q.lo68 ) ss_hi68++;
      if( vec[i] >= q.hi68 ) ss_lo68++;
      if( vec[i] >= q.lo95 ) ss_hi95++;
      if( vec[i] >= q.hi95 ) ss_lo95++;
    }
    TwoSidedMeasurement m;
    calculatePval( m.median , m.dmedian , ss_median , trials );
    calculatePval( m.lo68 , m.dlo68 , ss_lo68 , trials );
    calculatePval( m.hi68 , m.dhi68 , ss_hi68 , trials );
    calculatePval( m.lo95 , m.dlo95 , ss_lo95 , trials );
    calculatePval( m.hi95 , m.dhi95 , ss_hi95 , trials );
    return m;
  }

  static TwoSidedMeasurement getPvalFromQ( RooDataSet *qset , TString qname , TwoSidedMeasurement q ) {
    // same as above function, but with RooDataSet
    std::vector<Double_t> qvec;
    for( Int_t iq = 0 ; iq < qset->numEntries() ; iq++ ) qvec.push_back( qset->get(iq)->getRealValue(qname) );
    return getPvalFromQ( qvec , q );
  }

  static TwoSidedMeasurement getPvalFromQ( std::vector<Double_t> vec , Double_t qobs ) {
    Int_t trials( vec.size() );
    Int_t ss_median = 0;
    for( std::size_t i = 0 ; i < vec.size() ; i++ ) {
      if( vec[i] >= qobs ) ss_median++;
    }
    TwoSidedMeasurement m;
    Double_t dpval;
    calculatePval( m.median , dpval , ss_median , trials );
    m.lo68 = -666.;
    m.hi68 = -666.;
    m.lo95 = -666.;
    m.hi95 = -666.;
    return m;
  }

  void setupMinimizer( RooMinimizer &m , RooWorkspace *ws ) {
    m.setPrintLevel( -1 );
    m.setStrategy( 2 );
  }

  void dumpFitResults( RooAbsReal *nll , RooWorkspace *ws , TStopwatch &tw , const TString &channel , const Bool_t &quiet ) {

    // if the workspace has a variable to keep track of fit times, then set it
    if( ws->var("cputime") ) ws->var("cputime")->setVal( tw.CpuTime() );

    // dump best-fit parameters
    if( !quiet ) {
      report::info( tw , "fit summary" );
      if( ws->set("fitparams_"+channel) ) {
	TIterator *it = ws->set("fitparams_"+channel)->createIterator();
	RooAbsArg *v = 0;
	while( (v = (RooAbsArg*)it->Next()) ) {
	  if( ! v->isConstant() ) report::spaceprint( v , "a" );
	}
	delete it;
	v = 0;
      } else {
	assert( false );
      }
      report::spaceprint( "NLL = %g" , nll->getVal() );
      //report::spaceprint( nll , "a" );
    }

  }

  static Double_t fastFitNoErrors( RooAbsReal *nll , RooWorkspace *ws , const TString &channel = "all" , const Int_t &printlevel = -1 , const Int_t &strategy = 2 ) { //const Bool_t &quiet = kFALSE ) {

    // just use migrad & don't worry about getting the errors right

    // perform fit and time
    RooMinimizer m( *nll );
    setupMinimizer( m , ws );
    m.setPrintLevel( printlevel < -1 ? -1 : printlevel );
    m.setStrategy( strategy );
    TStopwatch tw;
    tw.Start();
    m.migrad();
    tw.Stop();
    Double_t result = nll->getVal();

    dumpFitResults( nll , ws , tw , channel , printlevel < -1 );
    return result;
  }

  static Double_t fitWithErrors( RooAbsReal *nll , RooWorkspace *ws , const TString &channel = "all" , const Bool_t &quiet = kFALSE ) {

    // perform fit and time
    RooMinimizer m( *nll );
    setupMinimizer( m , ws );
    TStopwatch tw;
    tw.Start();
    m.migrad();
    m.simplex();
    m.migrad();
    m.improve();
    tw.Stop();
    Double_t result = nll->getVal();

    dumpFitResults( nll , ws , tw , channel , quiet );
    return result;
  }

  
};


#endif
