#ifndef _EXTRAPLOTTOOLS_H_
#define _EXTRAPLOTTOOLS_H_


#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TKey.h>
#include <TLegend.h>

#include "Report.h"
#include "HistTools.h"
#include "Plotter.h"
#include "InputConfig.h"
#include "TreeReader.h"
#include "PhysicsProcess.h"
#include "Config.h"


namespace
plot
{

  //
  // Structs for keeping track of histogram configuratinos
  // 
  
  struct HistConfig1D {
    TString xname;
    TString xtitle;
    Int_t xbins;
    Double_t xmin;
    Double_t xmax;
    Double_t xunits;
    HistConfig1D( TString _xname , TString _xtitle , Int_t _xbins , Double_t _xmin , Double_t _xmax , Double_t _xunits )
      : xname( _xname )
      , xtitle( _xtitle )
      , xbins( _xbins )
      , xmin( _xmin )
      , xmax( _xmax )
      , xunits( _xunits )
    {}
  };

  /*
  struct HistConfig2D {
    TString name;
    TString xtitle;
    Int_t nbins;
    Double_t xmin;
    Double_t xmax;
    Double_t units;
    HistConfig1D( TString _name , TString _xtitle , Int_t _nbins , Double_t _xmin , Double_t _xmax , Double_t _units )
      : name( _name )
      , xtitle( _xtitle )
      , nbins( _nbins )
      , xmin( _xmin )
      , xmax( _xmax )
      , units( _units )
    {}
  };
  */  

  // For ordering hists by integral
  static bool TH1DIntegralCompareLT( TH1D *lhs , TH1D *rhs ) { return lhs->Integral() < rhs->Integral(); }
  static bool TH1DIntegralCompareGT( TH1D* lhs , TH1D *rhs ) { return lhs->Integral() > rhs->Integral(); }


  //
  // 
  static std::map<TString,TH1D*> getTH1Dmap( const PhysicsProcess &proc , InputConfig *in , TreeReader *tr , std::vector<HistConfig1D> hconfigs , Bool_t mcweights ) {

    report::info( "Filling histograms for process %s" , proc.getName().Data() );
    
    std::map<TString,TH1D*> hmap;

    for( HistConfig1D hconfig : hconfigs ) {
      hmap[hconfig.xname] = new TH1D( hconfig.xname+proc.getName() , TString::Format(";%s;Events / Bin",hconfig.xtitle.Data()) , hconfig.xbins , hconfig.xmin , hconfig.xmax );
    }

    Double_t xsec_weight , event_weight , val;

    TH1D *h_total	     = new TH1D( "h_total_"+proc.getName() , "" , 1 , 0 , 1 );
    TH1D *h_passing	     = new TH1D( "h_passing_"+proc.getName() , "" , 10 , 20 , 30 );
    TH1D *h_passing_weighted = new TH1D( "h_passing_weighted_"+proc.getName() , "" , 10 , 20 , 30 );

    tr->setupOutputTree( proc.getName() );
    tr->setSignalMode( proc.isSignal() );
    tr->initCutflow();
    
    Long64_t nev_total = 0;
    Long64_t pev_total = 0;
    
    // Loop over the input ROOT files corresponding to this process
    for( size_t ifile = 0 ; ifile < proc.numSamples() ; ifile++ ) {
      
      // load Tree for this sample
      TTree *tree_tmp = proc.getSampleTree(ifile);
      tr->setTree( tree_tmp );

      h_total->SetBinContent( 1 , h_total->GetBinContent(1) + proc.getSampleNevents(ifile) );

      Double_t sumw_file = 0.;
      
      // Loop over events in tree
      Long64_t nev = tree_tmp->GetEntries();
      Long64_t pev = 0;
      report::startProgressBar( proc.getSampleTag(ifile).Length() < 20 ? proc.getSampleTag(ifile) : "" );
      for( Long64_t iev = 0 ; iev < nev ; iev++ ) {

	report::updateProgressBar( Double_t(iev+1) , Double_t(nev) , TString::Format(" SumW=%0.1f ε=%0.3f",sumw_file,Double_t(pev)/Double_t(iev)) );
	
	if( ! tr->getEntry(iev) ) continue;
	
	if( ! tr->passesSelection() ) continue;

	//tr->debugInfo();
	//assert( false );

	pev++;
	
	// xsec weight
	xsec_weight = proc.getSampleWeight(ifile);

	// extra per-event weights
	event_weight = 1.0;
	if( mcweights ) {
	  event_weight *= (*tr)["weight_mc"] * (*tr)["weight_leptonSF"] * (*tr)["weight_bTagSF_77"];
	  if( config::UsePileupWeights )
	    event_weight *= (*tr)["weight_pileup"];
	}

	sumw_file += event_weight * xsec_weight;
	
	Double_t leadingLeptonPt = tr->leadingLeptonPt() * tth::GeV;
	leadingLeptonPt = TMath::Min( 29.999 , leadingLeptonPt );
	leadingLeptonPt = TMath::Max( 20.001 , leadingLeptonPt );
	h_passing->Fill( leadingLeptonPt , event_weight );
	h_passing_weighted->Fill( leadingLeptonPt , event_weight * xsec_weight );

	tr->fillOutputTree( ifile , iev , xsec_weight * event_weight );

	for( HistConfig1D hconfig : hconfigs ) {
	  if( !mcweights && hconfig.xname.Contains("weight_") ) continue;
	  tr->fillHist( hmap[hconfig.xname] , hconfig.xname , hconfig.xunits , hconfig.xname.Contains("weight_") ? xsec_weight : xsec_weight * event_weight );
	}
	
      } // end loop over tree entries
      delete tree_tmp;

      nev_total += nev;
      pev_total += pev;
      
      //report::updateProgressBar( nev , nev );
      //report::debug( "Integral = %g" , sumw_file );
      
    } // end loop over files

    // Fix the x-labels for inclusive histograms
    for( HistConfig1D hconfig : hconfigs ) {
      if( hconfig.xname.EndsWith("_INCL") ) {
	for( Int_t ibin = 1 ; ibin <= hmap[hconfig.xname]->GetNbinsX() ; ibin++ ) {
	  hmap[hconfig.xname]->GetXaxis()->SetBinLabel( ibin , TString::Format("#geq%g",hmap[hconfig.xname]->GetBinLowEdge(ibin)) );
	}
	hmap[hconfig.xname]->GetXaxis()->CenterLabels();
      }
    }

    hmap["passing"]	     = h_passing;
    hmap["passing_weighted"] = h_passing_weighted;
    hmap["total"]	     = h_total;

    report::info( "Total Efficiency = %li / %li = %g" , pev_total , nev_total , Double_t(pev_total) / Double_t(nev_total) );

    tr->printCutflow();
    
    return hmap;
  }

  
  //
  // Fill histogram corresponding to PhysicsProcss "proc"
  //
  static TH1D* getTH1D( const PhysicsProcess &proc , InputConfig *in , TreeReader *tr , HistConfig1D hconfig , Bool_t mcweights ) {
    std::vector<HistConfig1D> hconfigVec;
    hconfigVec.push_back( hconfig );
    return getTH1Dmap(proc,in,tr,hconfigVec,mcweights)[hconfig.xname];
  }
  
  
};

#endif
