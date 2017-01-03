#ifndef _STACKPLOT_H_
#define _STACKPLOT_H_


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
#include "ExtraPlotTools.h"


class
StackPlotter
{

private:

  InputConfig *in;
  Plotter *p;
  TreeReader *tr;

  std::vector<plot::HistConfig1D> queue1D;
  //std::vetcor<plot::HistConfig2D> queue2D;
 
public:

  StackPlotter( InputConfig *_in , Plotter *_p , TreeReader *_tr ) {
    in = _in;
    p  = _p;
    tr = _tr;
  }
  
  ~StackPlotter() {}

  void addToQueue( TString xname , TString xtitle , Int_t xbins , Double_t xmin , Double_t xmax , Double_t xunits = 1.0 ) { queue1D.push_back( plot::HistConfig1D(xname,xtitle,xbins,xmin,xmax,xunits) ); }

  void drawDataMC( const Bool_t &use_event_weights = kTRUE ) {

    //
    // get filled histograms
    //
    std::vector< std::map<TString,TH1D*> > h_data_all = getFilledDataHistograms1D();
    std::vector< std::map<TString,TH1D*> > h_bkg_all  = getFilledBackgroundHistograms1D( use_event_weights );
    std::vector< std::map<TString,TH1D*> > h_sig_all  = getFilledSignalHistograms1D( use_event_weights );

    //
    // Draw all the plots
    //
    for( plot::HistConfig1D hconfig : queue1D ) {

      std::vector<TH1D*> h_data_vec;
      std::vector<TH1D*> h_bkg_vec;
      std::vector<TH1D*> h_sig_vec;

      for( std::size_t i = 0 ; i < in->numData() ; ++i ) {
	h_data_vec.push_back( h_data_all[i][hconfig.xname] );
	h_data_vec.back()->SetName( in->getData(i).getRootTitle() );
      }

      for( std::size_t i = 0 ; i < in->numBackgrounds() ; i++ ) {
	h_bkg_vec.push_back( h_bkg_all[i][hconfig.xname] );
	h_bkg_vec.back()->SetName( in->getBackground(i).getRootTitle() );
      }

      for( std::size_t i = 0 ; i < in->numSignals() ; ++i ) {
	h_sig_vec.push_back( h_sig_all[i][hconfig.xname] );
	h_sig_vec.back()->SetName( in->getSignal(i).getRootTitle() );
      }

      genericDraw( h_data_vec , h_sig_vec , h_bkg_vec , hconfig );
      
    }

    
  }


  void drawData() {

    //
    // get filled histograms
    //
    std::vector< std::map<TString,TH1D*> > h_data_all = getFilledDataHistograms1D();

    //
    // Draw all the plots
    //
    for( plot::HistConfig1D hconfig : queue1D ) {

      std::vector<TH1D*> h_data_vec;
      std::vector<TH1D*> h_bkg_vec;
      std::vector<TH1D*> h_sig_vec;

      for( std::size_t i = 0 ; i < in->numData() ; ++i ) {
	h_data_vec.push_back( h_data_all[i][hconfig.xname] );
	h_data_vec.back()->SetName( in->getData(i).getRootTitle() );
      }

      genericDraw( h_data_vec , h_sig_vec , h_bkg_vec , hconfig );
      
    }

    
  }


  void drawMC( const Bool_t &use_event_weights = kTRUE ) {

    //
    // get filled histograms
    //
    std::vector< std::map<TString,TH1D*> > h_bkg_all  = getFilledBackgroundHistograms1D( use_event_weights );
    std::vector< std::map<TString,TH1D*> > h_sig_all  = getFilledSignalHistograms1D( use_event_weights );

    //
    // Draw all the plots
    //
    for( plot::HistConfig1D hconfig : queue1D ) {

      std::vector<TH1D*> h_data_vec;
      std::vector<TH1D*> h_bkg_vec;
      std::vector<TH1D*> h_sig_vec;

      for( std::size_t i = 0 ; i < in->numBackgrounds() ; i++ ) {
	h_bkg_vec.push_back( h_bkg_all[i][hconfig.xname] );
	h_bkg_vec.back()->SetName( in->getBackground(i).getRootTitle() );
      }

      for( std::size_t i = 0 ; i < in->numSignals() ; ++i ) {
	h_sig_vec.push_back( h_sig_all[i][hconfig.xname] );
	h_sig_vec.back()->SetName( in->getSignal(i).getRootTitle() );
      }

      genericDraw( h_data_vec , h_sig_vec , h_bkg_vec , hconfig );
      
    }
    
  }


  void overlayMC( const Bool_t &use_event_weights = kTRUE ) {

    //
    // get filled histograms
    //
    std::vector< std::map<TString,TH1D*> > h_bkg_all  = getFilledBackgroundHistograms1D( use_event_weights );
    std::vector< std::map<TString,TH1D*> > h_sig_all  = getFilledSignalHistograms1D( use_event_weights );

    //
    // Draw all the plots
    //
    for( plot::HistConfig1D hconfig : queue1D ) {

      std::vector<TH1D*> h_bkg_vec;
      std::vector<TH1D*> h_sig_vec;

      for( std::size_t i = 0 ; i < in->numBackgrounds() ; i++ ) {
	h_bkg_vec.push_back( h_bkg_all[i][hconfig.xname] );
	h_bkg_vec.back()->SetName( in->getBackground(i).getRootTitle() );
      }

      for( std::size_t i = 0 ; i < in->numSignals() ; ++i ) {
	h_sig_vec.push_back( h_sig_all[i][hconfig.xname] );
	h_sig_vec.back()->SetName( in->getSignal(i).getRootTitle() );
      }

      genericOverlay( h_sig_vec , h_bkg_vec , hconfig );
      
    }
    
  }
  


private:

  std::vector< std::map<TString,TH1D*> > getFilledDataHistograms1D() {
    std::vector< std::map<TString,TH1D*> > res;
    for( std::size_t i = 0 ; i < in->numData() ; ++i )
      res.push_back( plot::getTH1Dmap( in->getData(i) , in , tr , queue1D , kFALSE ) );
    return res;
  }


  
  std::vector< std::map<TString,TH1D*> > getFilledSignalHistograms1D( const Bool_t &use_event_weights ) {

    std::vector< std::map<TString,TH1D*> > res;

    // Also write some extra histograms to output ROOT file with event counts w.r.t pT
    TH1D *h_total	     = (TH1D*)0;
    TH1D *h_passing	     = (TH1D*)0;
    TH1D *h_passing_weighted = (TH1D*)0;

    for( std::size_t i = 0 ; i < in->numSignals() ; ++i ) {
      res.push_back( plot::getTH1Dmap( in->getSignal(i) , in , tr , queue1D , use_event_weights ) );
      if( res.back()["passing"] ) {
	htools::safe_add( h_passing , res.back()["passing"] , "h_passing_sig" );
	htools::safe_add( h_passing_weighted , res.back()["passing_weighted"] , "h_passing_weighted_sig" );
	htools::safe_add( h_total , res.back()["total"] , "h_total_sig" );
	delete res.back()["passing"];
	delete res.back()["passing_weighted"];
	delete res.back()["total"];
      }
    }

    if( h_passing ) {
      p->write( h_passing , "h_passing_sig" );
      p->write( h_passing_weighted , "h_passing_weighted_sig" );
      p->write( h_total , "h_total_sig" );
      delete h_passing;
      delete h_passing_weighted;
      delete h_total;
    }

    return res;

  }


  
  std::vector< std::map<TString,TH1D*> > getFilledBackgroundHistograms1D( const Bool_t &use_event_weights ) {

    std::vector< std::map<TString,TH1D*> > res;

    // Also write some extra histograms to output ROOT file with event counts w.r.t. pT
    TH1D *h_total	     = (TH1D*)0;
    TH1D *h_passing	     = (TH1D*)0;
    TH1D *h_passing_weighted = (TH1D*)0;
    
    for( std::size_t i = 0 ; i < in->numBackgrounds() ; ++i ) {
      res.push_back( plot::getTH1Dmap( in->getBackground(i) , in , tr , queue1D , use_event_weights ) );
      if( res.back()["passing"] ) {
	htools::safe_add( h_passing , res.back()["passing"] , "h_passing_bkg" );
	htools::safe_add( h_passing_weighted , res.back()["passing_weighted"] , "h_passing_weighted_bkg" );
	htools::safe_add( h_total , res.back()["total"] , "h_total_bkg" );
	p->write( res.back()["passing"] , "h_passing_bkg_"+in->getBackground(i).getName() );
	p->write( res.back()["passing_weighted"] , "h_passing_weighted_bkg_"+in->getBackground(i).getName() );
	p->write( res.back()["total"] , "h_total_bkg_"+in->getBackground(i).getName() );
	delete res.back()["passing"];
	delete res.back()["passing_weighted"];
	delete res.back()["total"];
      }
    }

    if( h_passing ) {
      p->write( h_passing , "h_passing_bkg" );
      p->write( h_passing_weighted , "h_passing_weighted_bkg" );
      p->write( h_total , "h_total_bkg" );
      delete h_passing;
      delete h_passing_weighted;
      delete h_total;
    }

    return res;
  }


  // Draw function for vectors of 1D histograms representing data/signal/backgrounds
  void genericDraw( std::vector<TH1D*> &h_data_vec , std::vector<TH1D*> &h_sig_vec , std::vector<TH1D*> &h_bkg_vec , plot::HistConfig1D hconfig ) {

    report::info( "drawing stack plot of %s" , hconfig.xname.Data() );

    // Stylize data histograms
    for( std::size_t i = 0 ; i < h_data_vec.size() ; ++i ) {
      htools::setup( h_data_vec[i] , in->getData(i).getColor() );
      h_data_vec[i]->SetMarkerStyle( 20 );
    }

    // Stylize background histograms
    for( std::size_t i = 0 ; i < h_bkg_vec.size() ; ++i )
      htools::setup( h_bkg_vec[i] , in->getBackground(i).getColor() , 1 , 1001 );
    if( in->getSortBackgrounds() ) 
      std::sort( h_bkg_vec.begin() , h_bkg_vec.end() , plot::TH1DIntegralCompareLT );

    // Stylize signal histograms
    for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i )
      htools::setup( h_sig_vec[i] , in->getSignal(i).getColor() );

    TLegend *legend = new TLegend();
    TLegend *yields = new TLegend();
    TPaveText *text = new TPaveText();

    text->InsertText( htools::AtlasDefault );
    text->InsertText( tr->getDressedSelectionTitle() );
    text->InsertText( "" );
    text->InsertText( TString::Format("#intLdt = %0.2f fb^{-1}, #sqrt{s} = 13 TeV",in->getLuminosity()/1000.0) );
    text->InsertText( "" );

    //
    // Fill legends and stylize
    //

    Double_t ymax ( 0 );
    
    for( std::size_t i = 0 ; i < h_data_vec.size() ; ++i ) {
      legend->AddEntry( h_data_vec[i] , h_data_vec[i]->GetName() , "lpe" );
      yields->AddEntry( (TH1D*)NULL , "Blinded" /*TString::Format("%0.1f",h_data->Integral())*/ , "" );
      ymax = TMath::Max( ymax , h_data_vec[i]->GetMaximum() );
    }

    for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) {
      legend->AddEntry( h_sig_vec[i] , h_sig_vec[i]->GetName() , "l" );
      yields->AddEntry( (TH1D*)NULL , TString::Format("%0.1f",h_sig_vec[i]->Integral()) , "" );
      ymax = TMath::Max( ymax , h_sig_vec[i]->GetMaximum() );
    }
    
    Double_t sum_bkg ( 0 );
    for( Int_t i = h_bkg_vec.size()-1 ; i >= 0 ; i-- ) {
      legend->AddEntry( h_bkg_vec[i] , h_bkg_vec[i]->GetName() , "f" );
      yields->AddEntry( (TH1D*)NULL , TString::Format("%0.1f",h_bkg_vec[i]->Integral()) , "" );
      sum_bkg += h_bkg_vec[i]->Integral();
      for( Int_t j = 0 ; j < i ; j++ ) h_bkg_vec[i]->Add( h_bkg_vec[j] );
      ymax = TMath::Max( ymax , h_bkg_vec[i]->GetMaximum() );
    }

    if( h_bkg_vec.size() > 0 ) {
      legend->AddEntry( (TH1D*)NULL , "sum(BKG)" , "" );
      yields->AddEntry( (TH1D*)NULL , TString::Format("%0.1f",sum_bkg) , "" );
    }
    
    htools::setup( legend , htools::LEFT , htools::TOP , 0.5 );
    htools::setup( yields , htools::CENTER , htools::TOP , 0.5 );
    htools::setup( text , htools::RIGHT , htools::TOP , 0.5 );

    // Use consistent y-max for all histograms
    ymax *= 1.8;
    for( std::size_t i = 0 ; i < h_data_vec.size() ; ++i ) h_data_vec[i]->SetMaximum( ymax );
    for( std::size_t i = 0 ; i < h_bkg_vec.size() ; ++i ) h_bkg_vec[i]->SetMaximum( ymax );
    for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) h_sig_vec[i]->SetMaximum( ymax );

    // Possibly use log scaling
    Bool_t logy = kFALSE;
    if( hconfig.xname.Contains("eta") || hconfig.xname.Contains("phi") || hconfig.xname==TString("mu") ) {
      // keep linear scale
    } else {
      logy = kTRUE;
      for( std::size_t i = 0 ; i < h_data_vec.size() ; ++i ) { h_data_vec[i]->SetMinimum(0.1); h_data_vec[i]->SetMaximum( 100 * ymax ); }
      for( std::size_t i = 0 ; i < h_bkg_vec.size() ; ++i ) { h_bkg_vec[i]->SetMinimum(0.1); h_bkg_vec[i]->SetMaximum( 100 * ymax ); }
      for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) { h_sig_vec[i]->SetMinimum(0.1); h_sig_vec[i]->SetMaximum( 100 * ymax ); }
    }
    
    
    
    //
    // Create histograms representing sum(bkg) and ratio w.r.t. data
    //
    
    TH1D *h_uncert     = (TH1D*)0;
    TH1D *h_nouncert   = (TH1D*)0;
    TH1D *h_ratio_bkg  = (TH1D*)0;
    TH1D *h_ratio_data = (TH1D*)0;
    TLine *line	       = (TLine*)0;

    if( h_bkg_vec.size() > 0 ) {

      h_uncert = (TH1D*) h_bkg_vec.back()->Clone("h_uncert");
      h_uncert->SetFillColor( kBlack );
      h_uncert->SetLineColor( kBlack );
      h_uncert->SetFillStyle( 3245 );
      
      h_nouncert = (TH1D*) h_bkg_vec.back()->Clone("h_nouncert");
      for( Int_t ibin = 1 ; ibin <= h_nouncert->GetNbinsX() ; ibin++ ) h_nouncert->SetBinError(ibin,0.0);
    

      // Only do ratio plots if there is a single data histogram
      if( h_data_vec.size() == 1 ) {

	// calculate chi-square where norm of MC is allowed to float
	Double_t chi2ndf = h_data_vec[0]->Chi2Test( h_uncert , "UW CHI2/NDF" );
	text->InsertText( TString::Format("#chi^{2}/ndf = %0.2f",chi2ndf) );

	h_ratio_bkg = (TH1D*) h_uncert->Clone("h_ratio_bkg");
	h_ratio_bkg->Divide( h_nouncert );
	h_ratio_bkg->SetFillColor( kYellow - 7 );
	h_ratio_bkg->SetLineColor( kBlack );
	h_ratio_bkg->SetFillStyle( 1001 );
	h_ratio_bkg->SetMinimum( 0.5 );
	h_ratio_bkg->SetMaximum( 1.5 );

	h_ratio_data = (TH1D*) h_data_vec[0]->Clone("h_ratio_data");
	h_ratio_data->Divide( h_nouncert );
	
	line = new TLine( hconfig.xmin , 1 , hconfig.xmax , 1 );
	line->SetLineWidth( 2 );
	line->SetLineStyle( kDashed );
	line->SetLineColor( kRed );

	h_uncert->GetXaxis()->SetTitle( "" );
	h_uncert->GetXaxis()->SetLabelSize( 0 );
	h_uncert->GetXaxis()->SetLabelOffset( 999 );
	
	h_bkg_vec.back()->GetXaxis()->SetTitle( "" );
	h_bkg_vec.back()->GetXaxis()->SetLabelSize( 0 );
	h_bkg_vec.back()->GetXaxis()->SetLabelOffset( 999 );

	h_ratio_bkg->GetYaxis()->SetNdivisions( 5 , 5 , 0 , kTRUE );
	h_ratio_bkg->GetYaxis()->SetTitle( "Data / Pred." );


	//
	// Blind bins with non-negligable S/B
	//
	
	if( h_sig_vec.size()==1 ) {
	  for( Int_t ibin = 1 ; ibin <= h_sig_vec[0]->GetNbinsX() ; ibin++ ) {
	    if( h_sig_vec[0]->GetBinContent(ibin) / h_uncert->GetBinContent(ibin) > 0.02 ) {
	      report::warn( "blinding bin %i in distribution %s" , ibin , hconfig.xname.Data() );
	      h_data_vec[0]->SetBinContent( ibin , -1.0 );
	      h_data_vec[0]->SetBinError( ibin , 999999 );
	      h_ratio_data->SetBinContent( ibin , -1.0 );
	      h_ratio_data->SetBinError( ibin , 999999 );
	    }
	  }
	} // nSig == 1
	
	// Draw the plots
	p->getCanvas()->Clear();
	TPad *padtop = new TPad( "padtop" , "" , 0 , 0 , 1 , 1 );
	TPad *padbottom = new TPad( "padbottom" , "" , 0 , 0 , 1 , 1 );
	htools::setupTop( padtop );
	htools::setupBottom( padbottom );
	padbottom->SetGridy();
	if( logy ) padtop->SetLogy();

	padtop->Draw();
	padtop->cd();
	h_uncert->Draw( "e2" );
	for( Int_t i = h_bkg_vec.size()-1 ; i >= 0 ; i-- ) h_bkg_vec[i]->Draw( "hist same" );
	h_uncert->Draw( "e2 same" );
	for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) h_sig_vec[i]->Draw( "hist e same" );
	for( std::size_t i = 0 ; i < h_data_vec.size() ; ++i ) h_data_vec[i]->Draw( "pe0 same" );
	legend->Draw();
	yields->Draw();
	text->Draw();

	p->getCanvas()->cd();
	padbottom->Draw();
	padbottom->cd();
	h_ratio_bkg->Draw( "e2" );
	line->Draw();
	h_ratio_data->Draw( "pe0 same" );
	p->print();

	p->getCanvas()->cd();
	p->getCanvas()->Clear( "D" );

	delete padtop;
	delete padbottom;
	
      } // nData == 1
      else {

	p->getCanvas()->SetLogy( logy );

	h_uncert->Draw( "e2" );
	for( Int_t i = h_bkg_vec.size()-1 ; i >= 0 ; i-- ) h_bkg_vec[i]->Draw( "hist same" );
	h_uncert->Draw( "e2 same" );
	for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) h_sig_vec[i]->Draw( "hist e same" );
	legend->Draw();
	yields->Draw();
	text->Draw();
	p->print();
	
      }
      
    } // nBkg > 0
    else {

      // just plot signal and/or data
      p->getCanvas()->SetLogy( logy );

      if( h_sig_vec.size() > 0 ) h_sig_vec[0]->Draw( "hist" );
      else if( h_data_vec.size() > 0 ) h_data_vec[0]->Draw( "pe0" );
      else assert( false ); // ==> No histograms found!!
      
      for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) h_sig_vec[i]->Draw( "hist same" );
      for( std::size_t i = 0 ; i < h_data_vec.size() ; ++i ) h_data_vec[i]->Draw( "pe0 same" );
      legend->Draw();
      yields->Draw();
      text->Draw();
      p->print();
      
    }

    //
    // Write histograms to output ROOT file
    //
    
    if( h_data_vec.size() > 0 ) {
      p->write( h_data_vec[0] , TString::Format("h_data_%s",hconfig.xname.Data()).Data() );
      for( std::size_t i = 0 ; i < h_data_vec.size() ; ++i ) p->write( h_data_vec[i] , TString::Format("h_data_comp%i_%s",Int_t(i),hconfig.xname.Data()).Data() );
    }

    if( h_sig_vec.size() > 0 ) {
      p->write( h_sig_vec[0] , TString::Format("h_sig_%s",hconfig.xname.Data()).Data() );
      for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) p->write( h_sig_vec[i] , TString::Format("h_sig_comp%i_%s",Int_t(i),hconfig.xname.Data()).Data() );
    }

    if( h_bkg_vec.size() > 0 ) {
      p->write( h_uncert , TString::Format("h_mc_%s",hconfig.xname.Data()).Data() );
      for( std::size_t i = 0 ; i < h_bkg_vec.size() ; i++ ) p->write( h_bkg_vec[i] , TString::Format("h_bkg_comp%i_%s",Int_t(i),hconfig.xname.Data()) );
    }

    htools::safe_delete( h_bkg_vec );
    htools::safe_delete( h_sig_vec );
    htools::safe_delete( h_data_vec );
    htools::safe_delete( h_uncert );
    htools::safe_delete( h_nouncert );
    htools::safe_delete( h_ratio_bkg );
    htools::safe_delete( h_ratio_data );
    delete line;
    delete legend;
    delete yields;
    delete text;

  }



  // Overlay function for vectors of 1D histograms representing data/signal/backgrounds
  void genericOverlay( std::vector<TH1D*> &h_sig_vec , std::vector<TH1D*> &h_bkg_vec , plot::HistConfig1D hconfig ) {

    report::info( "overlaying plots of %s" , hconfig.xname.Data() );

    // Stylize background histograms
    for( std::size_t i = 0 ; i < h_bkg_vec.size() ; ++i )
      htools::setup( h_bkg_vec[i] , in->getBackground(i).getColor() );
    if( in->getSortBackgrounds() ) 
      std::sort( h_bkg_vec.begin() , h_bkg_vec.end() , plot::TH1DIntegralCompareGT );

    // Stylize signal histograms
    for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i )
      htools::setup( h_sig_vec[i] , in->getSignal(i).getColor() );

    TLegend *legend = new TLegend();
    TLegend *yields = new TLegend();
    TPaveText *text = new TPaveText();

    text->InsertText( htools::AtlasDefault );
    text->InsertText( tr->getDressedSelectionTitle() );
    text->InsertText( "" );
    text->InsertText( TString::Format("#intLdt = %0.1f fb^{-1}, #sqrt{s} = 14 TeV",in->getLuminosity()/1000.0) );
    text->InsertText( "" );

   
    //
    // Fill legends, stylize, and clone histograms
    //

    std::vector<TH1D*> h_vec;
    std::vector<TH1D*> h_vec_norm;
    
    Double_t sum_bkg ( 0 );
    for( std::size_t i = 0 ; i < h_bkg_vec.size() ; ++i ) {
      legend->AddEntry( h_bkg_vec[i] , h_bkg_vec[i]->GetName() , "l" );
      yields->AddEntry( (TH1D*)NULL , TString::Format("%0.1f",h_bkg_vec[i]->Integral()) , "" );
      sum_bkg += h_bkg_vec[i]->Integral();
      h_vec.push_back( h_bkg_vec[i] );
      h_vec_norm.push_back( (TH1D*) h_bkg_vec[i]->Clone(TString::Format("%s_NORM",h_bkg_vec[i]->GetName())) );
    }

    for( std::size_t i = 0 ; i < h_sig_vec.size() ; ++i ) {
      legend->AddEntry( h_sig_vec[i] , h_sig_vec[i]->GetName() , "l" );
      yields->AddEntry( (TH1D*)NULL , TString::Format("%0.1f",h_sig_vec[i]->Integral()) , "" );
      h_vec.push_back( h_sig_vec[i] );
      h_vec_norm.push_back( (TH1D*) h_sig_vec[i]->Clone(TString::Format("%s_NORM",h_sig_vec[i]->GetName())) );
    }

    if( kFALSE && h_bkg_vec.size() > 0 ) {
      legend->AddEntry( (TH1D*)NULL , "sum(BKG)" , "" );
      yields->AddEntry( (TH1D*)NULL , TString::Format("%0.1f",sum_bkg) , "" );
    }
    
    assert( h_vec.size() == h_vec_norm.size() );
    
    Double_t ymax ( 0 );
    Double_t ymax_norm ( 0 );
    for( std::size_t i = 0 ; i < h_vec.size() ; ++i ) {
      htools::normalize( h_vec_norm[i] );
      h_vec_norm[i]->GetYaxis()->SetTitle( "Normalized" );
      ymax = TMath::Max( ymax , h_vec[i]->GetMaximum() );
      ymax_norm = TMath::Max( ymax_norm , h_vec_norm[i]->GetMaximum() );
    }
    
    htools::setup( legend , htools::LEFT , htools::TOP , 0.8 );
    htools::setup( yields , htools::CENTER , htools::TOP , 0.8 );
    htools::setup( text , htools::RIGHT , htools::TOP , 0.8 );

    // Use consistent y-max for all histograms
    ymax *= 1.8;
    ymax_norm *= 1.8;
    for( std::size_t i = 0 ; i < h_vec.size() ; ++i ) h_vec[i]->SetMaximum( ymax );
    for( std::size_t i = 0 ; i < h_vec_norm.size() ; ++i ) h_vec_norm[i]->SetMaximum( ymax_norm );

    // Possibly use log scaling
    Bool_t logy = kFALSE;
    if( hconfig.xname.Contains("eta") || hconfig.xname.Contains("phi") || hconfig.xname==TString("mu") ) {
      // keep linear scale
    } else {
      logy = kTRUE;
      for( std::size_t i = 0 ; i < h_vec.size() ; ++i ) { h_vec[i]->SetMinimum(0.1); h_vec[i]->SetMaximum( 100 * ymax ); }
      for( std::size_t i = 0 ; i < h_vec_norm.size() ; ++i ) { h_vec_norm[i]->SetMinimum(0.1); h_vec_norm[i]->SetMaximum( 100 * ymax_norm ); }
    }
    
    
    
    //
    // Create histograms ratio w.r.t. leading background (or last signal if there is no background)
    //
    
    TH1D *h_uncert	  = (TH1D*)0;
    TH1D *h_nouncert	  = (TH1D*)0;
    TH1D *h_uncert_norm	  = (TH1D*)0;
    TH1D *h_nouncert_norm = (TH1D*)0;
    std::vector<TH1D*> h_vec_ratio;
    std::vector<TH1D*> h_vec_norm_ratio;
    TLine *line	       = (TLine*)0;

    if( h_vec.size() > 0 ) {

      h_uncert = (TH1D*) h_vec.front()->Clone("h_uncert");
      h_uncert_norm = (TH1D*) h_vec_norm.front()->Clone("h_uncert_norm");
      
      h_nouncert = (TH1D*) h_vec.front()->Clone("h_nouncert");
      h_nouncert_norm = (TH1D*) h_vec_norm.front()->Clone("h_nouncert_norm");
      for( Int_t ibin = 1 ; ibin <= h_nouncert->GetNbinsX() ; ibin++ ) {
	h_nouncert->SetBinError(ibin,0.0);
	h_nouncert_norm->SetBinError(ibin,0.0);
      }

      h_vec_ratio.push_back( (TH1D*) h_uncert->Clone("h_ratio_0") );
      h_vec_ratio.back()->Divide( h_nouncert );
      h_vec_ratio.back()->SetFillColor( kYellow - 7 );
      h_vec_ratio.back()->SetLineColor( kBlack );
      h_vec_ratio.back()->SetFillStyle( 1001 );
      h_vec_ratio.back()->SetMinimum( 0 );
      h_vec_ratio.back()->SetMaximum( 3 );

      h_vec_norm_ratio.push_back( (TH1D*) h_uncert_norm->Clone("h_ratio_norm_0") );
      h_vec_norm_ratio.back()->Divide( h_nouncert_norm );
      h_vec_norm_ratio.back()->SetFillColor( kYellow - 7 );
      h_vec_norm_ratio.back()->SetLineColor( kBlack );
      h_vec_norm_ratio.back()->SetFillStyle( 1001 );
      h_vec_norm_ratio.back()->SetMinimum( 0 );
      h_vec_norm_ratio.back()->SetMaximum( 3 );
      
      for( std::size_t i = 1 ; i < h_vec.size() ; ++i ) {
	h_vec_ratio.push_back( (TH1D*) h_vec[i]->Clone(TString::Format("h_ratio_%i",Int_t(i))) );
	h_vec_norm_ratio.push_back( (TH1D*) h_vec_norm[i]->Clone(TString::Format("h_ratio_norm_%i",Int_t(i))) );
	h_vec_ratio.back()->Divide( h_nouncert );
	h_vec_norm_ratio.back()->Divide( h_nouncert_norm );
      }
      
      line = new TLine( hconfig.xmin , 1 , hconfig.xmax , 1 );
      line->SetLineWidth( 2 );
      line->SetLineStyle( kDashed );
      line->SetLineColor( h_uncert->GetLineColor() );

      h_vec[0]->GetXaxis()->SetTitle( "" );
      h_vec[0]->GetXaxis()->SetLabelSize( 0 );
      h_vec[0]->GetXaxis()->SetLabelOffset( 999 );

      h_vec_norm[0]->GetXaxis()->SetTitle( "" );
      h_vec_norm[0]->GetXaxis()->SetLabelSize( 0 );
      h_vec_norm[0]->GetXaxis()->SetLabelOffset( 999 );
      
      h_vec_ratio[0]->GetYaxis()->SetNdivisions( 5 , 5 , 0 , kTRUE );
      h_vec_ratio[0]->GetYaxis()->SetTitle( "Ratio" );

      h_vec_norm_ratio[0]->GetYaxis()->SetNdivisions( 5 , 5 , 0 , kTRUE );
      h_vec_norm_ratio[0]->GetYaxis()->SetTitle( "Ratio" );
      
      // Draw the plots
      p->getCanvas()->Clear();
      TPad *padtop = new TPad( "padtop" , "" , 0 , 0 , 1 , 1 );
      TPad *padbottom = new TPad( "padbottom" , "" , 0 , 0 , 1 , 1 );
      htools::setupTop( padtop );
      htools::setupBottom( padbottom );
      padbottom->SetGridy();
      if( logy ) padtop->SetLogy();

      padtop->Draw();
      padtop->cd();
      h_vec[0]->Draw( "hist e" );
      for( std::size_t i = 1 ; i < h_vec.size() ; ++i ) h_vec[i]->Draw( "hist e same" );
      legend->Draw();
      yields->Draw();
      text->Draw();

      p->getCanvas()->cd();
      padbottom->Draw();
      padbottom->cd();
      h_vec_ratio[0]->Draw( "e2" );
      line->Draw();
      for( std::size_t i = 1 ; i < h_vec_ratio.size() ; ++i ) h_vec_ratio[i]->Draw( "pe0 same" );
      p->print();

      padtop->cd();
      h_vec_norm[0]->Draw( "hist e" );
      for( std::size_t i = 1 ; i < h_vec.size() ; ++i ) h_vec_norm[i]->Draw( "hist e same" );
      legend->Draw();
      yields->Draw();
      text->Draw();

      padbottom->cd();
      h_vec_norm_ratio[0]->Draw( "e2" );
      line->Draw();
      for( std::size_t i = 1 ; i < h_vec_ratio.size() ; ++i ) h_vec_norm_ratio[i]->Draw( "pe0 same" );
      p->print();
      
      p->getCanvas()->cd();
      p->getCanvas()->Clear( "D" );
      
      delete padtop;
      delete padbottom;
      
    }

    //
    // Write histograms to output ROOT file
    //

    htools::safe_delete( h_vec );
    htools::safe_delete( h_vec_norm );
    htools::safe_delete( h_vec_ratio );
    htools::safe_delete( h_vec_norm_ratio );
    htools::safe_delete( h_uncert );
    htools::safe_delete( h_nouncert );
    htools::safe_delete( h_uncert_norm );
    htools::safe_delete( h_nouncert_norm );
    delete line;
    delete legend;
    delete yields;
    delete text;

  }

  
  
};

#endif


  
  
